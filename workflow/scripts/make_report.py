"""최종 VCF를 사람이 읽을 수 있는 표와 요약 리포트(Markdown)로 바꾼다.

필요한 이유:
  VCF는 기계용 형식이라 한 줄에 정보가 빽빽하게 들어 있다. 실제 검사실에서도
  분석자는 VCF를 직접 읽지 않고, 유전자·HGVS 표기·접합성·ClinVar 분류로 정리된 표를 보고 판독한다.

출력:
  1) variants.tsv : PASS 변이 전체 (유전자, 영향, HGVS, genotype, depth, VAF, ClinVar)
  2) pgx_sites.tsv: 약물유전체 핵심 자리(star allele 정의 변이)의 genotype
     - 변이가 "없다"는 것도 결과다. 단, depth가 충분할 때만 "정상(ref)"이라고 말할 수 있으므로
       VCF에 없는 자리는 BAM에서 depth를 직접 확인한다.
  3) report.md    : QC·정확도·주요 변이 요약
"""
import gzip
import subprocess

inp, out, cfg = snakemake.input, snakemake.output, snakemake.config
SAMPLE = cfg["sample"]
MIN_DEPTH = int(cfg["report"]["min_depth_for_ref_call"])


def read_vcf(path):
    """VCF -> dict 목록. INFO는 dict로, FORMAT은 첫 번째 샘플 기준 dict로 풀어 준다."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(kv.split("=", 1) if "=" in kv else (kv, "") for kv in f[7].split(";"))
            fmt = dict(zip(f[8].split(":"), f[9].split(":"))) if len(f) > 9 else {}
            rows.append(
                dict(chrom=f[0], pos=int(f[1]), ref=f[3], alt=f[4], filter=f[6], info=info, fmt=fmt)
            )
    return rows


def first_ann(info):
    """SnpEff ANN의 첫 항목(영향이 가장 큰 것부터 정렬됨)을 꺼낸다."""
    ann = info.get("ANN", "")
    a = ann.split(",")[0].split("|") if ann else []
    a += [""] * (11 - len(a))
    return dict(effect=a[1], impact=a[2], gene=a[3], transcript=a[6], hgvs_c=a[9], hgvs_p=a[10])


def zygosity(gt):
    alleles = gt.replace("|", "/").split("/")
    if "." in alleles:
        return "판정불가"
    return "hom" if len(set(alleles)) == 1 else "het"


def vaf(fmt):
    """VAF(variant allele fraction) = ALT read 수 / 전체 read 수. het이면 ~0.5, hom이면 ~1.0 기대."""
    try:
        ad = [int(x) for x in fmt["AD"].split(",")]
        return f"{ad[1] / sum(ad):.2f}" if sum(ad) else "."
    except (KeyError, ValueError, IndexError):
        return "."


final = read_vcf(inp.final)
annotated = read_vcf(inp.annotated)  # 필터 탈락 변이 포함

# ---------------------------------------------------------------- 1) 변이 표
with open(out.variants, "w") as fh:
    cols = ["chrom", "pos", "ref", "alt", "gene", "effect", "impact", "HGVS_c", "HGVS_p",
            "genotype", "zygosity", "depth", "VAF", "rsID", "ClinVar_significance",
            "ClinVar_review_status", "ClinVar_disease"]
    fh.write("\t".join(cols) + "\n")
    for v in final:
        a, i, g = first_ann(v["info"]), v["info"], v["fmt"]
        rs = "rs" + i["RS"] if i.get("RS") else "."
        fh.write("\t".join(str(x) for x in [
            v["chrom"], v["pos"], v["ref"], v["alt"], a["gene"] or ".", a["effect"] or ".",
            a["impact"] or ".", a["hgvs_c"] or ".", a["hgvs_p"] or ".",
            g.get("GT", "."), zygosity(g.get("GT", "./.")), g.get("DP", "."), vaf(g), rs,
            i.get("CLNSIG", "."), i.get("CLNREVSTAT", "."), i.get("CLNDN", "."),
        ]) + "\n")

# ---------------------------------------------------------------- 2) PGx 자리
# 좌표를 손으로 적지 않고 ClinVar에서 rsID로 찾는다(빌드 혼동·오타 방지).
clinvar_by_rs = {}
with gzip.open(inp.clinvar, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 8)
        for kv in f[7].split(";"):
            if kv.startswith("RS="):
                for rs in kv[3:].split("|"):
                    clinvar_by_rs.setdefault(rs, []).append((f[0], int(f[1]), f[3], f[4]))


def depth_at(chrom, pos):
    r = subprocess.run(["samtools", "depth", "-a", "-r", f"{chrom}:{pos}-{pos}", inp.bam],
                       capture_output=True, text=True, check=True)
    return int(r.stdout.split()[2]) if r.stdout.strip() else 0


# 독립 검증: GIAB truth(여러 시퀀싱 기술로 교차 확인된 HG001의 정답)에서 같은 자리를 조회한다.
#  - truth에 변이가 있으면 그 genotype
#  - truth에 변이가 없고 고신뢰 영역 안이면 "정답도 변이 없음(0/0)"이라고 말할 수 있다
#  - 고신뢰 영역 밖이면 정답 자체가 불확실 -> 검증 불가
truth = {(v["chrom"], v["pos"], v["ref"], v["alt"]): v for v in read_vcf(inp.truth)}
confident = []
with open(inp.eval_bed) as fh:
    for line in fh:
        c, s, e = line.split()[:3]
        confident.append((c, int(s), int(e)))


def truth_genotype(chrom, pos, ref, alt):
    v = truth.get((chrom, pos, ref, alt))
    if v:
        gt = "/".join(sorted(v["fmt"]["GT"].replace("|", "/").split("/")))
        return f"{gt} ({zygosity(gt)})"
    # BED는 0-based 반열린 구간: start < pos(1-based) <= end
    if any(c == chrom and s < pos <= e for c, s, e in confident):
        return "0/0 (변이 없음)"
    return "고신뢰 영역 밖(검증 불가)"


pgx_rows = []
for site in cfg["report"]["pgx_sites"]:
    ref, alt = site["change"].split(">")
    hits = [c for c in clinvar_by_rs.get(str(site["rsid"]), []) if c[2] == ref and c[3] == alt]
    if not hits:
        raise SystemExit(f"ClinVar에서 rs{site['rsid']} {site['change']} 를 찾지 못함")
    chrom, pos = hits[0][0], hits[0][1]
    key = (chrom, pos, ref, alt)
    passed = [v for v in final if (v["chrom"], v["pos"], v["ref"], v["alt"]) == key]
    called = [v for v in annotated if (v["chrom"], v["pos"], v["ref"], v["alt"]) == key]
    dp = depth_at(chrom, pos)
    if passed:
        g = passed[0]["fmt"]
        result = f"{g['GT']} ({zygosity(g['GT'])}, AD={g.get('AD', '.')})"
    elif called:
        result = f"검출됐으나 필터 탈락({called[0]['filter']}) - 수동 확인 필요"
    elif dp >= MIN_DEPTH:
        result = "0/0 (변이 없음)"
    else:
        result = f"판정불가 (depth {dp} < {MIN_DEPTH})"
    pgx_rows.append([site["gene"], site["allele"], f"rs{site['rsid']}", f"{chrom}:{pos}",
                     site["change"], site["function"], str(dp), result,
                     truth_genotype(chrom, pos, ref, alt)])

with open(out.pgx, "w") as fh:
    fh.write("gene\tstar_allele\trsID\tposition\tchange\tallele_function\tdepth\tresult\tGIAB_truth\n")
    for r in pgx_rows:
        fh.write("\t".join(r) + "\n")

# ---------------------------------------------------------------- 3) Markdown 요약
def md_table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
    lines += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(lines)


# coverage: mosdepth thresholds = 유전자별로 10x/20x/30x/50x 이상인 염기 수
cov_rows = []
with gzip.open(inp.thresholds, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        c, s, e, name, x10, x20, x30, x50 = line.split()
        size = int(e) - int(s)
        cov_rows.append([name, f"{size:,}", f"{int(x20) / size:.1%}", f"{int(x30) / size:.1%}",
                         f"{int(x50) / size:.1%}"])

mean_depth = "."
with open(inp.summary) as fh:
    for line in fh:
        if line.startswith("total_region"):
            mean_depth = line.split()[3]

with open(inp.benchmark) as fh:
    bench = [line.rstrip("\n").split("\t") for line in fh]

impact_n = {}
for v in final:
    k = first_ann(v["info"])["impact"] or "주석없음"
    impact_n[k] = impact_n.get(k, 0) + 1

BENIGN = {"Benign", "Likely_benign", "Benign/Likely_benign"}
notable = []
for v in final:
    a, i, g = first_ann(v["info"]), v["info"], v["fmt"]
    sig = i.get("CLNSIG", "")
    # 단백질을 바꾸는 변이(HIGH/MODERATE)이거나, ClinVar 분류가 단순 benign이 아닌 변이
    if a["impact"] in ("HIGH", "MODERATE") or (sig and sig not in BENIGN):
        notable.append([a["gene"], f"{v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}",
                        a["hgvs_c"], a["hgvs_p"] or "-", a["effect"],
                        zygosity(g.get("GT", "./.")), vaf(g),
                        "rs" + i["RS"] if i.get("RS") else "-", sig.replace("|", ", ") or "-"])

n_filtered = sum(1 for v in annotated if v["filter"] != "PASS")
md = f"""# {SAMPLE} 변이검출 결과 요약

> 자동 생성 파일 (workflow/scripts/make_report.py). 연구·학습용이며 임상 판독 자료가 아님.

## 1. Coverage (타깃 평균 depth {mean_depth}x)

{md_table(["유전자(±padding)", "영역 크기(bp)", "≥20x", "≥30x", "≥50x"], cov_rows)}

## 2. 변이 수

- 검출 {len(annotated)}개 → hard filter 통과(PASS) {len(final)}개, 탈락 {n_filtered}개
- 영향도(SnpEff impact)별: {", ".join(f"{k} {n}" for k, n in sorted(impact_n.items()))}

## 3. GIAB v4.2.1 truth 대비 정확도 (패널 ∩ 고신뢰 영역)

{md_table(bench[0], bench[1:])}

## 4. 눈여겨볼 변이 (단백질 변화 또는 ClinVar 분류가 benign이 아닌 것)

{md_table(["유전자", "위치", "HGVS c.", "HGVS p.", "영향", "접합성", "VAF", "rsID", "ClinVar"], notable)}

## 5. 약물유전체(PGx) 핵심 자리

{md_table(["유전자", "star allele", "rsID", "위치(GRCh38)", "변화", "대립유전자 기능", "depth", "결과", "GIAB truth"], pgx_rows)}
"""
with open(out.md, "w", encoding="utf-8") as fh:
    fh.write(md)
