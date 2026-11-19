"""GIAB truth 대비 변이검출 정확도(precision / recall / F1)를 계산한다.

입력: bcftools isec 결과 폴더 (truth vs 내 파이프라인 결과)
  0000.vcf : truth에만 있는 변이      -> 놓친 변이 (False Negative)
  0001.vcf : 내 결과에만 있는 변이    -> 잘못 검출한 변이 (False Positive)
  0002.vcf : 양쪽 모두 (truth 쪽 레코드)  -> True Positive
  0003.vcf : 양쪽 모두 (내 결과 쪽 레코드)

지표의 의미:
  recall(민감도)    = TP / (TP + FN) : 진짜 변이 중 내가 찾아낸 비율
  precision(정밀도) = TP / (TP + FP) : 내가 보고한 변이 중 진짜인 비율
  genotype 일치율   : 위치·염기가 맞은 변이 중 접합성(0/1 vs 1/1)까지 맞은 비율

한계: 위치·REF·ALT가 정확히 같은 경우만 일치로 본다(allele 단위 비교).
      hap.py 같은 haplotype 기반 비교 도구와 달리, 같은 변이를 다르게 표기한
      복잡한 indel은 불일치로 셀 수 있다 -> 양쪽 모두 bcftools norm으로 표기를 통일해 완화.
"""
import os

isec_dir = snakemake.input.isec
out_tsv = snakemake.output.tsv
out_miss = snakemake.output.discordant


def read_vcf(path):
    """VCF 본문을 (chrom, pos, ref, alt, genotype, filter) 목록으로 읽는다."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            fmt = f[8].split(":")
            gt = f[9].split(":")[fmt.index("GT")]
            # phased(0|1)와 unphased(0/1)를 같은 것으로 취급, 대립유전자 순서도 무시
            gt = "/".join(sorted(gt.replace("|", "/").split("/")))
            rows.append((f[0], int(f[1]), f[3], f[4], gt))
    return rows


def vtype(ref, alt):
    return "SNP" if len(ref) == 1 and len(alt) == 1 else "INDEL"


fn = read_vcf(os.path.join(isec_dir, "0000.vcf"))
fp = read_vcf(os.path.join(isec_dir, "0001.vcf"))
tp_truth = read_vcf(os.path.join(isec_dir, "0002.vcf"))
tp_query = read_vcf(os.path.join(isec_dir, "0003.vcf"))

with open(out_tsv, "w") as out:
    out.write("type\tTP\tFP\tFN\tprecision\trecall\tF1\tgenotype_concordance\n")
    for t in ("SNP", "INDEL", "ALL"):
        keep = lambda r: t == "ALL" or vtype(r[2], r[3]) == t
        n_tp = sum(1 for r in tp_truth if keep(r))
        n_fp = sum(1 for r in fp if keep(r))
        n_fn = sum(1 for r in fn if keep(r))
        gt_ok = sum(1 for a, b in zip(tp_truth, tp_query) if keep(a) and a[4] == b[4])
        prec = n_tp / (n_tp + n_fp) if n_tp + n_fp else float("nan")
        rec = n_tp / (n_tp + n_fn) if n_tp + n_fn else float("nan")
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else float("nan")
        conc = gt_ok / n_tp if n_tp else float("nan")
        out.write(f"{t}\t{n_tp}\t{n_fp}\t{n_fn}\t{prec:.4f}\t{rec:.4f}\t{f1:.4f}\t{conc:.4f}\n")

# 불일치 변이 목록: "왜 틀렸는가"를 IGV 등에서 하나씩 확인하기 위한 자료
with open(out_miss, "w") as out:
    out.write("class\tchrom\tpos\tref\talt\ttype\ttruth_GT\tquery_GT\n")
    for r in fn:
        out.write(f"FN\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{vtype(r[2], r[3])}\t{r[4]}\t.\n")
    for r in fp:
        out.write(f"FP\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{vtype(r[2], r[3])}\t.\t{r[4]}\n")
    for a, b in zip(tp_truth, tp_query):
        if a[4] != b[4]:
            out.write(f"GT_MISMATCH\t{a[0]}\t{a[1]}\t{a[2]}\t{a[3]}\t{vtype(a[2], a[3])}\t{a[4]}\t{b[4]}\n")
