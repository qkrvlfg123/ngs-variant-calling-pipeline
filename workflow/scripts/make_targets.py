"""GENCODE GTF에서 패널 유전자의 좌표를 뽑아 타깃 BED 파일을 만든다.

필요한 이유:
  - 임상 패널 검사는 분석 대상 영역을 BED 파일로 정의한다.
  - 좌표를 손으로 적으면 오타·빌드(GRCh37/38) 혼동 위험이 있으므로,
    공식 유전자 주석(GENCODE)에서 유전자 이름으로 자동 추출한다.

BED 형식: chrom, start(0-based), end, name
"""
import sys

gtf_path, bed_path = snakemake.input.gtf, snakemake.output.bed
genes = set(snakemake.config["genes"])
pad = int(snakemake.config["padding"])

rows = []
with open(gtf_path) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "gene":
            continue
        # attribute 컬럼에서 gene_name "XXX"; 추출
        attrs = dict(
            kv.strip().split(" ", 1) for kv in f[8].rstrip(";").split(";") if kv.strip()
        )
        name = attrs.get("gene_name", "").strip('"')
        if name in genes:
            # GTF는 1-based 닫힌 구간, BED는 0-based 반열린 구간 -> start에서 1을 뺀다
            start = max(0, int(f[3]) - 1 - pad)
            end = int(f[4]) + pad
            rows.append((f[0], start, end, name))

found = {r[3] for r in rows}
missing = genes - found
if missing:
    sys.exit(f"GTF에서 찾지 못한 유전자: {sorted(missing)}")

rows.sort(key=lambda r: (r[0], r[1]))
with open(bed_path, "w") as out:
    for chrom, start, end, name in rows:
        out.write(f"{chrom}\t{start}\t{end}\t{name}\n")
