#!/usr/bin/env bash
# WSL 홈에서 실행한 결과 중 "작고 공유할 가치가 있는 것"만 Windows 레포의 docs/example_output/으로 복사.
# BAM·FASTQ·레퍼런스 같은 대용량 파일은 git에 올리지 않는다(재실행하면 다시 만들어짐).
# 사용(WSL 안에서): bash /mnt/c/Users/<사용자>/ngs-variant-calling-pipeline/scripts/export_results.sh
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN="${1:-$HOME/ngs}"
OUT="$REPO/docs/example_output"
mkdir -p "$OUT"
cp "$RUN"/results/report/* "$OUT/"
cp "$RUN"/results/benchmark/*.tsv "$OUT/"
cp "$RUN"/results/qc/multiqc_report.html "$OUT/"
cp "$RUN"/results/qc/*.mosdepth.summary.txt "$RUN"/results/qc/*.flagstat.txt "$OUT/"
cp "$RUN"/results/vcf/*.final.vcf.gz "$RUN"/results/vcf/*.final.vcf.gz.tbi "$OUT/"
cp "$RUN"/resources/targets/panel.bed "$OUT/"
cp "$RUN"/workflow/envs/environment.lock.yaml "$REPO/workflow/envs/" 2>/dev/null || true
echo "[export] -> $OUT"
ls -la "$OUT"
