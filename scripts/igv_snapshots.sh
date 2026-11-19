#!/usr/bin/env bash
# IGV(Integrative Genomics Viewer) 배치 모드로 관심 변이 자리의 read 정렬 화면을 PNG로 저장한다.
# 이유: 변이검출기의 숫자(QD, MQ ...)만 보지 않고 "read가 실제로 어떻게 쌓여 있는지"를 눈으로 확인하는 것이
#     변이 판독의 마지막 단계다. 필터에서 탈락한 변이나 임상적으로 중요한 변이는 IGV로 직접 본다.
# 준비: sudo apt install xvfb ; micromamba create -n igv -c conda-forge -c bioconda igv
#       (IGV는 Java 21이 필요해 파이프라인 환경(ngs, Java 17)과 분리)
# 사용(WSL 안에서): bash scripts/igv_snapshots.sh [실행 폴더=~/ngs]
set -euo pipefail
RUN="${1:-$HOME/ngs}"
OUT="$RUN/results/igv"
mkdir -p "$OUT"

# 이름 | 중심 좌표 | 좌우 여유(bp)
SITES="
CYP2C19_star2_rs4244285|94781859|40
CYP2C9_star2_rs1799853|94942290|40
CYP2C19_filtered_lowMQ_94813990|94813990|40
CYP2C19_filtered_lowMQ_94817023_94817028|94817025|40
"

BATCH="$OUT/igv_batch.txt"
{
    echo "new"
    echo "genome $RUN/resources/ref/chr10.fa"
    echo "load $RUN/results/bam/HG001.markdup.bam"
    echo "load $RUN/results/vcf/HG001.annotated.vcf.gz"
    echo "load $RUN/resources/annot/chr10.mane.gtf"
    echo "snapshotDirectory $OUT"
    echo "maxPanelHeight 700"
    echo "colorBy NONE"
    for s in $SITES; do
        IFS='|' read -r name pos pad <<< "$s"
        echo "goto chr10:$((pos - pad))-$((pos + pad))"
        echo "sort base $pos"          # 해당 자리의 염기 기준으로 read 정렬 -> ALT read가 위로 모임
        echo "collapse"
        echo "snapshot $name.png"
    done
    echo "exit"
} > "$BATCH"

# xvfb = 화면 없는 서버에서 가상 디스플레이를 만들어 GUI 프로그램을 실행
xvfb-run --auto-servernum -s "-screen 0 1400x1000x24" \
    micromamba run -n igv igv -b "$BATCH"
ls -la "$OUT"
