#!/usr/bin/env bash
# Windows 쪽 레포(코드) -> WSL 홈(~/ngs)으로 코드만 복사한다.
# 이유: /mnt/c (특히 한글 경로)에서 직접 실행하면 I/O가 느리고 일부 도구가 경로를 못 읽는다.
#     코드는 Windows 레포가 원본, 실행·대용량 데이터는 WSL 홈에 둔다.
# 사용(WSL 안에서): bash /mnt/c/Users/<사용자>/ngs-variant-calling-pipeline/scripts/sync_to_wsl.sh
set -euo pipefail
SRC="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DST="${1:-$HOME/ngs}"
mkdir -p "$DST"
for d in config workflow scripts; do
    mkdir -p "$DST/$d"
    cp -r "$SRC/$d/." "$DST/$d/"
done
# Windows 줄바꿈(CRLF)이 섞이면 bash/Snakemake가 깨지므로 제거
find "$DST/config" "$DST/workflow" "$DST/scripts" -type f \
    \( -name '*.sh' -o -name '*.py' -o -name '*.yaml' -o -name 'Snakefile' \) \
    -exec sed -i 's/\r$//' {} +
echo "[sync] $SRC -> $DST"
