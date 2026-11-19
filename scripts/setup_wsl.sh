#!/usr/bin/env bash
# ============================================================
# WSL(Ubuntu) 안에서 1회 실행하는 환경 구축 스크립트
#  - micromamba(가벼운 conda 구현체) 설치
#  - bioconda 채널에서 파이프라인 도구 전체를 'ngs' 환경으로 설치
# 사용: bash scripts/setup_wsl.sh
# ============================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
BIN_DIR="$HOME/.local/bin"
mkdir -p "$BIN_DIR"

# 1) micromamba 단일 바이너리 내려받기 (관리자 권한 불필요)
if [ ! -x "$BIN_DIR/micromamba" ]; then
    echo "[setup] micromamba 설치"
    curl -fsSL https://micro.mamba.pm/api/micromamba/linux-64/latest \
        | tar -xj -C "$HOME/.local" bin/micromamba
fi

# 2) 환경 생성 (이미 있으면 업데이트)
if "$BIN_DIR/micromamba" env list | grep -qE '^\s*ngs\s'; then
    echo "[setup] 기존 ngs 환경 업데이트"
    "$BIN_DIR/micromamba" install -y -n ngs -f "$REPO_DIR/workflow/envs/environment.yaml"
else
    echo "[setup] ngs 환경 생성"
    "$BIN_DIR/micromamba" create -y -n ngs -f "$REPO_DIR/workflow/envs/environment.yaml"
fi

# 3) 설치된 버전 기록 -> 재현성 근거 (README에 인용)
"$BIN_DIR/micromamba" env export -n ngs > "$REPO_DIR/workflow/envs/environment.lock.yaml"
echo "[setup] 완료. 실행: micromamba run -n ngs snakemake --cores 8"
