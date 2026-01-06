#!/usr/bin/env bash
set -euo pipefail

DATASET="${1:-}"
shift || true

BASE="/sci/labs/orzuk/orzuk/projects/SingleCell/Data"
DO_EXTRACT=1
TS_MAX_FILES=0   # 0 = all files; set >0 to download only first N (debug)

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base) BASE="$2"; shift 2 ;;
    --no-extract) DO_EXTRACT=0; shift 1 ;;
    --ts-max-files) TS_MAX_FILES="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: $0 <dataset> [--base DIR] [--no-extract] [--ts-max-files N]"
      echo "Datasets: humanImmuneAging, ratCR, tabulaSapiens_v2"
      exit 0 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$DATASET" ]]; then
  echo "ERROR: missing dataset" >&2
  exit 1
fi

RAW_DIR="${BASE}/${DATASET}/raw"
META_DIR="${BASE}/${DATASET}/metadata"
LOG_DIR="${BASE}/${DATASET}/logs"
mkdir -p "$RAW_DIR" "$META_DIR" "$LOG_DIR"

LOG="${LOG_DIR}/download_$(date +%Y%m%d_%H%M%S).log"
touch "$LOG"
log(){ echo "[$(date +'%F %T')] $*" | tee -a "$LOG" >&2; }

need(){ command -v "$1" >/dev/null 2>&1 || { log "FATAL missing: $1"; exit 2; }; }
need wget
need tar
need gzip
need file
need python3

HAVE_UNZIP=0
command -v unzip >/dev/null 2>&1 && HAVE_UNZIP=1

WGET_COMMON=( --continue --tries=50 --timeout=60 --wait=1 --retry-connrefused --read-timeout=60 --no-verbose --progress=dot:giga )

download_url_to(){
  local url="$1"
  local out="$2"
  log "Downloading: $url"
  log " -> $out"
  mkdir -p "$(dirname "$out")"
  wget "${WGET_COMMON[@]}" -O "$out" "$url" 2>&1 | tee -a "$LOG"
  [[ -s "$out" ]] || { log "FATAL: empty download: $out"; exit 3; }
}

extract_any(){
  local path="$1"; local dest="$2"
  [[ $DO_EXTRACT -eq 1 ]] || return 0
  mkdir -p "$dest"
  local ft; ft="$(file -b "$path" || true)"
  log "Extract check: $(basename "$path") | $ft"
  if echo "$ft" | grep -qiE 'tar archive|gzip compressed data'; then
    tar -xf "$path" -C "$dest" 2>&1 | tee -a "$LOG" || true
  fi
  if echo "$ft" | grep -qiE 'Zip archive data'; then
    if [[ $HAVE_UNZIP -eq 1 ]]; then
      unzip -n "$path" -d "$dest" 2>&1 | tee -a "$LOG"
    else
      log "FATAL: unzip missing for zip file"
      exit 5
    fi
  fi
}

geo_scrape_suppl(){
  local gse="$1"; local regex="$2"
  local acc="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gse}"
  log "GEO scrape: $acc"
  local html="${RAW_DIR}/${gse}.acc.html"
  wget -q -O "$html" "$acc"
  local links="${RAW_DIR}/${gse}.links.txt"
  grep -oE 'https://ftp\.ncbi\.nlm\.nih\.gov/geo/series/[^"]+/suppl/[^"]+' "$html" \
    | sed 's/&amp;/\&/g' | sort -u > "$links" || true
  [[ -s "$links" ]] || { log "FATAL: no suppl links found"; exit 8; }
  local filtered="${RAW_DIR}/${gse}.filtered.txt"
  grep -E "$regex" "$links" > "$filtered" || true
  [[ -s "$filtered" ]] || { log "FATAL: no links match $regex"; exit 9; }
  log "GEO: downloading $(wc -l < "$filtered") files"
  while IFS= read -r url; do
    [[ -z "$url" ]] && continue
    local out="${RAW_DIR}/$(basename "$url")"
    download_url_to "$url" "$out"
    extract_any "$out" "$RAW_DIR"
  done < "$filtered"
}

figshare_write_manifest(){
  local article_id="$1"
  local manifest="${META_DIR}/figshare_files_${article_id}.tsv"
  log "Figshare article=$article_id"
  log "Manifest -> $manifest"

  python3 - <<PY 2>&1 | tee -a "$LOG"
import json, urllib.request

article_id="${article_id}"
manifest="${manifest}"

def jget(url):
    req=urllib.request.Request(url, headers={"User-Agent":"download_aging_data/4.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.loads(r.read().decode("utf-8"))

meta=jget(f"https://api.figshare.com/v2/articles/{article_id}")
files=meta.get("files", [])
if not files:
    files=jget(f"https://api.figshare.com/v2/articles/{article_id}/files")

with open(manifest,"w") as w:
    w.write("name\tsize\tid\tdownload_url\n")
    for f in files:
        name=f.get("name","")
        size=str(f.get("size",""))
        fid=str(f.get("id",""))
        url=f.get("download_url") or (f"https://figshare.com/ndownloader/files/{fid}" if fid else "")
        w.write(f"{name}\t{size}\t{fid}\t{url}\n")

print("Wrote manifest with", len(files), "files:", manifest)
PY

  echo "$manifest"
}

figshare_download_from_manifest_with_wget(){
  local manifest="$1"
  log "Downloading Figshare files via wget from manifest: $manifest"

  # Skip header; columns: name, size, id, download_url
  local n=0
  tail -n +2 "$manifest" | while IFS=$'\t' read -r name size fid url; do
    [[ -z "$fid" ]] && continue
    [[ -z "$url" ]] && url="https://figshare.com/ndownloader/files/${fid}"

    # Create a safe filename
    local fname="$name"
    if [[ -z "$fname" ]]; then
      fname="file_${fid}"
    fi
    # Avoid weird slashes
    fname="${fname//\//_}"

    local out="${RAW_DIR}/${fname}"
    if [[ -s "$out" ]]; then
      log "exists: $fname"
    else
      download_url_to "$url" "$out"
      extract_any "$out" "$RAW_DIR"
    fi

    n=$((n+1))
    if [[ "$TS_MAX_FILES" -gt 0 && "$n" -ge "$TS_MAX_FILES" ]]; then
      log "Stopping early due to --ts-max-files $TS_MAX_FILES"
      break
    fi
  done
}

log "=== download_aging_data start ==="
log "dataset=$DATASET"
log "raw_dir=$RAW_DIR"
log "log=$LOG"

case "$DATASET" in
  humanImmuneAging)
    geo_scrape_suppl "GSE299043" '\.h5ad$'
    ;;
  ratCR)
    geo_scrape_suppl "GSE137869" 'filtered_feature_bc_matrix\.tar\.gz$|\.tar\.gz$|\.tgz$|\.mtx\.gz$|features\.tsv\.gz$|barcodes\.tsv\.gz$'
    ;;
  tabulaSapiens_v2)
    manifest="$(figshare_write_manifest "27921984")"
    figshare_download_from_manifest_with_wget "$manifest"
    ;;
  *)
    log "FATAL: unknown dataset $DATASET"
    exit 1
    ;;
esac

log "Top-level raw listing:"
ls -lh "$RAW_DIR" | head -n 60 | tee -a "$LOG" >/dev/null || true
log "=== download_aging_data DONE dataset=$DATASET ==="
