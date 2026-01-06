#!/usr/bin/env bash
set -euo pipefail

# Unified downloader for aging datasets (scRNA + proteomics).
#
# Usage:
#   ./download_aging_data.sh <dataset> [--base DIR] [--no-extract] [--ts-max-files N] [--overwrite]
#   ./download_aging_data.sh ALL [--base DIR] ...   # submits sbatch jobs for all datasets
#
# Design goals:
# - Same CLI for all datasets.
# - Dataset-specific messy details live inside this script.
# - Robust logging and resumable downloads.

# List of all available datasets (used by ALL mode)
ALL_DATASETS=(
  humanImmuneAging
  ratCR
  tabulaSapiens_v2
  macaque_30tissues_PXD066108
  mouse_41organs_8tp_PR_PXD053154
  mouse_10organs_4to20mo_PXD047296
  rat_brain_liver_transcriptome_proteome_PXD002467
  mouse_aging_lung_scRNA_proteome_PXD012307
  human_plasma_age_PXD016199
  human_plasma_age_pred_PXD028281
  human_skin_young_old_PXD018430
  mouse_8organs_lifestages_PXD058684_jpost
)

# Handle --help early, before any directory operations
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  cat <<'EOF'
Usage: download_aging_data.sh <dataset> [OPTIONS]
       download_aging_data.sh ALL [OPTIONS]   # submit sbatch jobs for all datasets

Options:
  --base DIR        Output directory (default: /sci/labs/orzuk/orzuk/projects/SingleCell/Data)
  --no-extract      Skip automatic extraction of archives
  --ts-max-files N  Download only first N files (for testing large datasets)
  --overwrite       Re-download files even if they exist (default: skip existing)

Core scRNA datasets:
  humanImmuneAging
  ratCR
  tabulaSapiens_v2

Proteomics aging datasets (PRIDE / ProteomeXchange unless noted):
  macaque_30tissues_PXD066108
  mouse_41organs_8tp_PR_PXD053154
  mouse_10organs_4to20mo_PXD047296
  rat_brain_liver_transcriptome_proteome_PXD002467
  mouse_aging_lung_scRNA_proteome_PXD012307
  human_plasma_age_PXD016199
  human_plasma_age_pred_PXD028281
  human_skin_young_old_PXD018430
  mouse_8organs_lifestages_PXD058684_jpost (jPOST entry JPST003472; very large)

Special:
  ALL               Submit sbatch jobs for all datasets (runs in parallel on cluster)

Notes:
- Many proteomics datasets are huge. Use --ts-max-files N to sanity check the downloader first.
- By default, existing files are skipped. Use --overwrite to re-download.
EOF
  exit 0
fi

DATASET="${1:-}"
shift || true

BASE="/sci/labs/orzuk/orzuk/projects/SingleCell/Data"
DO_EXTRACT=1
TS_MAX_FILES=0   # 0 = all files; set >0 to download only first N (debug/safety)
OVERWRITE=0      # 0 = skip existing files; 1 = re-download everything

# Collect remaining args for forwarding to sbatch in ALL mode
FORWARD_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base) BASE="$2"; FORWARD_ARGS+=("--base" "$2"); shift 2 ;;
    --no-extract) DO_EXTRACT=0; FORWARD_ARGS+=("--no-extract"); shift 1 ;;
    --ts-max-files) TS_MAX_FILES="$2"; FORWARD_ARGS+=("--ts-max-files" "$2"); shift 2 ;;
    --overwrite) OVERWRITE=1; FORWARD_ARGS+=("--overwrite"); shift 1 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$DATASET" ]]; then
  echo "ERROR: missing dataset" >&2
  exit 1
fi

# ---------- Handle ALL mode: submit sbatch jobs for each dataset ----------
if [[ "$DATASET" == "ALL" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  SBATCH_SCRIPT="${SCRIPT_DIR}/download_aging_data.sbatch"

  if [[ ! -x "$SBATCH_SCRIPT" ]]; then
    echo "ERROR: sbatch script not found: $SBATCH_SCRIPT" >&2
    exit 1
  fi

  echo "Submitting sbatch jobs for ${#ALL_DATASETS[@]} datasets..."
  for ds in "${ALL_DATASETS[@]}"; do
    echo "  sbatch $ds"
    sbatch "$SBATCH_SCRIPT" "$ds" "${FORWARD_ARGS[@]}"
  done
  echo "Done. Use 'squeue -u \$USER' to monitor jobs."
  exit 0
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
HAVE_LFTP=0
command -v lftp >/dev/null 2>&1 && HAVE_LFTP=1

# wget tuned for unreliable networks + large files
WGET_COMMON=( --continue --tries=50 --timeout=60 --wait=1 --retry-connrefused --read-timeout=60 --no-verbose --progress=dot:giga )


download_url_to(){
  local url="$1"
  local out="$2"

  # Skip if file exists and overwrite is disabled
  if [[ $OVERWRITE -eq 0 && -s "$out" ]]; then
    log "SKIP (exists): $(basename "$out")"
    return 0
  fi

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
  grep -oE '(https?|ftp)://ftp\.ncbi\.nlm\.nih\.gov/geo/series/[^"]+/suppl/[^"]+' "$html" \
    | sed 's/&amp;/\&/g' \
    | python3 -c 'import sys, urllib.parse; [print(urllib.parse.unquote(l.strip())) for l in sys.stdin]' \
    | sort -u > "$links" || true
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
    req=urllib.request.Request(url, headers={"User-Agent":"download_aging_data/5.0"})
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

download_from_manifest_tsv(){
  # TSV columns: name size id download_url
  local manifest="$1"
  log "Downloading files via wget from manifest: $manifest"
  local n=0
  tail -n +2 "$manifest" | while IFS=$'\t' read -r name size fid url; do
    [[ -z "$fid" && -z "$url" ]] && continue
    [[ -z "$url" ]] && url="https://figshare.com/ndownloader/files/${fid}"

    local fname="$name"
    if [[ -z "$fname" ]]; then
      fname="file_${fid:-unknown}"
    fi
    fname="${fname//\//_}"

    local out="${RAW_DIR}/${fname}"
    if [[ $OVERWRITE -eq 0 && -s "$out" ]]; then
      log "SKIP (exists): $fname"
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

# ---------- PRIDE / ProteomeXchange helpers ----------

pride_mirror_ftp_dir(){
  # Mirrors an FTP directory into RAW_DIR (or a subdir). Works with very large datasets.
  # Uses lftp mirror if available, else wget recursive.
  local ftp_dir="$1"     # e.g. ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2025/11/PXD066108
  local dest="${2:-$RAW_DIR}"

  mkdir -p "$dest"
  log "PRIDE mirror:"
  log "  src=$ftp_dir"
  log "  dest=$dest"
  log "  tool=$([[ $HAVE_LFTP -eq 1 ]] && echo lftp || echo wget)"
  log "  overwrite=$([[ $OVERWRITE -eq 1 ]] && echo yes || echo no)"

  if [[ $HAVE_LFTP -eq 1 ]]; then
    if [[ "$TS_MAX_FILES" -gt 0 ]]; then
      log "WARN: --ts-max-files is not enforced for lftp mirror (depends on lftp version)."
    fi
    # --only-newer: skip files that exist and are same size (unless --overwrite)
    local lftp_mirror_opts="--parallel=4"
    [[ $OVERWRITE -eq 0 ]] && lftp_mirror_opts+=" --only-newer"
    # Disable SSL for PRIDE FTP (plain FTP)
    lftp -e "set ftp:ssl-allow no; set net:max-retries 50; set net:timeout 60; mirror ${lftp_mirror_opts} ${ftp_dir#ftp://ftp.pride.ebi.ac.uk} \"$dest\"; bye" ftp.pride.ebi.ac.uk \
      2>&1 | tee -a "$LOG"
  else
    local args=( -r -np -nH --cut-dirs=4 -R "index.html*" )
    # -nc (no-clobber): skip existing files unless --overwrite
    [[ $OVERWRITE -eq 0 ]] && args+=( -nc )
    wget "${WGET_COMMON[@]}" "${args[@]}" -P "$dest" "$ftp_dir/" 2>&1 | tee -a "$LOG"
  fi
}

pride_download_pxd(){
  local pxd="$1"
  local ym="$2" # yyyy/mm
  local ftp_dir="ftp://ftp.pride.ebi.ac.uk/pride/data/archive/${ym}/${pxd}"
  pride_mirror_ftp_dir "$ftp_dir" "$RAW_DIR"
}

# ---------- jPOST helper (HTML scrape -> manifest) ----------

jpost_write_manifest(){
  # Scrape jPOST entry page for file download URLs.
  # This is intentionally heuristic; jPOST may change HTML.
  local jpst="$1"  # e.g. JPST003472
  local manifest="${META_DIR}/jpost_${jpst}_files.tsv"
  local url="https://repository.jpostdb.org/entry/${jpst}"
  log "jPOST entry: $url"
  log "Manifest -> $manifest"

  python3 - <<PY 2>&1 | tee -a "$LOG"
import re, urllib.request
jpst="${jpst}"
url="${url}"
manifest="${manifest}"

req=urllib.request.Request(url, headers={"User-Agent":"download_aging_data/5.0"})
html=urllib.request.urlopen(req, timeout=60).read().decode("utf-8", errors="ignore")

candidates=set()

# Absolute URLs
for m in re.findall(r'https?://[^\s"\']+', html):
    if 'jpostdb.org' in m and ('download' in m or 'file' in m):
        candidates.add(m.split('#')[0].rstrip(')'))

# Relative endpoints that look like downloads
for m in re.findall(r'/(?:download|download_file|downloadFile|file|files)[^"\'\s<>]+', html):
    if 'entry' in m:
        continue
    candidates.add("https://repository.jpostdb.org" + m.split('#')[0])

# hrefs
for m in re.findall(r'href="([^"]+)"', html):
    if ('download' in m or 'file' in m) and ('jpostdb' in m or m.startswith('/')):
        if m.startswith('/'):
            m="https://repository.jpostdb.org"+m
        candidates.add(m.split('#')[0])

bad_substrings=('twitter.com','facebook.com','mailto:','/entry/','/help','/about')
candidates=[c for c in candidates if not any(b in c for b in bad_substrings)]

with open(manifest,'w') as w:
    w.write("name\turl\n")
    for i,u in enumerate(sorted(set(candidates))):
        w.write(f"file_{i}\t{u}\n")

print(f"Wrote {len(candidates)} candidate links -> {manifest}")
if len(candidates)==0:
    print("WARN: jPOST scraping found 0 links. You may need to download manually from the web UI.")
PY

  echo "$manifest"
}

jpost_download_from_manifest(){
  local manifest="$1"
  log "Downloading jPOST candidate links from: $manifest"
  local n=0
  tail -n +2 "$manifest" | while IFS=$'\t' read -r name url; do
    [[ -z "$url" ]] && continue
    local out="${RAW_DIR}/${name}"
    if [[ $OVERWRITE -eq 0 && -s "$out" ]]; then
      log "SKIP (exists): $name"
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
log "overwrite=$([[ $OVERWRITE -eq 1 ]] && echo yes || echo no)"
log "log=$LOG"

case "$DATASET" in
  # ---- existing scRNA datasets ----
  humanImmuneAging)
    geo_scrape_suppl "GSE299043" '\.h5ad$'
    ;;
  ratCR)
    geo_scrape_suppl "GSE137869" 'filtered_feature_bc_matrix\.tar\.gz$|\.tar\.gz$|\.tgz$|\.mtx\.gz$|features\.tsv\.gz$|barcodes\.tsv\.gz$'
    ;;
  tabulaSapiens_v2)
    manifest="$(figshare_write_manifest "27921984")"
    download_from_manifest_tsv "$manifest"
    ;;

  # ---- proteomics: PRIDE / ProteomeXchange ----
  macaque_30tissues_PXD066108)
    pride_download_pxd "PXD066108" "2025/11"
    ;;
  mouse_41organs_8tp_PR_PXD053154)
    pride_download_pxd "PXD053154" "2025/10"
    ;;
  mouse_10organs_4to20mo_PXD047296)
    pride_download_pxd "PXD047296" "2025/08"
    ;;
  rat_brain_liver_transcriptome_proteome_PXD002467)
    pride_download_pxd "PXD002467" "2015/10"
    ;;
  mouse_aging_lung_scRNA_proteome_PXD012307)
    pride_download_pxd "PXD012307" "2019/01"
    ;;
  human_plasma_age_PXD016199)
    pride_download_pxd "PXD016199" "2020/10"
    ;;
  human_plasma_age_pred_PXD028281)
    pride_download_pxd "PXD028281" "2023/07"
    ;;
  human_skin_young_old_PXD018430)
    pride_download_pxd "PXD018430" "2020/08"
    ;;

  # ---- proteomics: jPOST (very large) ----
  mouse_8organs_lifestages_PXD058684_jpost)
    # PXD058684 is hosted on jPOST as JPST003472 (public; multi-organ aging life stages).
    # This can be multi-terabyte; use --ts-max-files for testing.
    manifest="$(jpost_write_manifest "JPST003472")"
    jpost_download_from_manifest "$manifest"
    ;;

  *)
    log "FATAL: unknown dataset $DATASET"
    exit 1
    ;;
esac

log "Top-level raw listing:"
ls -lh "$RAW_DIR" | head -n 80 | tee -a "$LOG" >/dev/null || true
log "=== download_aging_data DONE dataset=$DATASET ==="
