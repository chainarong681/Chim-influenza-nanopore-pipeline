#!/usr/bin/env bash
# =============================================================================
# flu_pipeline.sh - Influenza A (H5N1) Nanopore analysis, fully automatic
# Run inside Ubuntu (WSL2).  Needs flu_tools.py in the same folder.
#
#   bash flu_pipeline.sh -s                                  # one-time setup (install conda env + tools)
#   bash flu_pipeline.sh -i /mnt/c/Chim_nanopore/barcode21
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS="$SCRIPT_DIR/flu_tools.py"

IN_DIR="."
ENV_NAME="flu_analysis"
MODEL="r1041_e82_400bps_sup_v5.2.0"
THREADS="$(nproc)"
CAND=""
SEL_REF=""
# อีเมลสำหรับ NCBI Entrez (ตั้งค่าครั้งเดียวตรงนี้ ไม่ต้องใส่ -e ทุกครั้ง; ใส่ -e เพื่อ override ได้)
EMAIL="chim.lab@example.com"
NC_DATASET="community/moncla-lab/iav-h5/ha/all-clades"
MIN_LEN=200
MIN_Q=9
MIN_DEPTH=20
FORCE=0
SETUP_ONLY=0

usage() {
cat <<'EOF'
การใช้งาน: bash flu_pipeline.sh [ตัวเลือก]

  -i DIR    โฟลเดอร์ barcode ที่มี *.fastq.gz (รับทั้ง /mnt/c/... และ C:/...)   [ค่าเริ่มต้น: โฟลเดอร์ปัจจุบัน]
  -e EMAIL  อีเมลสำหรับ NCBI Entrez (ไม่บังคับ มีค่าเริ่มต้นในสคริปต์แล้ว)
  -c FILE   ใช้ candidate_refs.fasta ที่มีอยู่แล้ว (ข้ามการดาวน์โหลด)
  -R FILE   ระบุ reference ที่เลือกเอง (ข้ามการเลือกอัตโนมัติ) -> คัดลอกเป็น selected_reference.fasta
  -m MODEL  โมเดล Medaka                     [ค่าเริ่มต้น: r1041_e82_400bps_sup_v5.2.0]
  -t N      จำนวน threads                    [ค่าเริ่มต้น: ทุก core]
  -d NAME   ชื่อ Nextclade dataset           [ค่าเริ่มต้น: community/moncla-lab/iav-h5/ha/all-clades]
  -D N      depth ขั้นต่ำสำหรับรายงาน coverage [ค่าเริ่มต้น: 20]
  -f        รันใหม่ทุกขั้นตอน (ไม่ข้ามไฟล์ที่มีอยู่แล้ว)
  -s        ติดตั้ง Miniconda + environment + เครื่องมือทั้งหมด (ทำครั้งเดียว) แล้วจบ
  -h        แสดงวิธีใช้
EOF
}

while getopts "i:e:c:R:m:t:d:D:fsh" opt; do
  case "$opt" in
    i) IN_DIR="$OPTARG" ;;
    e) EMAIL="$OPTARG" ;;
    c) CAND="$OPTARG" ;;
    R) SEL_REF="$OPTARG" ;;
    m) MODEL="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    d) NC_DATASET="$OPTARG" ;;
    D) MIN_DEPTH="$OPTARG" ;;
    f) FORCE=1 ;;
    s) SETUP_ONLY=1 ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

log()  { echo -e "\n[$(date +%H:%M:%S)] $*"; }
warn() { echo "[WARNING] $*" >&2; }
die()  { echo "[ERROR] $*" >&2; exit 1; }
done_already() { [[ $FORCE -eq 0 && -s "$1" ]]; }
to_linux_path() {
  local p="${1//\\//}"
  if [[ "$p" =~ ^[A-Za-z]:/ ]]; then wslpath -u "$p"; else echo "$p"; fi
}

# ---------------------------------------------------------------- conda helpers
find_conda_sh() {
  local base
  if command -v conda >/dev/null 2>&1; then base="$(conda info --base)"; else base="$HOME/miniconda3"; fi
  echo "$base/etc/profile.d/conda.sh"
}

do_setup() {
  if [[ ! -f "$(find_conda_sh)" ]]; then
    log "ไม่พบ conda - กำลังติดตั้ง Miniconda"
    wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh
    bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
  fi
  set +u
  # shellcheck disable=SC1090
  source "$(find_conda_sh)"
  conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main || true
  conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r || true
  if ! conda env list | grep -q "^$ENV_NAME "; then
    log "สร้าง environment $ENV_NAME"
    conda create -y -n "$ENV_NAME" -c bioconda -c conda-forge minimap2 samtools medaka nanofilt ncbi-datasets-cli
  fi
  conda activate "$ENV_NAME"
  log "ติดตั้งเครื่องมือเพิ่มเติม"
  conda install -y -c bioconda -c conda-forge seqkit nanoplot nextclade bcftools gsl biopython "matplotlib<3.9" pandas
  set -u
  for t in minimap2 samtools medaka_consensus seqkit NanoPlot nextclade python; do
    command -v "$t" >/dev/null 2>&1 && echo "OK   $t" || warn "ไม่พบ $t"
  done
  log "ติดตั้งเสร็จ - รันวิเคราะห์ได้ด้วย: bash flu_pipeline.sh -i <โฟลเดอร์ barcode> -e <อีเมล>"
}

if [[ $SETUP_ONLY -eq 1 ]]; then do_setup; exit 0; fi

# ------------------------------------------------------------------ preparation
[[ -f "$TOOLS" ]] || die "ไม่พบ flu_tools.py ข้างไฟล์นี้ ($TOOLS)"
[[ -n "$CAND" ]]    && CAND="$(realpath -m "$(to_linux_path "$CAND")")"
[[ -n "$SEL_REF" ]] && SEL_REF="$(realpath -m "$(to_linux_path "$SEL_REF")")"
IN_DIR="$(to_linux_path "$IN_DIR")"
[[ -d "$IN_DIR" ]] || die "ไม่พบโฟลเดอร์ $IN_DIR"
cd "$IN_DIR"
CAND="${CAND:-$PWD/candidate_refs.fasta}"

exec > >(tee -a pipeline.log) 2>&1
log "เริ่มวิเคราะห์ที่: $PWD  (threads=$THREADS, model=$MODEL)"

CONDA_SH="$(find_conda_sh)"
[[ -f "$CONDA_SH" ]] || die "ไม่พบ conda - รัน: bash flu_pipeline.sh -s"
set +u
# shellcheck disable=SC1090
source "$CONDA_SH"
conda activate "$ENV_NAME" || die "ไม่พบ environment $ENV_NAME - รัน: bash flu_pipeline.sh -s"
set -u
for t in minimap2 samtools medaka_consensus seqkit NanoPlot nextclade python; do
  command -v "$t" >/dev/null 2>&1 || die "ไม่พบโปรแกรม $t - รัน: bash flu_pipeline.sh -s"
done

# ----------------------------------------------------------------------- steps
log "[1/11] รวมไฟล์ FASTQ"
if done_already all_raw_reads.fastq; then
  echo "ข้าม (มี all_raw_reads.fastq แล้ว)"
else
  mapfile -d '' FQ < <(find . -type f -name "*.fastq.gz" -print0)
  (( ${#FQ[@]} > 0 )) || die "ไม่พบไฟล์ *.fastq.gz ใน $PWD"
  echo "พบ ${#FQ[@]} ไฟล์"
  zcat "${FQ[@]}" > all_raw_reads.fastq
fi

log "[2/11] NanoPlot QC ก่อนกรอง"
if done_already nanoplot_raw_out/NanoPlot-report.html; then
  echo "ข้าม"
else
  NanoPlot --fastq all_raw_reads.fastq -o nanoplot_raw_out --plots hex -t "$THREADS" || warn "NanoPlot ล้มเหลว (ไปขั้นตอนต่อไป)"
fi

log "[3/11] กรอง read (>= ${MIN_LEN} bp, Q >= ${MIN_Q})"
if done_already filtered.fastq; then
  echo "ข้าม"
else
  seqkit seq -m "$MIN_LEN" -Q "$MIN_Q" all_raw_reads.fastq > filtered.fastq
fi
[[ -s filtered.fastq ]] || die "filtered.fastq ว่างเปล่า - ไม่มี read ผ่านเกณฑ์"
seqkit stats -a filtered.fastq | tee filtered_stats.txt
python "$TOOLS" readlen --fastq filtered.fastq --min-len "$MIN_LEN" --out read_length_qc.png

log "[4/11] เตรียม candidate references"
if [[ -s "$CAND" ]]; then
  echo "ใช้ $CAND ($(grep -c '^>' "$CAND") sequences)"
else
  python "$TOOLS" fetch-refs --email "$EMAIL" --out "$CAND"
fi

log "[5/11] Minimap2 screening กับ candidate refs"
if done_already full_check.sorted.bam.bai && done_already segment_distribution_stats.txt; then
  echo "ข้าม"
else
  # -F 0x904 = ตัด unmapped/secondary/supplementary ออก นับเฉพาะ primary alignment (1 read = 1 reference)
  minimap2 -ax map-ont --secondary=no -t "$THREADS" "$CAND" filtered.fastq \
    | samtools view -@ "$THREADS" -b -F 0x904 \
    | samtools sort -@ "$THREADS" -o full_check.sorted.bam
  samtools index full_check.sorted.bam
  samtools idxstats full_check.sorted.bam > segment_distribution_stats.txt
fi

log "[6/11] เลือก best reference ต่อเซกเมนต์"
if [[ -n "$SEL_REF" ]]; then
  cp "$SEL_REF" selected_reference.fasta
  echo "ใช้ reference ที่ระบุเอง: $SEL_REF"
elif done_already selected_reference.fasta; then
  echo "ข้าม (ใช้ selected_reference.fasta เดิม - แก้ไขด้วยมือได้แล้วรันซ้ำ)"
else
  python "$TOOLS" select-ref --idxstats segment_distribution_stats.txt --candidates "$CAND" --out-ref selected_reference.fasta
fi

log "[7/11] Alignment + coverage กับ selected_reference"
if done_already final_alignment.sorted.bam.bai && done_already coverage_depth.txt; then
  echo "ข้าม"
else
  minimap2 -ax map-ont --secondary=no -t "$THREADS" selected_reference.fasta filtered.fastq \
    | samtools sort -@ "$THREADS" -o final_alignment.sorted.bam
  samtools index final_alignment.sorted.bam
  samtools flagstat final_alignment.sorted.bam | tee flagstat.txt
  samtools depth -a final_alignment.sorted.bam > coverage_depth.txt
fi
python "$TOOLS" coverage --depth coverage_depth.txt --ref selected_reference.fasta --min-depth "$MIN_DEPTH" --out coverage_report.html

log "[8/11] Medaka consensus"
if done_already medaka_out/consensus.fasta; then
  echo "ข้าม"
else
  medaka_consensus -i filtered.fastq -d selected_reference.fasta -o medaka_out -m "$MODEL" -t "$THREADS" \
    || die "Medaka ล้มเหลว (ตรวจชื่อโมเดลด้วย -m หรือแพ็กเกจ bcftools/gsl: รัน -s ใหม่)"
fi
[[ -s medaka_out/consensus.fasta ]] || die "ไม่พบ medaka_out/consensus.fasta"

log "[9/11] QC consensus และโครงสร้าง"
python "$TOOLS" consensus-qc --ref selected_reference.fasta --consensus medaka_out/consensus.fasta \
  --depth coverage_depth.txt --min-depth "$MIN_DEPTH" --model "$MODEL" --out structural_qc_report.html

log "[10/11] Typing / Clade assignment (HA + Nextclade)"
python "$TOOLS" extract-ha --ref selected_reference.fasta --consensus medaka_out/consensus.fasta --out HA_segment_consensus.fasta \
  || warn "สกัด HA ไม่สำเร็จ ข้าม Nextclade"
if [[ -s HA_segment_consensus.fasta ]]; then
  if [[ ! -d h5_dataset || $FORCE -eq 1 ]]; then
    nextclade dataset get --name "$NC_DATASET" --output-dir h5_dataset || warn "ดาวน์โหลด Nextclade dataset ไม่สำเร็จ (ต้องใช้อินเทอร์เน็ต)"
  fi
  if [[ -d h5_dataset ]]; then
    nextclade run --input-dataset h5_dataset --output-all nextclade_out/ HA_segment_consensus.fasta \
      && python "$TOOLS" nextclade-html --tsv nextclade_out/nextclade.tsv --out nextclade_report.html \
      || warn "Nextclade ล้มเหลว"
  fi
fi

log "[11/11] ประเมินความรุนแรง (HA cleavage, PB2, NA)"
python "$TOOLS" virulence --ref selected_reference.fasta --consensus medaka_out/consensus.fasta --out virulence_report.html

log "เสร็จสมบูรณ์ ✅  ผลลัพธ์อยู่ที่ $PWD"
ls -1 ./*.html ./*.png ./*.txt medaka_out/consensus.fasta HA_segment_consensus.fasta 2>/dev/null || true
