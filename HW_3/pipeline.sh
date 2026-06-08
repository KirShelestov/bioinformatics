#!/usr/bin/env bash

set -euo pipefail

SAMPLE="SRR2584863"           # Accession из NCBI SRA (E.coli K-12)
THREADS=4
MAPPED_THRESHOLD=90           # Порог % картирования

REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz"
REF_GZ="ref/ecoli_k12.fa.gz"
REF="ref/ecoli_k12.fa"

FASTQ_1="data/${SAMPLE}_1.fastq.gz"
FASTQ_2="data/${SAMPLE}_2.fastq.gz"

SAM="results/${SAMPLE}.sam"
BAM="results/${SAMPLE}.bam"
BAM_SORTED="results/${SAMPLE}.sorted.bam"
FLAGSTAT="results/${SAMPLE}.flagstat.txt"
VCF="results/${SAMPLE}.vcf"
QC_DIR="results/fastqc_${SAMPLE}"
LOG="logs/pipeline_$(date +%Y%m%d_%H%M%S).log"

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
info()    { echo -e "${GREEN}[INFO]${NC}  $*" | tee -a "$LOG"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*" | tee -a "$LOG"; }
error()   { echo -e "${RED}[ERROR]${NC} $*" | tee -a "$LOG"; exit 1; }

check_deps() {
    info "Проверка зависимостей..."
    for cmd in fastqc bwa samtools freebayes prefetch fasterq-dump; do
        command -v "$cmd" &>/dev/null || error "Не найден: $cmd. Установи его и повтори."
    done
    info "Все зависимости найдены."
}

setup_dirs() {
    mkdir -p data ref results logs "$QC_DIR"
}

# --------------- Шаг 1: Скачивание ридов ------------------------------------
download_reads() {
    if [[ -f "$FASTQ_1" && -f "$FASTQ_2" ]]; then
        info "FASTQ уже скачаны, пропускаем."
        return
    fi
    info "Скачиваем риды $SAMPLE с NCBI SRA..."
    prefetch "$SAMPLE" --output-directory data/
    fasterq-dump "data/${SAMPLE}/${SAMPLE}.sra" \
        --outdir data/ \
        --threads "$THREADS" \
        --split-files
    pigz -p "$THREADS" "data/${SAMPLE}_1.fastq" "data/${SAMPLE}_2.fastq" 2>/dev/null \
        || gzip "data/${SAMPLE}_1.fastq" "data/${SAMPLE}_2.fastq"
    info "Риды скачаны: $FASTQ_1, $FASTQ_2"
}

# --------------- Шаг 2: Скачивание референса ---------------------------------
download_ref() {
    if [[ -f "$REF" ]]; then
        info "Референс уже есть, пропускаем."
        return
    fi
    info "Скачиваем референсный геном E.coli K-12..."
    wget -q -O "$REF_GZ" "$REF_URL"
    gunzip -k "$REF_GZ"
    mv "${REF_GZ%.gz}" "$REF"
    info "Референс: $REF"
}

# --------------- Шаг 3: FastQC -----------------------------------------------
run_fastqc() {
    info "Запускаем FastQC..."
    fastqc "$FASTQ_1" "$FASTQ_2" \
        --outdir "$QC_DIR" \
        --threads "$THREADS" \
        >> "$LOG" 2>&1

    for f in "$QC_DIR"/*.html; do
        mv "$f" "${f%.html}_QC-report.html"
    done
    info "FastQC готов: $QC_DIR"
}

# --------------- Шаг 4: Индексирование референса ----------------------------
index_ref() {
    if [[ -f "${REF}.bwt" ]]; then
        info "Индекс BWA уже существует, пропускаем."
        return
    fi
    info "Индексируем референс с BWA..."
    bwa index "$REF" >> "$LOG" 2>&1
    info "Индекс создан."
}

# --------------- Шаг 5: Выравнивание (BWA MEM) ------------------------------
run_bwa() {
    if [[ -f "$SAM" ]]; then
        info "SAM уже существует, пропускаем BWA."
        return
    fi
    info "Картируем риды на референс..."
    bwa mem -t "$THREADS" "$REF" "$FASTQ_1" "$FASTQ_2" \
        > "$SAM" 2>> "$LOG"
    info "Картирование завершено: $SAM"
}

# --------------- Шаг 6: Конвертация SAM -> BAM -------------------------------
run_samtools_view() {
    info "Конвертируем SAM -> BAM..."
    samtools view -bS -@ "$THREADS" "$SAM" \
        | samtools sort -@ "$THREADS" -o "$BAM"
    samtools index "$BAM"
    info "BAM готов: $BAM"
}

# --------------- Шаг 7: samtools flagstat + парсинг --------------------------
run_flagstat() {
    info "Запускаем samtools flagstat..."
    samtools flagstat "$BAM" > "$FLAGSTAT"
    cat "$FLAGSTAT" | tee -a "$LOG"
    info "Flagstat сохранён: $FLAGSTAT"
}

parse_mapped_percent() {
    MAPPED_PCT=$(grep "mapped (" "$FLAGSTAT" \
        | grep -oP '\d+\.\d+(?=%)' \
        | head -1)
    echo "$MAPPED_PCT"
}

# --------------- Шаг 8: Алгоритм оценки качества картирования ---------------
check_mapping_quality() {
    info "=========================================="
    info "  Оценка качества картирования"
    info "=========================================="

    MAPPED_PCT=$(parse_mapped_percent)

    if [[ -z "$MAPPED_PCT" ]]; then
        error "Не удалось получить % картирования из $FLAGSTAT"
    fi

    info "Процент картированных ридов: ${MAPPED_PCT}%"

    IS_OK=$(awk -v pct="$MAPPED_PCT" -v thr="$MAPPED_THRESHOLD" \
        'BEGIN { print (pct+0 > thr+0) ? "yes" : "no" }')

    if [[ "$IS_OK" == "yes" ]]; then
        echo -e "${GREEN}✓ Mapping quality: OK (${MAPPED_PCT}% > ${MAPPED_THRESHOLD}%)${NC}" | tee -a "$LOG"
        return 0
    else
        echo -e "${RED}✗ Mapping quality: not OK (${MAPPED_PCT}% ≤ ${MAPPED_THRESHOLD}%)${NC}" | tee -a "$LOG"
        return 1
    fi
}

# --------------- Шаг 9: Сортировка BAM --------------------------------------
run_sort() {
    info "Сортируем BAM..."
    samtools sort -@ "$THREADS" "$BAM" -o "$BAM_SORTED"
    samtools index "$BAM_SORTED"
    info "Отсортированный BAM: $BAM_SORTED"
}

# --------------- Шаг 10: Вариант-коллинг FreeBayes --------------------------
run_freebayes() {
    info "Запускаем FreeBayes..."
    freebayes -f "$REF" "$BAM_SORTED" > "$VCF" 2>> "$LOG"
    info "VCF готов: $VCF"
    VARIANT_COUNT=$(grep -v "^#" "$VCF" | wc -l)
    info "Найдено вариантов: $VARIANT_COUNT"
}

# =============================================================================
# MAIN
# =============================================================================
main() {
    info "=========================================="
    info "  Запуск пайплайна: $(date)"
    info "  Образец: $SAMPLE"
    info "=========================================="

    check_deps
    setup_dirs
    download_reads
    download_ref
    run_fastqc
    index_ref
    run_bwa
    run_samtools_view
    run_flagstat

    if check_mapping_quality; then
        run_sort
        run_freebayes
        info "=========================================="
        info "  Пайплайн завершён успешно!"
        info "  VCF: $VCF"
        info "  Лог: $LOG"
        info "=========================================="
    else
        warn "Пайплайн остановлен: низкое качество картирования."
        warn "Проверь качество ридов (FastQC) и параметры BWA."
        exit 1
    fi
}

main "$@"
