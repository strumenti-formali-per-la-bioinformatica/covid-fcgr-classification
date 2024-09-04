#!/bin/bash

# Percorsi necessari
VCF_DIR="/home/carmens01/new_vcf"  # Directory contenente i file VCF
FASTA_REF="/home/carmens01/hg19.fa"  # Percorso al file hg19.fasta
OUTPUT_DIR="/home/carmens01/output_fasta"  # Directory principale per i file di output
SORTED_VCF_DIR="/home/carmens01/sorted_vcf"  # Directory per i file VCF ordinati

# Si assicura che la directory di output esista
mkdir -p "$OUTPUT_DIR"
mkdir -p "$SORTED_VCF_DIR"

# Itera su ogni file VCF nella directory specificata
for VCF in "$VCF_DIR"/*.vcf; do
    # Determina il nome base del file per l'output
    BASENAME=$(basename "$VCF" .vcf)
    
    # Costruisce il percorso del file VCF di output ordinato
    SORTED_FILE="$SORTED_VCF_DIR/${BASENAME}_sorted.vcf"

    # Ordina il file VCF
    bcftools sort "$VCF" -o "$SORTED_FILE"
    
    # Indicizza il file VCF ordinato con GATK
    gatk IndexFeatureFile -I "$SORTED_FILE"

    # Estrae i cromosomi mutati dal file VCF e crea un file gli  intervalli
    INTERVALS_FILE="$SORTED_VCF_DIR/${BASENAME}_chr.intervals"
    grep '^chr' "$VCF" | cut -f 1 | uniq > "$INTERVALS_FILE"
    
    # Crea una directory per ogni paziente nella directory di output (per ottenere delle sottocartelle per ciascun paziente)
    PATIENT_FASTA_DIR="$OUTPUT_DIR/$BASENAME"
    mkdir -p "$PATIENT_FASTA_DIR"
    
    # Genera il FASTA , tenendo conto della  lista degli intervals creata precedentemente
    OUTPUT_FASTA="$PATIENT_FASTA_DIR/${BASENAME}_alternate.fasta"
    gatk FastaAlternateReferenceMaker -R "$FASTA_REF" -O "$OUTPUT_FASTA" -V "$SORTED_FILE" -L "$INTERVALS_FILE"
    cd "$PATIENT_FASTA_DIR"
    # Estrae il file fasta per ciascun cromosoma
    faidx --split-files "${BASENAME}_alternate.fasta" -e "lambda x: x.split()[1].split(':')[0]" --long-names

    # Per ciascun cromosoma esegue lo script python per la creazione della rappresentazione FCGR 
    PYTHON_SCRIPT="/home/carmens01/fcgr_script.py"
    python "$PYTHON_SCRIPT" "$PATIENT_FASTA_DIR"
done
