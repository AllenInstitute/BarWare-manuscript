#!/bin/bash

# Time statement function
stm() {
    local ts=$(date +"%Y-%m-%d %H:%M:%S")
    echo "["$ts"] "$1
}

# profile single well worth of KITE processing

# profile_KITE_helper.sh NAME

# assign variables
INDEX="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/Cell_Hash_Mismatch.idx"
FASTQ="/mnt/disks/barcode-tender-manuscript/X017_fastq"
NAME=$1
OUTDIR="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/${NAME}"
WHITELIST="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/3M-february-2018.txt"
T2G="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/Cell_Hash_Mismatch.t2g"

# create arrays to store fastq filenames
READ1=($( ls -d $FASTQ/* | grep "${NAME}-HTO" | grep _R1_ | sort ))
READ2=($( ls -d $FASTQ/* | grep "${NAME}-HTO" | grep _R2_ | sort ))
BARCODES="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/kite/${NAME}_barcodes.tsv"

# create zipped list of fastq files
unset READ_PAIRS
for ((i = 0; i < ${#READ1[*]}; i++)); do
    READ_PAIRS+=("${READ1[$i]}" "${READ2[$i]}")
done


# pseudoalign reads
echo $( stm "Aligning reads for ${NAME}" )
COMMAND="/shared/apps/kallisto/kallisto bus -i $INDEX -o $OUTDIR -x 10xv3 "${READ_PAIRS[*]}""
$( $COMMAND )
echo $( stm "Alignment for ${NAME} complete" )

# bustools correct
echo $( stm "Correcting barcodes for ${NAME}" )
$( /shared/apps/bustools/bustools correct -w $WHITELIST "${OUTDIR}/output.bus" -o "${OUTDIR}/output_corrected.bus" )
echo $( stm "Correction for ${NAME} complete" )

# bustools sort
echo $( stm "Sorting for ${NAME}" )
$( /shared/apps/bustools/bustools sort -o "${OUTDIR}/output_sorted.bus" "${OUTDIR}/output_corrected.bus" )
echo $( stm "Sorting for ${NAME} complete" )

# bustools count
echo $( stm "Counting ${NAME}" )
mkdir "${OUTDIR}/featurecounts"
$( /shared/apps/bustools/bustools count -o "${OUTDIR}/featurecounts/featurecounts" --genecounts -g $T2G -e "${OUTDIR}/matrix.ec" -t "${OUTDIR}/transcripts.txt" "${OUTDIR}/output_sorted.bus" )
echo $( stm "Counting for ${NAME} complete" )