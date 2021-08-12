#!/bin/bash

# profile processing time, peak memory usage of CITE-Seq-Count
printf "script %s is being run by %s with PID %s\n" $0 $USER $$

# read list of wells
MANIFESTS="/mnt/disks/barcode-tender-manuscript/barcode_manifests"
WELLS=$( ls $MANIFESTS | grep .tsv.gz )
TAGLIST="${MANIFESTS}/taglist.csv"
FASTQ="/mnt/disks/barcode-tender-manuscript/X017_fastq"
OUTDIR="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/CITE-seq_Count_NoCorrection"
WHITELISTS="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/CITE-seq_Count_NoCorrection/whitelists"

# initialize variables to NULL
NAME=""
READ1=""
READ2=""
BARCODES=""

# call command for each library
for W in $WELLS; do
    NAME=$( echo $W | sed s'/_barcodes.*//' )

    # create arrays to store fastq filenames
    READ1=($( ls -d $FASTQ/* | grep "${NAME}-HTO" | grep _R1_ | sort ))
    READ2=($( ls -d $FASTQ/* | grep "${NAME}-HTO" | grep _R2_ | sort ))
    BARCODES="${WHITELISTS}/${NAME}_barcodes.tsv"

    # set command
    COMMAND="/shared/apps/miniconda3/bin/CITE-seq-Count -R1 $( IFS=,; echo "${READ1[*]}" ) -R2 $( IFS=,; echo "${READ2[*]}" ) -t $TAGLIST -o ${OUTDIR}/${NAME} -wl ${BARCODES} --no_umi_correction -cbf 1 -cbl 16 -umif 17 -umil 28 -cells 40000"

    echo "Runing CITE-seq-Count on ${NAME}"
    printf "\n Command: \n\n$COMMAND\n"

    # run and time command
    $( command time -v -o "/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/CITE-seq_Count_NoCorrection/CITE-seq_Count_NoCorrection_${NAME}_profile.txt" $COMMAND )

    echo "counting complete for $NAME"
    echo

done

printf "CITE-seq-Count (no UMI correction) profiling complete\n"
