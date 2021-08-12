#!/bin/bash

# profile processing time, peak memory usage of BarCounter
printf "script %s is being run by %s with PID %s\n" $0 $USER $$

# read list of wells
MANIFESTS="/mnt/disks/barcode-tender-manuscript/barcode_manifests"
WELLS=$( ls $MANIFESTS | grep .tsv.gz )

TAGLIST="${MANIFESTS}/taglist.csv"
FASTQ="/mnt/disks/barcode-tender-manuscript/X017_fastq"
OUTDIR="/mnt/disks/barcode-tender-manuscript/ADT_count_benchmarking/BarCounter"

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
    BARCODES="${MANIFESTS}/${NAME}_barcodes.tsv.gz"


    # set command
    COMMAND="/shared/apps/barcounter -w ${BARCODES} -t ${TAGLIST} -1 $( IFS=,; echo "${READ1[*]}" ) -2 $( IFS=,; echo "${READ2[*]}" ) -o ${OUTDIR}"

    echo "Runing BarCounter on ${NAME}"
    printf "\n Command: \n\n$COMMAND\n"

    # run and time command
    $( command time -v -o "${OUTDIR}/BarCounter_${NAME}_profile.txt" $COMMAND )

    echo "counting complete for $NAME"
    echo

done

printf "BarCounter profiling complete\n"
