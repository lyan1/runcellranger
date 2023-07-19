#!/bin/bash

# Developed by Jielin Yu and Le Yan
# Version 0.1

FASTA_LINK=$1
GTF_LINK=$2
#JOB_NAME=$3
JOB_DIR=$3
GENOME_NAME=$4
USE_EXISTING_REF_GENOME=$5
REF_GENOME=$6 
#num_cpus=$6 
USE_MULTI=$7
JOB_NAME=$8

#check if job dir exists
#if [ -d "$JOB_DIR" ] 
#then
#    echo "Job direxctory exists." 
#else
#    echo "Error: Job directory does not exists. Please specify a valid path of your job directory."
#    exit 1
#fi

#cd to job_dir and mkdir cellranger, exit if cellranger folder exists
cd $JOB_DIR

if [ -d "cellranger" ]
then
    echo "Cellranger folder exists. If you want to start a new run, please delete the cellranger folder from the previous run."
    exit 1
fi

mkdir cellranger
export WORK_DIR=$JOB_DIR"/"cellranger"/"

#preperation for mkref

if [ $USE_EXISTING_REF_GENOME = "No" ]
then

MKREF_DIR=$JOB_DIR"/"cellranger"/"mkref"/"
mkdir -p $MKREF_DIR

cd $MKREF_DIR

#Download files
if [ $FASTA_LINK = emp ]
then
	FASTA_FILE=`find $JOB_DIR -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fna" \)`
else
	cd $MKREF_DIR
	wget $FASTA_LINK
	gunzip `ls -t | head -1`
	FASTA_FILE=`find -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fna" \)| cut -d"/" -f2`
fi

if [ $GTF_LINK = emp ]
then
	GTF_FILE=`find $JOB_DIR -maxdepth 1 -type f \( -name "*.gtf" \)`
	GTF=`echo $GTF_FILE | rev | cut -d"." -f2-  | cut -d"/" -f1 | rev`
	GTF_FILTERED=$GTF.filtered.gtf
else
	cd $MKREF_DIR
	wget $GTF_LINK
	gunzip `find *.gtf.*`
	GTF_FILE=`find *.gtf`
	GTF=`echo $GTF_FILE | rev | cut -d"." -f2-  | rev`
	GTF_FILTERED=$GTF.filtered.gtf
fi

if [ $GENOME_NAME = emp ]
then
	GENOME_NAME=$GTF
fi

fi

#Functions

declare -A mkfastq
declare -A multi

read_mkfastq() {
    local i=0
    local lines
    local j
    read first_line
    while read -r lines; do
        j=0
        for v in $(echo -e $lines | sed "s/,/ /g" | sed "s/*/any/g"); do
            mkfastq[$i,$j]="$v"
            j=$((j+1))
        done
        i=$((i+1))
    done
}

read_multi() {
    local i=0
    local lines
    local j
    while read -r lines; do
        j=0
        for v in $(echo -e $lines | sed "s/,/ /g" | sed "s/*/any/g"); do
            multi[$i,$j]="$v"
            j=$((j+1))
        done
        i=$((i+1))
    done
}

run_mkgtf () {
    echo " =================="
    echo "  CELLRANGER mkgtf "
    echo " =================="
    cellranger mkgtf $GTF_FILE $GTF_FILTERED               
    if [ $? -eq 0 ]
    then
        echo "mkgtf successful"
        #rsync -av --progress $JOB_DIR/$ID_NAME $DATA_DIR > rsync_mkfastq.out 2>&1 
        #cp rsync_mkfastq.out $DATA_DIR/.
    fi
}

run_mkref () {
    echo " =================="
    echo "  CELLRANGER mkref "
    echo " =================="
    cellranger mkref --genome=$GENOME_NAME --fasta=$FASTA_FILE --genes=$GTF_FILTERED
    if [ $? -eq 0 ]
    then
        echo "mkref successful"
        #rsync -av --progress $JOB_DIR/$ID_NAME $DATA_DIR > rsync_mkfastq.out 2>&1
        #cp rsync_mkfastq.out $DATA_DIR/.
    fi
}

run_mkfastq () {
    
    cd $WORK_DIR
    #CSV_FILE=`find $JOB_DIR -maxdepth 1 -type f \( -name "*.csv" \)`

    echo " ========================"
    echo "    CELLRANGER mkfastq   "
    echo " ========================"

    cellranger mkfastq --run=$JOB_DIR --csv=${JOB_DIR}/${JOB_NAME}_mkfastq.csv --jobmode=local --localcores=16 --localmem=29 --id=mkfastq
    if [ $? -eq 0 ]
    then
        echo "mkfastq successful"
    fi
}

run_count () {
    
    cd $WORK_DIR
    mkdir count
    cd count
    for files in "${WORK_DIR}/mkfastq/outs/fastq_path"/*
    do
        if [[ `basename $files` != "Reports" ]] && [[ `basename $files` != "Stats" ]] &&[[ `basename $files` != *"Undetermined"* ]]
        then
            echo $files
            FASTQ_FOLDER_NAME=`basename $files`
        fi
    done

    echo " ======================"
    echo "    CELLRANGER count   "
    echo " ======================"
    for fastq_sample in "${WORK_DIR}/mkfastq/outs/fastq_path/${FASTQ_FOLDER_NAME}"/*
    do
    cellranger count --fastqs=${WORK_DIR}/mkfastq/outs/fastq_path/ --sample=`basename ${fastq_sample}` --transcriptome=${REF_GENOME} --id=`basename ${fastq_sample}` --jobmode=local --maxjobs=16 --localmem=29
    done

    if [ $? -eq 0 ]
    then
        echo "count successful"
    fi
}

run_multi () {
    cd $WORK_DIR
    mkdir multi
    cd multi

    for files in "${WORK_DIR}/mkfastq/outs/fastq_path"/*
    do
        if [[ `basename $files` != "Reports" ]] && [[ `basename $files` != "Stats" ]] &&[[ `basename $files` != *"Undetermined"* ]]
        then
            echo $files
            export FASTQ_FOLDER_NAME=`basename $files`
        fi
    done

    while read -r lines;do
	samples=`echo -e $lines | cut -d, -f1`
	echo "[gene-expression]" > ${samples}.csv
	echo "reference,${REF_GENOME}" >> ${samples}.csv
	echo "[vdj]" >> ${samples}.csv
	echo "reference,${REF_GENOME}/vdj" >> ${samples}.csv
	echo "[libraries]" >> ${samples}.csv
	echo "fastq_id,fastqs,lanes,feature_types" >> ${samples}.csv
	echo $lines | envsubst >> ${samples}.csv

    echo " ============================================"
    echo "    CELLRANGER multi for sample ${samples}   "
    echo " ============================================"
    
    cellranger multi --id=${samples} --csv=${samples}.csv

    done < ${JOB_DIR}/${JOB_NAME}_multi.csv

    if [ $? -eq 0 ]
    then
        echo "multi successful"
    fi

}

run_multi_aggr () {

    cd $WORK_DIR
    mkdir aggr
    cd aggr
    cat ${JOB_DIR}/${JOB_NAME}_aggr.csv | envsubst > aggr.csv

    echo " ====================="
    echo "    CELLRANGER aggr   "
    echo " ====================="
    cellranger aggr --id=aggr --csv=aggr.csv --jobmode=local --maxjobs=16 --localmem=29

    if [ $? -eq 0 ]
    then
        echo "aggregate successful"
    fi

}

run_aggr () {
    
    cd $WORK_DIR
    mkdir aggr
    cd aggr
    touch aggr.csv
    echo "sample_id,molecule_h5" >> aggr.csv

    for count_result in "${WORK_DIR}/count"/*
    do
        if [[ `basename $count_result` != *"_"* ]]
        then
            echo "`basename $count_result`,${count_result}/outs/molecule_info.h5" >> aggr.csv
        fi
    done

    
    echo " ====================="
    echo "    CELLRANGER aggr   "
    echo " ====================="
    cellranger aggr --id=aggregate --csv=${WORK_DIR}/aggr/aggr.csv --jobmode=local --maxjobs=16 --localmem=29

    if [ $? -eq 0 ]
    then
        echo "aggregate successful"
    fi
}


echo " ================="
echo "  BEGIN EXECUTION "
echo " ================="

if [ $USE_EXISTING_REF_GENOME = "No" ]
then
run_mkgtf
run_mkref
fi

run_mkfastq

if [ $USE_MULTI = "Yes" ]
then
run_multi
run_multi_aggr
cd $WORK_DIR
tar -cf cellranger.tar aggr mkfastq multi
else
run_count
run_aggr
cd $WORK_DIR
tar -cf cellranger.tar aggr mkfastq count
fi

echo "================="
echo "  END EXECUTION  "
echo "================="

#ls -lR $DATA_DIR > $BCL_NAME.filelists.txt 
#cp $DATA_DIR/$BCL_NAME.filelists.txt .
#rm -rf $JOB_DIR/*
#cp $DATA_DIR/$BCL_NAME.filelists.txt $JOB_DIR/.

echo " ==============================="
echo " JOB DIRECTORY IS = $JOB_DIR" 
echo " ==============================="

