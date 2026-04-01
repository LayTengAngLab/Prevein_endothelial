#!/bin/bash

#OUTPUTDIR="/mnt/data/output/fame-output/result.v1/"
OUTPUTDIR="/mnt/data/output/fame-output/result.v1/sherry"
#DATADIR="/mnt/data/test.v1/progress/*"
DATADIR="/mnt/data/cutnrunsherry/rawdata/01.RawData/*"
FASTQ_SCREEN_GENOMES="/home/elliexi/miniconda3/envs/fameenv/share/fastq-screen-0.15.3-0/download_genomes/regular_genomes_config_file/fastq_screen.conf"

mkdir -p $OUTPUTDIR/QC
mkdir -p $OUTPUTDIR/QC/ATAQV
mkdir -p $OUTPUTDIR/processed
mkdir -p $OUTPUTDIR/bigwigs
mkdir -p $OUTPUTDIR/peaks 

function submit_job() {
experiment=$1
d=$2
FASTQ1_t="$d"/*1.fq.gz
FASTQ2_t="$d"/*2.fq.gz
BAM="$OUTPUTDIR/processed/${experiment}/${experiment}_sorted_dedup.bam"

#!/bin/bash
if [[ "$CONDA_DEFAULT_ENV" != "fameenv" ]]; then
    source /home/elliexi/miniconda3/etc/profile.d/conda.sh
    conda activate fameenv || { echo "Failed to activate Conda environment"; exit 1; }
fi

mkdir -p $OUTPUTDIR/QC/${experiment}
mkdir -p $OUTPUTDIR/processed/${experiment}
mkdir -p $OUTPUTDIR/peaks/${experiment}

echo $FASTQ1_t
echo $FASTQ2_t
echo $experiment
echo $BAM

cat $FASTQ1_t > $OUTPUTDIR/processed/${experiment}/${experiment}_1.fq.gz
cat $FASTQ2_t > $OUTPUTDIR/processed/${experiment}/${experiment}_2.fq.gz

fastqc -o $OUTPUTDIR/QC/${experiment} $OUTPUTDIR/processed/${experiment}/${experiment}_1.fq.gz $OUTPUTDIR/processed/${experiment}/${experiment}_2.fq.gz
fastq_screen --aligner bowtie2 --conf $FASTQ_SCREEN_GENOMES --outdir $OUTPUTDIR/QC/${experiment} $OUTPUTDIR/processed/${experiment}/${experiment}_1.fq.gz $OUTPUTDIR/processed/${experiment}/${experiment}_2.fq.gz

NGmerge -a -1 $OUTPUTDIR/processed/${experiment}/${experiment}_1.fq.gz -2 $OUTPUTDIR/processed/${experiment}/${experiment}_2.fq.gz -o $OUTPUTDIR/processed/${experiment}/NGMerge -v -n 8

fastqc -o $OUTPUTDIR/QC/${experiment} $OUTPUTDIR/processed/${experiment}/NGMerge_1.fastq.gz $OUTPUTDIR/processed/${experiment}/NGMerge_2.fastq.gz

bowtie2 -x /mnt/data/genomes/Human/GRCh38_noalt_as/GRCh38_noalt_as \
    -1 $OUTPUTDIR/processed/${experiment}/NGMerge_1.fastq.gz -2 $OUTPUTDIR/processed/${experiment}/NGMerge_2.fastq.gz \
    --very-sensitive -I 10 -X 2000 \
    --no-unal -p 8 \
    -S $OUTPUTDIR/processed/${experiment}/${experiment}_aligned.sam \
    > $OUTPUTDIR/QC/${experiment}/${experiment}.bowtie2_error.log 2> $OUTPUTDIR/QC/${experiment}/${experiment}.bowtie2_aligned.log

samblaster -i $OUTPUTDIR/processed/${experiment}/${experiment}_aligned.sam --removeDups | samtools view -bS - | samtools sort - -o $BAM -@ 32
samtools index $BAM -@ 32

bamCoverage -b $BAM -o $OUTPUTDIR/bigwigs/${experiment}.small.bw --binSize 5 -p max --normalizeUsing RPKM --maxFragmentLength 120
bamCoverage -b $BAM -o $OUTPUTDIR/bigwigs/${experiment}.bw --binSize 5 -p max --normalizeUsing RPKM 

macs2 callpeak -t $BAM -g 2728206152 --outdir $OUTPUTDIR/peaks/${experiment} -n ${experiment}

ataqv --peak-file $OUTPUTDIR/peaks/${experiment}/${experiment}_peaks.narrowPeak \
--name ${experiment} \
--metrics-file $OUTPUTDIR/QC/ATAQV/${experiment}.ataqv.json.gz \
--mitochondrial-reference-name "chrM" \
--tss-file /mnt/data/reffiles/refTSS_v4.1_human_coordinate.hg38.bed.txt \
human \
$OUTPUTDIR/processed/${experiment}/${experiment}_sorted_dedup.bam > $OUTPUTDIR/QC/ATAQV/${experiment}.ataqv.out
END
}

for d in $DATADIR
do
    experiment=$(basename $d)
    submit_job $experiment $d
done