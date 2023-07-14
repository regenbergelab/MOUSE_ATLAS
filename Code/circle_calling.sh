# Circle calling using Circle-Map pipeline

#----------------- Creat mouse reference genome index -----------------

$ mkdir -p /dellfsqd2/ST_LBI/USER/liangxue/reference/mm10/
$ cd /dellfsqd2/ST_LBI/USER/liangxue/reference/mm10/
$ wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz .
$ gzip -d mm10.fa.gz
$ nohup bwa index mm10.fa &
$ nohup samtools faidx mm10.fa &

#-----------------  Run BWA-MEM -----------------

echo "#!/bin/bash
/dellfsqd2/ST_LBI/USER/liangxue/app/bwa-0.7.17/bwa mem -q -t 8 \
/dellfsqd2/ST_LBI/USER/liangxue/reference/mm10/mm10.fa \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/01.rawdata/FS$1/FS$1_1.clean.fq.gz \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/01.rawdata/FS$1/FS$1_2.clean.fq.gz > \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/02.alignData/FS$1/FS$1_circle.sam" > step2_bwa_$1.sh

echo "/dellfsqd2/ST_LBI/USER/liangxue/app/bwa-0.7.17/bwa mem -q -t 8 \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/plasmid_ref/plasmid.fasta \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/01.rawdata/FS$1/FS$1_1.clean.fq.gz \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/01.rawdata/FS$1/FS$1_2.clean.fq.gz > \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/02.alignData_align2plasmid/FS$1/FS$1_circle.sam" > step2_bwa_align2plasmid_$1.sh

#----------------- Run samtools -----------------

echo "#!/bin/bash
source /dellfsqd2/ST_LBI/USER/liangxue/app/miniconda3/bin/activate /dellfsqd2/ST_LBI/USER/liangxue/app/miniconda3/eccDNA2 

/dellfsqd2/ST_LBI/USER/liangxue/app/samtools-1.9/samtools sort -@ 10 -n -o \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/03.preData/FS$1/qname_FS$1_circle.bam \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/02.alignData/FS$1/FS$1_circle.sam

/dellfsqd2/ST_LBI/USER/liangxue/app/samtools-1.9/samtools sort -@ 10 -o \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/03.preData/FS$1/sorted_FS$1_circle.bam \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/02.alignData/FS$1/FS$1_circle.sam

#----------------- Run Circle-Map ReadExtractor -----------------

echo "#!/bin/bash
source /dellfsqd2/ST_LBI/USER/liangxue/app/miniconda3/bin/activate /dellfsqd2/ST_LBI/USER/liangxue/app/miniconda3/eccDNA2 

Circle-Map ReadExtractor -i ./qname_FS$1_circle.bam -o ./FS$1_circular_read_candidates.bam

/dellfsqd2/ST_LBI/USER/liangxue/app/samtools-1.9/samtools sort -@ 10 -o \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/03.preData/FS$1/FS$1_sort_circular_read_candidates.bam \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/03.preData/FS$1/FS$1_circular_read_candidates.bam

/dellfsqd2/ST_LBI/USER/liangxue/app/samtools-1.9/samtools index -@ 10 \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/03.preData/FS$1/FS$1_sort_circular_read_candidates.bam

/dellfsqd2/ST_LBI/USER/liangxue/app/samtools-1.9/samtools index -@ 10 \
/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/03.preData/FS$1/sorted_FS$1_circle.bam" > step3_samtools_$1.sh

#----------------- Run Circle-Map Realign -----------------

echo "#!/bin/bash
source /dellfsqd2/ST_LBI/USER/liangxue/app/miniconda3/bin/activate /dellfsqd2/ST_LBI/USER/liangxue/app/miniconda3/eccDNA2 

export TMPDIR=/dellfsqd2/ST_LBI/USER/liangxue/project/Circular_Atlas/tmp_dir/FS$1
mkdir -p \$TMPDIR

Circle-Map Realign -t 10 -i \
../../03.preData/FS$1/FS$1_sort_circular_read_candidates.bam -qbam \
../../03.preData/FS$1/qname_FS$1_circle.bam -sbam \
../../03.preData/FS$1/sorted_FS$1_circle.bam -fasta \
../../../../reference/mm10/mm10.fa -o \
../../04.FinallyData/FS$1/FS$1_raw_circle.bed" > step4_circle-map_$1.sh
