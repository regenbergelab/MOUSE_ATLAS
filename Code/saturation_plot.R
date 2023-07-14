#!/bin/bash
local_path="/dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/"
sample=$1

echo "#!/usr/bin/bash"
echo "date"

echo "## randomly select 20% of the reads from the input FASTQ file"
echo "/dellfsqd2/ST_LBI/USER/liangxue/app/seqtk/seqtk sample -s100 ${local_path}/fastq/${sample}/${sample}_1.clean.fq.gz 0.2 | gzip > ${local_path}/fastq/${sample}/${sample}_1.clean.20.fq.gz"
echo "/dellfsqd2/ST_LBI/USER/liangxue/app/seqtk/seqtk sample -s100 ${local_path}/fastq/${sample}/${sample}_2.clean.fq.gz 0.2 | gzip > ${local_path}/fastq/${sample}/${sample}_2.clean.20.fq.gz"

echo "## randomly select 40% of the reads from the input FASTQ file"
echo "/dellfsqd2/ST_LBI/USER/liangxue/app/seqtk/seqtk sample -s100 ${local_path}/fastq/${sample}/${sample}_1.clean.fq.gz 0.4 | gzip > ${local_path}/fastq/${sample}/${sample}_1.clean.40.fq.gz"
echo "/dellfsqd2/ST_LBI/USER/liangxue/app/seqtk/seqtk sample -s100 ${local_path}/fastq/${sample}/${sample}_2.clean.fq.gz 0.4 | gzip > ${local_path}/fastq/${sample}/${sample}_2.clean.40.fq.gz"

echo "## randomly select 60% of the reads from the input FASTQ file"
echo "/dellfsqd2/ST_LBI/USER/liangxue/app/seqtk/seqtk sample -s100 ${local_path}/fastq/${sample}/${sample}_1.clean.fq.gz 0.6 | gzip > ${local_path}/fastq/${sample}/${sample}_1.clean.60.fq.gz"
echo "/dellfsqd2/ST_LBI/USER/liangxue/app/seqtk/seqtk sample -s100 ${local_path}/fastq/${sample}/${sample}_2.clean.fq.gz 0.6 | gzip >  ${local_path}/fastq/${sample}/${sample}_2.clean.60.fq.gz"

echo "## randomly select 80% of the reads from the input FASTQ file"
echo "/dellfsqd2/ST_LBI/USER/liangxue/app/seqtk/seqtk sample -s100 ${local_path}/fastq/${sample}/${sample}_1.clean.fq.gz 0.8 | gzip >  ${local_path}/fastq/${sample}/${sample}_1.clean.80.fq.gz"
echo "/dellfsqd2/ST_LBI/USER/liangxue/app/seqtk/seqtk sample -s100 ${local_path}/fastq/${sample}/${sample}_2.clean.fq.gz 0.8 | gzip > ${local_path}/fastq/${sample}/${sample}_2.clean.80.fq.gz"

echo "date"
echo "echo '${sample} Downsampling Done'"
