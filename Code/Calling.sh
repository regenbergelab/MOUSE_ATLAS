echo "#!/bin/bash
echo "bwa mem -q -t 8 mm10.fa FS$1/FS$1_1.clean.fq.gz FS$1_2.clean.fq.gz > FS$1_circle.sam" > step2_bwa_$1.sh
echo "bwa mem -q -t 8 plasmid.fasta FS$1_1.clean.fq.gz FS$1_2.clean.fq.gz > FS$1_circle.sam" > step2_bwa_align2plasmid_$1.sh
samtools sort -@ 10 -n -o qname_FS$1_circle.bam FS$1_circle.sam
samtools sort -@ 10 -o sorted_FS$1_circle.bam FS$1_circle.sam
Circle-Map ReadExtractor -i qname_FS$1_circle.bam -o FS$1_circular_read_candidates.bam
samtools sort -@ 10 -o FS$1_sort_circular_read_candidates.bam FS$1_circular_read_candidates.bam
samtools index -@ 10 FS$1_sort_circular_read_candidates.bam
samtools index -@ 10 sorted_FS$1_circle.bam" > step3_samtools_$1.sh
echo "export TMPDIR=/tmp_dir/FS$1; mkdir -p \$TMPDIR
Circle-Map Realign -t 10 -i FS$1_sort_circular_read_candidates.bam -qbam qname_FS$1_circle.bam -sbam sorted_FS$1_circle.bam -fasta mm10.fa -o FS$1_raw_circle.bed" > step4_circle-map_$1.sh
