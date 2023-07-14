# Generate random eccDNA datasets

#======================== Requirements ========================
# Input file A: Here we take the output BED files from the Circle-Map as input file; you can also take any one tab delimited file with chromosome ID, start coordinates, end coordinates as the first there
columns.

# Input file B: Download the mouse chromosome size file named `mm10.chrom.sizes` from this link(https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes)
and then remove unknown chromosomes and just keep known 21 chromosomes with the `grep` command.

$ grep -v 'chrM' mm10.chrom.sizes |grep -v '_' mm10.chrom.sizes > final_mm10.chrom.sizes

#======================== Getting started ========================
# Here is the step-by-step guidance to run the pipeline.

# Step1: Import the output BED file from the Circle-Map.
data <- read.csv("my_CircleMap_output.bed")

# Step2: Randomly select 5,000 circles from all detected real circles
data_sub5000 <- data[sample(nrow(data), 5000),]

# Step3: Calculated the length of each random circle
data_sub5000$size_sub5000 <- data_sub5000$V3 - data_sub5000$V2
size <- as.data.frame(data_sub5000$size_sub5000)

# Step4: Save the size of 5000 random circles into a text file
write.table(size, file = "size_1.txt", quote = F, col.names = F, row.names = F)

# Note: Repeat Step1 - Step4 10 times for generating 10 independent random datasets.

# Step5: Create the shell script for generating random regions using bedtools random.
$ vim Generate_random_dataset.sh
for size in `cat $1`;
do
bedtools random -g final_mm10.chrom.sizes -n 1 -l $size >> Random_dataset_$2.csv;
done

# Note:
The `-n` option: specify the number of regions to generate. Here we want to generate 1
random circle for each real circle.
The `-l` option: specify the length of regions to generate. Here we want to generate random
circles of the same length as each real circle.

# Step6: Perform above shell script 10 times for generating 10 independent random datasets.
$ sh Generate_random_dataset.sh size_1.txt 1
$ sh Generate_random_dataset.sh size_1.txt 1
$ sh Generate_random_dataset.sh size_2.txt 2
$ sh Generate_random_dataset.sh size_3.txt 3
$ sh Generate_random_dataset.sh size_4.txt 4
$ sh Generate_random_dataset.sh size_5.txt 5
$ sh Generate_random_dataset.sh size_7.txt 7
$ sh Generate_random_dataset.sh size_8.txt 8
$ sh Generate_random_dataset.sh size_9.txt 9
$ sh Generate_random_dataset.sh size_10.txt 10

# Step7: check one of the random datasets.
$ less Random_dataset_1.csv | head -5
