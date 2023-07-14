echo "#!/bin/bash

# Filter circles

# Obtain no missing column circles
awk '{if(NF==11) print}' /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/RawBed/FS$1_raw_circle.bed > /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/NoMissing/FS$1_NoMissing_circle.bed

# Obtain low conf circles
awk '{if(\$4>=1 && \$5==1)print}' /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/NoMissing/FS$1_NoMissing_circle.bed > /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/LowConf/FS$1_low_conf_circle.bed1
awk '{if(\$4>=0 && \$5>=2)print}' /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/NoMissing/FS$1_NoMissing_circle.bed > /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/LowConf/FS$1_low_conf_circle.bed2

cat /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/LowConf/FS$1_low_conf_circle.bed1 /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/LowConf/FS$1_low_conf_circle.bed2 > /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/LowConf/FS$1_low_conf_circle.bed

rm /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/LowConf/FS$1_low_conf_circle.bed1
rm /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/LowConf/FS$1_low_conf_circle.bed2

# Obtain high conf circles
awk '{if(\$11<=0.05)print}' /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/LowConf/FS$1_low_conf_circle.bed > /dellfsqd2/ST_LBI/USER/liangxue/project/eccDNA/atlas/output/MediumConf/FS$1_medium_conf_circle.bed

echo 'finish FS$1' " > Filter_circles_FS$1.sh
