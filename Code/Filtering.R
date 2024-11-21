echo "#!/bin/bash
awk '{if(NF==11) print}' FS$1_raw_circle.bed > FS$1_NoMissing_circle.bed
# Obtain low conf circles
awk '{if(\$4>=1 && \$5==1)print}' FS$1_NoMissing_circle.bed > FS$1_low_conf_circle.bed1
awk '{if(\$4>=0 && \$5>=2)print}' FS$1_NoMissing_circle.bed > FS$1_low_conf_circle.bed2
cat FS$1_low_conf_circle.bed1 FS$1_low_conf_circle.bed2 > FS$1_low_conf_circle.bed
rm FS$1_low_conf_circle.bed1
rm FS$1_low_conf_circle.bed2
awk '{if(\$11<=0.05)print}' FS$1_low_conf_circle.bed > FS$1_medium_conf_circle.bed
echo 'finish FS$1' " > Filter_circles_FS$1.sh
