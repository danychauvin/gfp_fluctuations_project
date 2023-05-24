#!/bin/bash
PREPROC_DIR_TPL=$1
RAW_PATH=$2
FLATFIELD_PATH=$3

source $4

GL_DETECTION_TEMPLATE_PATH=$5


ROTATION=$6
TMAX=$7
IMAGE_REGISTRATION_METHOD=$8
FRAMES_TO_IGNORE=()
 
NORMALIZATION_CONFIG_PATH="$9"
NORMALIZATION_REGION_OFFSET=${10}
 
#### DO NOT EDIT BELOW THIS LINE ####
ROI_BOUNDARY_OFFSET_AT_MOTHER_CELL=0
 
# generate an array of same length as position array, with every element containing the ROTATION scalar value
ROTATIONS=$ROTATION
for f in `seq ${#POS_NAMES[*]}`
do
        f=`echo $f - 1 | bc`
        ROTATIONS[$f]=$ROTATION
done
 
module purge
module load MMPreproc
 
source mm_dispatch_preprocessing.sh
mm_dispatch_preprocessing
