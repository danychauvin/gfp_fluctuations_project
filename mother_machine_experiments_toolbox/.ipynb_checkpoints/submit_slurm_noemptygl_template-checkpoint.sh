#!/bin/bash
PREPROC_DIR_TPL='/scicore/home/nimwegen/rocasu25/MM_Data/Dany/20230228/20230228_ace_preproc_noemptygl'
RAW_PATH="/scicore/home/nimwegen/rocasu25/MM_Data/Dany/20230228/20230228_ace/20230228_ace_4/20230228_ace_4_MMStack.ome.tif"
FLATFIELD_PATH="/scicore/home/nimwegen/rocasu25/MM_Data/Dany/20230228/20230228_gfp_flatfield/20230228_gfp_flatfield_1/20230228_gfp_flatfield_1_MMStack.ome.tif"

POS_NAMES=(Pos12 Pos13 Pos14 Pos15)

GL_DETECTION_TEMPLATE_PATH="/scicore/home/nimwegen/rocasu25/MM_Data/Dany/20230228/TEMPLATE_vng124"


ROTATION=89.8
TMAX=400
IMAGE_REGISTRATION_METHOD=1
FRAMES_TO_IGNORE=()
 
NORMALIZATION_CONFIG_PATH="true"
NORMALIZATION_REGION_OFFSET=100

FORCED_INTENSITY_NORMALIZATION_RANGE=(632,4763)
 
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
