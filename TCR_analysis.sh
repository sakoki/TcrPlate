#!/bin/bash

# Calcutate runtimedata
STARTTIME=$(date +%s)

# Set up arguments
currentDate=`date +%m-%d-%y`
inputDir="./data/"
outputDir="./data/output/${currentDate}/"
echo Date: $currentDate

# Run analysis
Rscript ./R/TCR_clonality_figures.R \
--input "${inputDir}fullTCRseq.csv" \
--output $outputDir \
--select_pt_loc "Patient_S, Patient_D, Patient_U, Patient_AA, Patient_AB, Patient_AC" \
--select_pt_sev "Patient_D, Patient_U, Patient_AA, Patient_Y"

ENDTIME=$(date +%s)
RUNTIME=$(expr $ENDTIME - $STARTTIME)

echo "*****Done*****"
echo Task completed in $RUNTIME s
