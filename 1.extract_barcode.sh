#!/bin/bash

###################################################################
# Extract barcode
# ID is the sample name.
# LV_distance is the 
# INPUT is the sample folder. For example, the location of the input sample Test.fastq is 0.raw/Test_fastq/Test.fastq.
# OUTPUT is the results folder.


  
  ID="Test"	
  LV_distance="1"
  INPUT="0.raw/${ID}_fastq"
  OUTPUT="1.extract_barcode/${ID}"

  if [[ ! -d $OUTPUT ]]; then
    mkdir -p $OUTPUT
  fi
  seqkit split ${INPUT}/${ID}.fastq -p 100 -f 

  for i in $(seq -w 1 100); do
     python ./1.1_extract_barcode.py \
     -i ${INPUT}/${ID}.fastq.split/${ID}.part_${i}.fastq \
     -o ${OUTPUT} \
     -d ${LV_distance} \
     -prefix ${i}_${ID} &
  done

