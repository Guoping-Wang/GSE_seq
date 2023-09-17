#!/bin/bash

set -e


  ID="Test"
  LV_distance="1"
  Bc_num="30"
  INPUT="1.extract_barcode/${ID}"
  OUTPUT="2.cluster_barcode/${ID}"

  if [[ ! -d $OUTPUT ]]; then
    mkdir -p $OUTPUT
  fi
  
  cat ${INPUT}/*_${ID}_bc.fq >${INPUT}/${ID}_total_bc.fq
  cat ${INPUT}/*_${ID}.fq >${INPUT}/${ID}_total.fq
  
 
  python ./2.1_cluster_barcode.py \
         -i ${INPUT}/${ID}_total.fq \
         -bc ${INPUT}/${ID}_total_bc.fq \
         -bc_d ${LV_distance} \
	 -bc_num ${Bc_num} \
         -o ${OUTPUT} \
         -prefix ${ID} \
         -t 60 

