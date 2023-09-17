# -*- coding: utf-8 -*-

import os
from fuzzysearch import find_near_matches
import pyfastx
import collections
import argparse
from argparse import RawTextHelpFormatter

def complement(sequence):
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return (sequence.upper())

def reverse(sequence):
    sequence = sequence.upper()
    return (sequence[::-1])

def revCom(sequence):
    sequence = complement(sequence)
    sequence = reverse(sequence)
    return (sequence)
def extract_bc():
    primer1 = "CTGCCGTATGCTGAGGTGGTC"
    primer2 = "CACTGACCTCAAGTCTGCACT"
    primer3 = revCom(primer2)
    primer4 = revCom(primer1)
    fq1 = pyfastx.Fastq(args.i)
    
    for i in range(len(fq1)):
        bc_position_pos1 = find_near_matches(primer1, str(fq1[i]), max_l_dist=args.d)
        bc_position_pos2 = find_near_matches(primer2, str(fq1[i]), max_l_dist=args.d)


        len_seq=len(fq1[i].seq)

        bc_position_neg3 = find_near_matches(primer3, str(fq1[i]), max_l_dist=args.d)
        bc_position_neg4 = find_near_matches(primer4, str(fq1[i]), max_l_dist=args.d)
        barcode_pos=[]
        barcode_neg=[]
        barcode_set=[]
        barcode_pos_qual=[]
        barcode_neg_qual=[]
        
    # for positive stand
        if len(bc_position_pos1)>0 and len(bc_position_pos2)>0:
            for b1 in range(len(bc_position_pos1)):
                for b2 in range(len(bc_position_pos2)):
                    if bc_position_pos2[b2].start-bc_position_pos1[b1].end>17 and bc_position_pos2[b2].start-bc_position_pos1[b1].end <21:
                        barcode_pos_tmp=fq1[i].seq[bc_position_pos1[b1].end:bc_position_pos2[b2].start]
                        barcode_pos.append(barcode_pos_tmp)
                        barcode_pos_qual_tmp=fq1[i].qual[bc_position_pos1[b1].end:bc_position_pos2[b2].start]
                        barcode_pos_qual.append(barcode_pos_qual_tmp)


    # for negative stand
        if len(bc_position_neg3)>0 and len(bc_position_neg4)>0:
            for b1 in range(len(bc_position_neg3)):
                for b2 in range(len(bc_position_neg4)):
                    if bc_position_neg4[b2].start-bc_position_neg3[b1].end>17 and bc_position_neg4[b2].start-bc_position_neg3[b1].end <21:
                        barcode_neg_tmp=fq1[i].seq[bc_position_neg3[b1].end:bc_position_neg4[b2].start]
                        barcode_neg.append(barcode_neg_tmp)
                        barcode_neg_qual_tmp=fq1[i].qual[bc_position_neg3[b1].end:bc_position_neg4[b2].start]
                        barcode_neg_qual.append(barcode_neg_qual_tmp)

 
    #output the barcode file
        with open("{0}/{1}_bc.fq".format(out_dir, args.prefix), "a+", encoding = "utf-8") as fout:
            for bc in range(len(barcode_pos)):
                fout.write(f'@{fq1[i].name}\n{barcode_pos[bc]}\n+\n{barcode_pos_qual[bc]}\n')
            for bc in range(len(barcode_neg)):
                fout.write(f'@{fq1[i].name}\n{revCom(barcode_neg[bc])}\n+\n{barcode_neg_qual[bc][::-1]}\n')
        fout.close()
        
    #remove the adaptor and output
        total=[]
        total=sorted(bc_position_pos1+bc_position_pos2+bc_position_neg3+bc_position_neg4)
        result=dict()
        out_fq = open("{0}/{1}.fq".format(out_dir, args.prefix), 'a')
        if len(total)>=2:
            for idx in range(len(total)-1):
                result[idx]=total[idx+1].start-total[idx].end
            max_idx=max(result,key=lambda x:result[x])
            max_len=total[max_idx+1].start-total[max_idx].end
            if total[0].start >=len_seq-total[-1].end:
                if max_len >=total[0].start:
                    out_fq.write(f'@{fq1[i].name}\n{fq1[i].seq[total[max_idx].end:total[max_idx+1].start]}\n+\n{fq1[i].qual[total[max_idx].end:total[max_idx+1].start]}\n')   
                else:
                    out_fq.write(f'@{fq1[i].name}\n{fq1[i].seq[:total[0].start]}\n+\n{fq1[i].qual[:total[0].start]}\n')   
            else:
                if max_len >=len_seq-total[-1].end:
                    out_fq.write(f'@{fq1[i].name}\n{fq1[i].seq[total[max_idx].end:total[max_idx+1].start]}\n+\n{fq1[i].qual[total[max_idx].end:total[max_idx+1].start]}\n')   
                else:
                    out_fq.write(f'@{fq1[i].name}\n{fq1[i].seq[total[-1].end:]}\n+\n{fq1[i].qual[total[-1].end:]}\n')               
        elif len(total)==1:
            if len(fq1[i].seq)/2 > total[0].start:
                out_fq.write(f'@{fq1[i].name}\n{fq1[i].seq[total[0].end:]}\n+\n{fq1[i].qual[total[0].end:]}\n')   
            else:
                out_fq.write(f'@{fq1[i].name}\n{fq1[i].seq[:total[0].start]}\n+\n{fq1[i].qual[:total[0].start]}\n')   
        out_fq.close()           

if __name__ == '__main__':
    ''' '''
    parser = argparse.ArgumentParser(prog='barcode.py', formatter_class=RawTextHelpFormatter,
    description='''

    ''')
    parser.add_argument('-i', type=str, default="", help="input fastq (gzipped or un-gzipped)")
    parser.add_argument('-o', type=str, default="", help="output fold")
    parser.add_argument('-d', type=int, default=1, help="Primer allowed Levenshtein Distance [1]")
    parser.add_argument('-prefix', type=str, default="barcode", help="output barcode prefix [barcode]")
    args = parser.parse_args()
    
    out_dir = args.o
    if not os.path.isdir(out_dir): os.mkdir(out_dir) 
    if os.path.exists("{0}/{1}.fq".format(out_dir, args.prefix)):os.remove("{0}/{1}.fq".format(out_dir, args.prefix))
    if os.path.exists("{0}/{1}_bc.fq".format(out_dir, args.prefix)):os.remove("{0}/{1}_bc.fq".format(out_dir, args.prefix))

    extract_bc()
