# -*- coding: utf-8 -*-

import os
import collections
import argparse
from argparse import RawTextHelpFormatter
import pyfastx
import itertools
from collections import defaultdict
import subprocess

def STARCODE(inFile, OutFile, dist, thread):

    starcodeRNN = """starcode -d %d -t %d -i %s -o %s -s --seq-id  """ % (dist, thread, inFile, OutFile)
    _ = subprocess.call(starcodeRNN, shell = True)

if __name__ == '__main__':
    ''' '''
    parser = argparse.ArgumentParser(prog='barcode.py', formatter_class=RawTextHelpFormatter,
    description=''' ''')
    parser.add_argument('-i', type=str, default="", help="input fastq (gzipped or un-gzipped)")
    parser.add_argument('-bc', type=str, default="", help="input barcode fasta (gzipped or un-gzipped)")
    parser.add_argument('-bc_d', type=int, default=4, help="barcode allowed Levenshtein Distance [4]")
    parser.add_argument('-bc_num', type=int, default=100, help="barcode nunber cutoff [100]")
    parser.add_argument('-o', type=str, default="", help="output fold")    
    parser.add_argument('-prefix', type=str, default="barcode", help="output barcode prefix [barcode]")
    parser.add_argument('-t', type=int, default=8, help="number of thread [8]")
    args = parser.parse_args()
  
    out_dir = args.o
    if not os.path.isdir(out_dir): os.mkdir(out_dir) 
    out_dir_prefix ="{0}/{1}".format(args.o, args.prefix )
    if not os.path.isdir(out_dir_prefix): os.mkdir(out_dir_prefix) 
        
# -- Cluster the barcode based on the Levenshtein Distance 
    STARCODE(args.bc, "{0}/{1}.{2}".format(args.o,args.prefix, 'clr'), dist=args.bc_d, thread=args.t)
    
# -- read barcode fasta and index by pyfastx
    fa_bc = pyfastx.Fastq(args.bc)
    fq1 = pyfastx.Fastq(args.i)
    group_bc = dict()
    bc_culster=dict()

    with open( "{0}/{1}.{2}".format(args.o,args.prefix, 'clr'), 'r') as clstr:
        for line in clstr:
            if len(line.strip().split())==3:
                key,num,idx = line.strip().split() 
                idx = [int(i) - 1 for i in idx.strip().split(',') ] # 1-base to 0-base
                group_bc[key] = [fa_bc[i] for i in idx] 
                if len(group_bc[key]) >args.bc_num:
                    with open("{0}/{2}-{1}_bc.fa".format(out_dir_prefix, key, len(group_bc[key]) ), "w", encoding = "utf-8") as fout:
                        for i in range(len(group_bc[key])):
                            fout.write(f'>{group_bc[key][i].name}\n{group_bc[key][i].seq}\n')
                    fout.close()
        # -- save the reads with total barcode number - barcode - real reads number
                    with open("{0}/{2}-{1}.fq".format(out_dir_prefix,key,len(group_bc[key]) ),"w",encoding = "utf-8") as fout:
                        bc_culster_tmp=list()
                        for i in range(len(group_bc[key])):    
                            bc_culster_tmp.append(group_bc[key][i].name)
        #save name k[ey]=m64189e_210722_192328/3/ccs  
                        bc_culster[key]=list()
                        bc_culster[key].extend(list(set(bc_culster_tmp)))
                        for i in range(len(bc_culster[key])):
                            try:
                                fout.write(f'@{bc_culster[key][i]}\n{fq1[bc_culster[key][i]]}\n+\n{fq1[bc_culster[key][i]].qual}\n')
                            except KeyError:
                                pass
                    fout.close()  
