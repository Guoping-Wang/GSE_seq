Scripts used for data analysis in GSE-seq under MIT license.
## Project introduction

Herein, we present the Generic Single Entity Sequencing platform (GSE-Seq), which boasts ample versatility, high throughput, and high coverage, and is enabled by an innovative workflow, addressing the critical challenges in single entities sequencing: (1) one-step manufacturing of massive barcodes, (2) degradable hydrogel-based in situ sample processing and whole genome amplification, (3) integrated in-drop library preparation, (4) compatible long-read sequencing. By GSE-Seq, we have achieved a significant milestone by enabling high-throughput, long-read single-entity profiling of dsDNA and ssDNA from single virus sequencing (SV-seq) and single bacteria sequencing (SB-seq) of the human gut and marine sediment for the first time.

## System Requirements
The customized scripts only require a standard computer with sufficient RAM to support the in-memory operations.

### OS Requirements & Hardware Requirements
This package is supported for Linux. The package has been tested on the following systems: Ubuntu 18.04
RAM: 256 GB
CPU: 128 cores

### Dependencies
fuzzysearch(version 0.7.3), pyfastx(version 2.0.0), seqkit (version 0.15.0), starcode (version 1.4)

The custom scripts supports Python >=3.6. Installation through bioconda is the recommended way to install the dependencies as it resolves non-python dependencies and uses pre-configured compilation options.
Make sure you have installed the conda packaging manager (e.g. miniconda or anaconda) before starting.

1. Create a new environment:
   conda create -n GSE_Seq python=3.6
2. Activate the conda environment using the following code:
   conda activate GSE_Seq
3. Instal the  dependencies:
   conda install -c conda-forge -c bioconda fuzzysearch=0.7.3 pyfastx=2.0.0 seqkit=0.15.0 starcode=1.4


You can also install them one by one using the following links:

https://github.com/taleinat/fuzzysearch
https://github.com/lmdu/pyfastx
https://github.com/shenwei356/seqkit
https://github.com/gui11aume/starcode

## Usage
First, please use the following code to set permissions:
chmod 775 ./1.extract_barcode.sh
chmod 775 ./2.cluster_barcode.sh

1. Extract the barcodes from the fastq file after performing basecalling on the raw Pacbio sequencing. You can edit the location of the fastq file.
	bash ./1.extract_barcode.sh

2. Cluster the barcodes.
	bash ./2.cluster_barcode.sh

## Citation
Wang, G., Zhao, L., Shi, Y., Qu, F., Ding, Y., Liu, W., Liu, C., Luo, G., Li, M., Bai, X. and Li, L., 2023. High-throughput generic single-entity sequencing using droplet microfluidics. bioRxiv, pp.2023-08.

## Author contributions
G.W., L.Z., Y.S., and F.Q. conceived the idea and methodology for the study and conducted the experiments and data collection. L.Z., G.W., Y.S., Y.D., W.L., C.L., and M.L. spearheaded the computational analysis. G.L., L.L., and X.B. contributed to the study's methodology. G.W. and Y.S. led in manuscript drafting. L.Z., F.Q., Y.H., and J.Y. revised the manuscript. J.Y. and Y.H. supervised the entire study. All authors read and approved the final manuscript.

## Acknowledgement
We would like to express our gratitude to the Yu and Megan lab members for their contribution of the project. We also appreciate the valuable contributions of Chengchen Zhang, Xia Xu, Yong Feng, Fenfen Ji, Yufeng Lin, Danyu Chen, Guangyao Chen, Eagle Chu, and Zheng Zhang. This project was supported by National Key R&D Program of China (No. 2020YFA0509200/2020YFA0509203); Research Talent Hub-Innovation and Technology Fund Hong Kong (ITS/177/21FP); RGC Research Impact Fund Hong Kong (R4032-21F); Shenzhen-Hong Kong-Macao Science and Technology Program (Category C) Shenzhen (SGDX20210823103535016); Vice-Chancellor’s Discretionary Fund Chinese University of Hong Kong; Research Grants Council of the Hong Kong Special Administrative Region, China (project no.: CUHK 14207121 and 14219922), as well as the VC Discretionary Fund, provided by the Chinese University of Hong Kong (project #: 8601014).

