#### LFDP Reference library demultiplexing and DADA2 processing

#############################################################################
############################################################################
# NOTE: The first section of this script is demultiplexing, and primer trimming
# and correct naming of files with LFDP species codes in Linux. After this 
# the rest of the processing occurs in R (clearly marked below)

##### IMPORTANT !!!!! #################
# The file paths for the LINUX part of the script MUST be changed as appropriate
# for the new user. The relative paths for the R project "eDNA-LFDP.Rproj" and 
# the github LFDP-eDNA repository only function correctly for the R scripts. 



###### IN THE LINUX TERMINAL, NAVIGATE TO THE LFDP-eDNA repository and the LFDP_Reference_Library Directory
# Make sure you have access to a unix system
# install conda (or miniconda etc)
# create an environment in conda, call it a name and install:

# cutadapt
conda install -c bioconda cutadapt=1.18
# sabre
conda install -c bioconda sabre
# dos2unix
conda install -c conda-forge dos2unix
# R
conda install -c r r

# Open the terminal and navigate to your sequencing raw_data folder. For me it is here (this will be different for you)
cd /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library

echo ">>> Extracting TAR file: $TAR_FILE"
tar -xf X204SC23104261-Z01-F001.tar

## data consists of 16 folders, each with F/R sequencing files consiting of the 8 samples in that metasample
## make sure its unix
dos2unix pcr1barcodes.txt

##using sabre to demuyltiplex to sample - need to formulate all the barcode text files for each metasample - ms1 to ms16
#ms1
cd X204SC23104261-Z01-F001/01.RawData/ms1

echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms1.sh

# run shell script
bash ms1.sh


#ms2
cd ..
cd ms2
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms2.sh

bash ms2.sh


#ms3
cd ..
cd ms3
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms3.sh

bash ms3.sh

#ms4
cd ..
cd ms4
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms4.sh

bash ms4.sh

#ms5
cd ..
cd ms5
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms5.sh

bash ms5.sh

#ms6
cd ..
cd ms6
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms6.sh

bash ms6.sh

#ms7
cd ..
cd ms7
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms7.sh

bash ms7.sh

#ms8
cd ..
cd ms8
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms8.sh

bash ms8.sh

#ms9
cd ..
cd ms9
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms9.sh

bash ms9.sh

#ms10
cd ..
cd ms10
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms10.sh

bash ms10.sh

#ms11
cd ..
cd ms11
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms11.sh

bash ms11.sh

#ms12
cd ..
cd ms12
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms12.sh

bash ms12.sh

#ms13
cd ..
cd ms13
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms13.sh

bash ms13.sh

#ms14
cd ..
cd ms14
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms14.sh

bash ms14.sh

#ms15
cd ..
cd ms15
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/pcr1barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms15.sh

bash ms15.sh

#ms16
cd ..
cd ms16
echo 'for i in *1.fq.gz; do bn=${i/1.fq.gz};
sabre pe -f ${bn}1.fq.gz -r ${bn}2.fq.gz -b /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/ms16barcodes.txt -u ${bn}unassigned1.fq -w ${bn}unassigned.fq;
done' > ms16.sh

bash ms16.sh

## All samples now demultiplexed from all metasamples
## Using cut adapt to remove primers - Doing replicate by replicate and make reverse primer an exact match to the corresponding reverse barcode that was not removed by sabre
## obviously you will need to adapt these sequences for the specific index used for each replicate of each plate..  rc = reverse compliment

cd ..
cd ms1
#### Trnl base primers (i.e. tag IDs etc) are:
# Trnlg: GGGCAATCCTGAGCCAA
# Trnlh: CCATTGAGTCTCTGCACCTATC

#Tube 1 of 8 reverse barcode is GAGTGG
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGCCACTC
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): GAGTGGCCATTGAGTCTCTGCACCTATC

cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq

#Tube 2 of 8 reverse barcode is CCACGTC
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGGACGTGG
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): CCACGTCCCATTGAGTCTCTGCACCTATC

cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq

#Tube 3 of 8 reverse barcode is TTCTCAGC
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGGCTGAGAA
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): TTCTCAGCCCATTGAGTCTCTGCACCTATC

cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq

#Tube 4 of 8 reverse barcode is CTAGG
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGCCTAG
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): CTAGGCCATTGAGTCTCTGCACCTATC

cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq

#Tube 5 of 8 reverse barcode is TGCTTA
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGTAAGCA
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): TGCTTACCATTGAGTCTCTGCACCTATC

cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq

#Tube 6 of 8 reverse barcode is GCGAAGT
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGACTTCGC
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): GCGAAGTCCATTGAGTCTCTGCACCTATC

cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq

#Tube 7 of 8 reverse barcode is AATCCTAT
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGATAGGATT
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): AATCCTATCCATTGAGTCTCTGCACCTATC

cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq

#Tube 8 of 8 reverse barcode is ATCTG
# Trnlg (fwd): GGGCAATCCTGAGCCAA
# Trnlh (rcRv): GATAGGTGCAGAGACTCAATGGCAGAT
# Trnlg (rcfwd): TTGGCTCAGGATTGCCC
# Trnlh (Rv): ATCTGCCATTGAGTCTCTGCACCTATC

cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms2
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms3
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms4
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms5
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms6
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms7
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms8
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms9
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms10
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms11
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms12
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms13
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms14
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms15
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGACGTGG -G CCACGTCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube2.out1.fq.gz --untrimmed-paired-output tube2.out2.fq.gz -o tube2.trim1.fq.gz -p tube2.trim2.fq.gz tube2_f.fq tube2_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGGCTGAGAA -G TTCTCAGCCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube3.out1.fq.gz --untrimmed-paired-output tube3.out2.fq.gz -o tube3.trim1.fq.gz -p tube3.trim2.fq.gz tube3_f.fq tube3_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCTAG -G CTAGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube4.out1.fq.gz --untrimmed-paired-output tube4.out2.fq.gz -o tube4.trim1.fq.gz -p tube4.trim2.fq.gz tube4_f.fq tube4_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGTAAGCA -G TGCTTACCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube5.out1.fq.gz --untrimmed-paired-output tube5.out2.fq.gz -o tube5.trim1.fq.gz -p tube5.trim2.fq.gz tube5_f.fq tube5_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGACTTCGC -G GCGAAGTCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube6.out1.fq.gz --untrimmed-paired-output tube6.out2.fq.gz -o tube6.trim1.fq.gz -p tube6.trim2.fq.gz tube6_f.fq tube6_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGATAGGATT -G AATCCTATCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube7.out1.fq.gz --untrimmed-paired-output tube7.out2.fq.gz -o tube7.trim1.fq.gz -p tube7.trim2.fq.gz tube7_f.fq tube7_r.fq
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCAGAT -G ATCTGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube8.out1.fq.gz --untrimmed-paired-output tube8.out2.fq.gz -o tube8.trim1.fq.gz -p tube8.trim2.fq.gz tube8_f.fq tube8_r.fq

cd ..
cd ms16
cutadapt -g GGGCAATCCTGAGCCAA...GATAGGTGCAGAGACTCAATGGCCACTC -G GAGTGGCCATTGAGTCTCTGCACCTATC...TTGGCTCAGGATTGCCC --untrimmed-output tube1.out1.fq.gz --untrimmed-paired-output tube1.out2.fq.gz -o tube1.trim1.fq.gz -p tube1.trim2.fq.gz tube1_f.fq tube1_r.fq

## now all tubes trimmed of primers and only the adapter combinations Per PCR tube are retained
#an example
cd ..
cd ..
cd ..
cp ms1.txt X204SC23104261-Z01-F001/01.RawData/ms1/ms1.txt
cp ms2.txt X204SC23104261-Z01-F001/01.RawData/ms2/ms2.txt
cp ms3.txt X204SC23104261-Z01-F001/01.RawData/ms3/ms3.txt
cp ms4.txt X204SC23104261-Z01-F001/01.RawData/ms4/ms4.txt
cp ms5.txt X204SC23104261-Z01-F001/01.RawData/ms5/ms5.txt
cp ms6.txt X204SC23104261-Z01-F001/01.RawData/ms6/ms6.txt
cp ms7.txt X204SC23104261-Z01-F001/01.RawData/ms7/ms7.txt
cp ms8.txt X204SC23104261-Z01-F001/01.RawData/ms8/ms8.txt
cp ms9.txt X204SC23104261-Z01-F001/01.RawData/ms9/ms9.txt
cp ms10.txt X204SC23104261-Z01-F001/01.RawData/ms10/ms10.txt
cp ms11.txt X204SC23104261-Z01-F001/01.RawData/ms11/ms11.txt
cp ms12.txt X204SC23104261-Z01-F001/01.RawData/ms12/ms12.txt
cp ms13.txt X204SC23104261-Z01-F001/01.RawData/ms13/ms13.txt
cp ms14.txt X204SC23104261-Z01-F001/01.RawData/ms14/ms14.txt
cp ms15.txt X204SC23104261-Z01-F001/01.RawData/ms15/ms15.txt
cp ms16.txt X204SC23104261-Z01-F001/01.RawData/ms16/ms16.txt

cd X204SC23104261-Z01-F001/01.RawData/ms1
dos2unix ms1.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms1.txt

cd ..
cd ms2
dos2unix ms2.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms2.txt

cd ..
cd ms3
dos2unix ms3.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms3.txt

cd ..
cd ms4
dos2unix ms4.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms4.txt

cd ..
cd ms5
dos2unix ms5.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms5.txt

cd ..
cd ms6
dos2unix ms6.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms6.txt

cd ..
cd ms7
dos2unix ms7.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms7.txt

cd ..
cd ms8
dos2unix ms8.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms8.txt

cd ..
cd ms9
dos2unix ms9.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms9.txt

cd ..
cd ms10
dos2unix ms10.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms10.txt

cd ..
cd ms11
dos2unix ms11.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms11.txt

cd ..
cd ms12
dos2unix ms12.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms12.txt

cd ..
cd ms13
dos2unix ms13.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms13.txt

cd ..
cd ms14
dos2unix ms14.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms14.txt

cd ..
cd ms15
dos2unix ms15.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms15.txt

cd ..
cd ms16
dos2unix ms16.txt          
while read BADNAME GOODNAME
do
mv $BADNAME $GOODNAME
done < ms16.txt

## now extracting all proper demultiplexed, indexed, labeled & primer/index trimmed reads
## from each sample into another directory for DADA2
cd ..
mkdir PRtreesDADA

cd ms1
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms2
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms3
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms4
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms5
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms6
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms7
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms8
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms9
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms10
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms11
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms12
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms13
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms14
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms15
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

cd ms16
mv **/*F.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
mv **/*R.fq.gz /Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA
rm tube*
rm *.fq
cd .. 

## Moving onto dada analysis

##########################################################
################## SWTICHING TO R ########################
R

library(dada2)
packageVersion("dada2")

## Making filepath based on where the trimmed files are. Now is a good point to follow the DADA2 tutorial while trying these steps

### IF RUNNING THE TERMINAL OUTSIDE OF THE eDNA-LFDP.Rproj the FULL relative path MUST be used:
setwd("/Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA")
path <- "/Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA"

### IF RUNNING FROM INSIDE eDNA-LFDP.Rproj the relative path can be used:
setwd("/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA")
path <- "/LFDP_Reference_Library/X204SC23104261-Z01-F001/01.RawData/PRtreesDADA"

fnFs <- sort(list.files(path, pattern="F.fq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R.fq.gz", full.names = TRUE))
###tidying up sample names and replicates
sample.names <-sapply(strsplit(basename(fnFs), "_"), function(x){paste(x[[1]])})

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

## checking some quality plots 
plotQualityProfile(fnFs[70:80]) # 
plotQualityProfile(fnRs[70:80])  # 

## Filtering with standard parameters (did not truncate reads as they are already very short - no quality issues on the tails)
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,  #trimRight=100,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

head(out, n=50)

## Learning error rates on these data
## There is a problem with these error functions with NOVASEQ data. Thus we need to hack the function fitter
## Ok we have an issue here with the error models due to NOvaSeq's binning of quality scores. Implementing the solution of 
## JacobRPrice, Here: https:/github.com/benjjneb/dada2/issues/1307
## by modifying the error model, his trial 1: alter loess arguments (weights and span) & enforce monotonicity

library(magrittr)
library(dplyr)


loessErrfun_mod <- function(trans) {
  qq <- as.numeric(colnames(trans))
  est <- matrix(0, nrow=0, ncol=length(qq))
  for(nti in c("A","C","G","T")) {
    for(ntj in c("A","C","G","T")) {
      if(nti != ntj) {
        errs <- trans[paste0(nti,"2",ntj),]
        tot <- colSums(trans[paste0(nti,"2",c("A","C","G","T")),])
        rlogp <- log10((errs+1)/tot)  # 1 psuedocount for each err, but if tot=0 will give NA
        rlogp[is.infinite(rlogp)] <- NA
        df <- data.frame(q=qq, errs=errs, tot=tot, rlogp=rlogp)
        
        # original
        # ###! mod.lo <- loess(rlogp ~ q, df, weights=errs) ###!
        # mod.lo <- loess(rlogp ~ q, df, weights=tot) ###!
        # #        mod.lo <- loess(rlogp ~ q, df)
        
        # Gulliem Salazar's solution
        # https:/github.com/benjjneb/dada2/issues/938
        mod.lo <- loess(rlogp ~ q, df, weights = log10(tot),span = 2)
        
        pred <- predict(mod.lo, qq)
        maxrli <- max(which(!is.na(pred)))
        minrli <- min(which(!is.na(pred)))
        pred[seq_along(pred)>maxrli] <- pred[[maxrli]]
        pred[seq_along(pred)<minrli] <- pred[[minrli]]
        est <- rbind(est, 10^pred)
      } # if(nti != ntj)
    } # for(ntj in c("A","C","G","T"))
  } # for(nti in c("A","C","G","T"))
  
  # HACKY
  MAX_ERROR_RATE <- 0.25
  MIN_ERROR_RATE <- 1e-7
  est[est>MAX_ERROR_RATE] <- MAX_ERROR_RATE
  est[est<MIN_ERROR_RATE] <- MIN_ERROR_RATE
  
  # enforce monotonicity
  # https:/github.com/benjjneb/dada2/issues/791
  estorig <- est
  est <- est %>%
    data.frame() %>%
    mutate_all(funs(case_when(. < X40 ~ X40,
                              . >= X40 ~ .))) %>% as.matrix()
  rownames(est) <- rownames(estorig)
  colnames(est) <- colnames(estorig)
  
  # Expand the err matrix with the self-transition probs
  err <- rbind(1-colSums(est[1:3,]), est[1:3,],
               est[4,], 1-colSums(est[4:6,]), est[5:6,],
               est[7:8,], 1-colSums(est[7:9,]), est[9,],
               est[10:12,], 1-colSums(est[10:12,]))
  rownames(err) <- paste0(rep(c("A","C","G","T"), each=4), "2", c("A","C","G","T"))
  colnames(err) <- colnames(trans)
  # Return
  return(err)
}

errF <- learnErrors(
  filtFs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

errR <- learnErrors(
  filtRs,
  multithread = TRUE,
  errorEstimationFunction = loessErrfun_mod,
  verbose = TRUE
)

## Having a look at new error plots with altered learn errors function
oldFerrs <- plotErrors(errF, nominalQ=TRUE)
oldFerrs
#dev.off() 
#pdf("rv_error.pdf")
oldRerrs <- plotErrors(errR, nominalQ=TRUE)
oldRerrs

### applying dada2 core inference algorithm - using default of all libraries processed seperately - no pooling
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

## Merging the paired-ends - no pooling
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

## Making sequence table
leu_p1_seqtab <- makeSequenceTable(mergers)

## keeping track of pipeline loss up to merging tables
#### Tracking read loss through the pipeline - non pooled
getN <- function(x) sum(getUniques(x))
leu_p1_track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN))
colnames(leu_p1_track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged")
head(leu_p1_track)


write.csv(leu_p1_seqtab, "/Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/PRlibraryOTU.csv")
str(leu_p1_seqtab)

prlibotu <- leu_p1_seqtab
prlibotu1 <- as.data.frame(prlibotu)
origcols <- ncol(prlibotu1) 
prlibotu1$highestOTU<-names(prlibotu1)[apply(prlibotu1,1,which.max)]
prlibotu1$totseq<-rowSums(prlibotu1[1:121,1:origcols])
prlibotu1$tophap <- apply(prlibotu1[1:121,1:origcols], 1, max)
prlibotu1$tophapprop <- prlibotu1$tophap/prlibotu1$totseq

## making routine to find the second most abundant haplotype to compare to the proportion of the most abundant
## (Closer the proportion, more spurious the most abundant haplotype is for being thw true endogenous haplotype)
p2 <- prlibotu1
tail(prlibotu1)
dim(p2)
p2 <- p2[,1:origcols] - p2$tophap
p2[p2 == 0] <- NA
p2$secondhap <- apply(p2[,1:origcols], 1, max,na.rm=TRUE)
p2$secondhap1 <- p2$secondhap+prlibotu1$tophap

prlibotu1$secondhap <- p2$secondhap1
prlibotu1$secondhapprop <- prlibotu1$secondhap/prlibotu1$totseq

#### Load these libraries if needed:

#install.packages("remotes")
#devtools::install_github("rodrigarc/RepertoiR")

### install and load the scifer library
#BiocManager::install("scifer")
library(scifer)
df_to_fasta(row.names(prlibotu1), prlibotu1$highestOTU, file_name = "/Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/PRfreflib.fasta")

write.csv(prlibotu1,"/Users/glennd/Documents/GitHub/eDNA-LFDP/LFDP_Reference_Library/PRreflib_haplotype_data.csv")
