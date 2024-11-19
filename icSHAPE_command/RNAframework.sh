######################RNAframework

Tool	Description
rf-index	Automatically queries UCSC genome database and builds the transcriptome Bowtie reference index for the RF Map module
rf-map	Performs reads pre-processing and mapping
rf-count	Calculates per-base RT-stops/mutations and coverage from transcriptome-level SAM/BAM files, and generates mutation map files for processing with DRACO
rf-count-genome	Calculates per-base RT-stops/mutations and coverage from genome-level SAM/BAM files
rf-norm	Performs whole-transcriptome normalization of structure probing data
rf-fold	Produces secondary structures for the analyzed transcripts using structure probing data to guide folding
rf-structextract	Extracts substructures based on certain criteria (e.g. median SHAPE reactivity, median Shannon entropy, stability higher than expected by chance, etc.)
rf-duplex	Analyzes direct RNA-RNA interaction mapping experiments (i.e. COMRADES, SPASH, PARIS, etc.)
rf-compare	Compares secondary structures inferred by rf-fold with a set of reference structures, and computes PPV and sensitivity
rf-jackknife	Iteratively optimize slope and intercept parameters to maximize PPV and sensitivity using a set of reference structures
rf-modcall	Performs analysis of Ψ-seq/Pseudo-seq and 2OMe-seq data
rf-peakcall	Performs peak calling of RNA immunoprecipitation (RIP) experiments
rf-motifdiscovery	Discovers significantly enriched sequence motifs in RIP peaks
rf-combine	Combines results of multiple experiments into a single profile
rf-correlate	Calculates pairwise correlation of structure probing experiments
rf-wiggle	Produces WIGGLE track files from RC or XML input files
rf-rctools	Allows manipulating RC files
rf-mutate	Designs structure disrupting (and compensatory) mutations
rf-json2rc	Post-processes DRACO JSON output files into RC files


--properly-paired

#######################rf-count###########################
#########mutation
rf-count -p 80 -wt 70 -a -o /data1/bioinfo/zhoutong/data/RASP/PRJNA608297_icSHAPE2021b/01_count_mut/mouse -t /data1/bioinfo/zhoutong/data/RASP/PRJNA608297_icSHAPE2021b/01_count_mut/mouse/Xist_tmp  -f /data1/bioinfo/zhoutong/data/RASP/Mus_Musculus_annotation_files/mouse_Xist.fasta -m --orc  /data1/bioinfo/zhoutong/data/RASP/PRJNA608297_icSHAPE2021b/mapping_known/mouse/mESN_Mus_musculus_umifastp.sam
-cc  推荐shape-map使用


#只考虑 mismatch
rf-count -p 40 -wt 50 -a -o /home/tmp/a375/count_mut/NAIN3_293ft_n1/only_mismatch  -ow  -t /home/tmp/a375/count_mut/tmp_293ft  -f /home/tmp/data/RNAframework/chr22/human_rRNA_tRNA_mtRNA.fa -m  /home/tmp/a375/18s_sam/293ft_n1.sam -nd -ni

########stop
#mf :mask file
rf-count -p 80 -wt 50 -a -o /home/tmp/data/a375/01_count_mut/a375_n1_output   -t /home/tmp/data/a375/01_count_mut/a375_n1_output/tmp_N  -f  /home/tmp/data/a375/human_rRNA_tRNA_mtRNA.fa   /home/tmp/data/a375/18s_sam/a375_n1.sam    -ow

#######查看rc文件
rf-rctools view /home/tmp/data/RNAframework/mouse_treat/18S_2_DMS_4min.rc
#转成列表
rf-rctools view  -t /home/tmp/a375/count_stop/NAIN3_293ft_n1/293ft_n1.rc > /home/tmp/data/RNAframework/mouse_untreat/mouse_18s_mutation/extract.txt
sed 's/\t/,/g' /home/tmp/data/RNAframework/mouse_untreat/mouse_18s_mutation/extract.txt > /home/tmp/data/RNAframework/mouse_untreat/mouse_18s_mutation/extract.csv

##########rf-norm
rf-norm -u /home/tmp/qq/fastp_3_16/dmso/egfp/dmso_5min.rc -t /home/tmp/qq/fastp_3_16/5min/egfp/egfp_5min.rc -o /home/tmp/qq/fastp_3_16/rawactivity_egfp_mut -ow -sm 3 -nm 2
#sm  1.Ding et al., 2014
# 2. Rouskin et al., 2014
# 3. Siegfried et al., 2014
# 4. Zubradt et al., 2016
Method for signal normalization (1-3, Default: 1):
# 1. 2-8% Normalization
# 2. 90% Winsorizing
# 3. Box-plot Normalization

rf-norm -t /home/tmp/data/single_cell/PRJNA946308/02_count_mut/single_cell_transcripts_in_hESC_NAIN3_RHE1552_b1_output/single_cell_transcripts_in_hESC_NAIN3_RHE1552_b1.rc -o /home/tmp/data/RASP/41_PRJNA646706/SRP_7SK/01_count_mut/1M4_HEK293_in_vivo_SSII_Homo_sapiens_fastp_output/rawactivity_mut_nonorm  -sm 4 -nm 2 -r 

# no norm

rf-norm -t /home/tmp/data/single_cell/PRJNA946308/02_count_mut/single_cell_transcripts_in_hESC_NAIN3_RHE1552_b1_output/single_cell_transcripts_in_hESC_NAIN3_RHE1552_b1.rc -o /home/tmp/data/single_cell/PRJNA946308/02_count_mut/single_cell_transcripts_in_hESC_NAIN3_RHE1552_b1_output/mut_score1  -sm 4 -nm 2

sm 1 stop
sm 2 stop 不考虑对照
sm 3 mutation
sm 4 mutation 不考虑对照
-r no norm Reports raw reactivities (skips data normalization)

vim human_small.xml
:1,38d
:33,35d

sed 's/\t/,/g'  human_small.xml >human_small.txt
vim human_small.txt
:%s/,,,//g

awk '{printf "%s", $0} END {print ""}' human_small.txt > human_small_new.txt




sed 's/\t/,/g'  /home/tmp/a375/count_stop/NAIN3_293ft_n1/rawactivity_stop/human_small.xml >/home/tmp/a375/count_stop/NAIN3_293ft_n1/shape_score/human_small.csv

awk -v OFS=',' '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,}' /home/zhoutong/dmsseq_dmsmapseq/sample_data/mismatches_for_SVM_training.txt > /home/zhoutong/dmsseq_dmsmapseq/sample_data/mismatches_for_SVM_training.csv


#rf-fold
rf-fold   -o /home/tmp/data/RASP/41_PRJNA646706/SRP_7SK/01_count_mut/1M4_HEK293_in_vivo_SSII_Homo_sapiens_fastp_output/rf_fold -ow -ct --folding-method 2  -p 50    /home/tmp/data/RASP/41_PRJNA646706/SRP_7SK/01_count_mut/1M4_HEK293_in_vivo_SSII_Homo_sapiens_fastp_output/rawactivity_mut_nonorm/SRP_00150.dp.xml -rs /home/tmp/RNAstructure -w











#得到AUC
icSHAPE-pipe evaluateSHAPE \
    -i /home/tmp/a375/count_mut/merge.out \
    -s /home/tmp/data/RNAframework/chr22/norm_shape/human_small.dot \
    -o /home/tmp/a375/count_mut/merge_18s.pdf


#第一次测试，stop实验方法，RNAframework计算stop——score，raw reactivity 使用1，norm使用1 ，算AUC=0.52
#第二次测试，stop实验方法，RNAframework计算mutation-score，
#第三次测试，stop实验方法，merge rc，计算

#merge
rf-combine -p 40 -o /home/tmp/data/RNAframework/merge/ -ow /home/tmp/data/RNAframework/shape_react_norm/stop/NR_003278.3.xml /home/tmp/data/RNAframework/shape_react_norm/mutation/NR_003278.3.xml -s --min-correlation


/home/tmp/data/RNAframework/shape_react_norm/stop/stop.out



rf-norm -nm 2 /home/tmp/data/RNAframework/shape_react_norm/stop/stop.out > /home/tmp/data/RNAframework/shape_react_norm/stop/normstop.out

 

################MOUSE mouse_5s
bowtie2-build /home/tmp/data/RNAframework/mouse_5s/sequence.fasta /home/tmp/data/RNAframework/mouse_5s/bowtie/MUS5s
shapermapper2 --name mouse_5s --target fa/18s.fa --amplicon --primers stm_primer.txt --overwrite --min-depth 100 --modified --R1 trim/super_1_rmdup_trim.fastq --R2 trim/super_2_rmdup_trim.fastq --out super --verbose --nproc 60 --min-depth 100

bowtie2  -x /home/tmp/data/a375/srp_7sk/b_index/srp_7sk /home/icshape/a375_293ft_11_8/2.trim/293ft_n1_rmdup_trimmed.fastq  -S /home/tmp/data/a375/srp_7sk/mapping/srp_7sk.sam  -p 80  2>Align.summary




/home/zhoutong/data/RASP/PRJNA257310/fastp_data/v65_polyA_plus_icSHAPE_in_vivo_NAI_N3_Biological_Replicate_1_fastp.fq.gz



/home/bioinfo/03_software/00_env/anaconda3/envs/icz/bin/python









