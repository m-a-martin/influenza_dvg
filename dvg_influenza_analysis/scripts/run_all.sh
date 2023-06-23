# first download all data
python3 dvg_influenza_analysis/dvg_influenza_analysis/scripts/mccrone_runtable.py

# now run nextflow jobs
# in practice these were run on a cluster
nextflow run -c Influenza-virus-DI-identification-pipeline/new-config-files/vic_2012_2013.conf  \
	Influenza-virus-DI-identification-pipeline/full-pipeline-customizable-mm-final.nf
nextflow run -c Influenza-virus-DI-identification-pipeline/new-config-files/perth_2010_2012.conf  \
	Influenza-virus-DI-identification-pipeline/full-pipeline-customizable-mm-final.nf
nextflow run -c Influenza-virus-DI-identification-pipeline/new-config-files/hk_2014_2015.conf  \
	Influenza-virus-DI-identification-pipeline/full-pipeline-customizable-mm-final.nf


# align and get reference maps
# to do automate seperation of reference sequence
for seg in pb2 pb1 pa ha np na m ns
do 
	mafft --auto data/${seg}_ref.fasta > data/${seg}_ref_aln.fasta
	python3 dvg_influenza_analysis/scripts/get_pos_map.py \
		--seqs data/${seg}_ref_aln.fasta \
		--name $(echo "$seg" | tr '[:lower:]' '[:upper:]')

done

# process data
# first tabulate from runtable and ensure downloaded data is correct 
python3 dvg_influenza_analysis/scripts/get_n_sample_stats.py \
	--dataDir "dvg_influenza_analysis/*_output/Virema/*.par" \
	--libs HK2 HK7 perth perth_2 HK1 HK8 HK6 vic_2 vic \
	--runTable dvg_influenza_analysis/data/SraRunTable.txt \
	--dateDat dvg_influenza_analysis/data/elife-35962-fig1-data1-v3.csv

# this combines data across runs
# output file has one row per specid
# also adds mapped start and stop locations
python3 dvg_influenza_analysis/scripts/process_dat.py \
	--runTable dvg_influenza_analysis/data/SraRunTable.txt \
	--dataDir "dvg_influenza_analysis/*_output/Virema/*.par" \
	--mapDir "dvg_influenza_analysis/data/*_map.tsv" \
	--dpDir "dvg_influenza_analysis/*_output/depth/*_dp.tsv" \
	--segmentAcc dvg_influenza_analysis/data/acc_segment.tsv


# we want to filter the above data
# filtering: 
# 1. remove any DVGs not present in both technical replicates
# 2. remove any DVGs with <10 supporting reads
# 3. remove any DVGs present in plasmid controls
python3 dvg_influenza_analysis/scripts/filter_dat.py \
	--dat dvg_influenza_analysis/output/parsed_dvgs.tsv \
	--minReads 10 \
	--onlyShared \
	--filterPlasmid

# alternative filtering where we
# just look at DVGs in PB2, PB1, PA segments
python3 dvg_influenza_analysis/scripts/filter_dat.py \
	--dat dvg_influenza_analysis/output/parsed_dvgs.tsv \
	--minReads 10 \
	--onlyShared \
	--segments PB2 PB1 PA \
	--filterPlasmid

# alternative filtering where we additionally just look at DVGs 
# with >500nt deleted and only in PB2, PB1, PA segments
python3 dvg_influenza_analysis/scripts/filter_dat.py \
	--dat dvg_influenza_analysis/output/parsed_dvgs.tsv \
	--minReads 10 \
	--onlyShared \
	--segments PB2 PB1 PA \
	--minDel 500 \
	--filterPlasmid

python3 dvg_influenza_analysis/scripts/filter_dat.py \
	--dat dvg_influenza_analysis/output/parsed_dvgs.tsv \
	--minReads 10 \
	--onlyShared \
	--minDel 500 \
	--filterPlasmid

# figures
python3 dvg_influenza_analysis/scripts/plot_figure_1.py \
	--dataDir 'dvg_influenza_analysis/*_output/Virema/*.par' \
	--allRuns dvg_influenza_analysis/data/all_runs.tsv

python3 dvg_influenza_analysis/scripts/plot_figure_2.py \
	--dat dvg_influenza_analysis/output/parsed_dvgs_10_True_True_0_all.tsv \
	--mapDir "dvg_influenza_analysis/data/*_map.tsv" \
	--padSize 210 \
	--repSpecid HS1530 \
	--allRuns dvg_influenza_analysis/data/all_runs.tsv 

python3 dvg_influenza_analysis/scripts/plot_figure_3.py \
	--dat dvg_influenza_analysis/output/parsed_dvgs_10_True_True_500_PB2_PB1_PA.tsv \
	--padSize 210 \
	--mapDir "dvg_influenza_analysis/data/*_map.tsv" \
	--allRuns dvg_influenza_analysis/data/all_runs.tsv \
	--repEnrollID 50319

python3 dvg_influenza_analysis/scripts/plot_figure_4.py \
	--dat dvg_influenza_analysis/output/parsed_dvgs_10_True_True_500_PB2_PB1_PA.tsv \
	--allRuns dvg_influenza_analysis/data/all_runs.tsv \
	--pairDat dvg_influenza_analysis/data/elife-35962-fig3-data1-v3.csv

python3 dvg_influenza_analysis/scripts/plot_figure_s1.py \
	--dat dvg_influenza_analysis/output/parsed_dvgs_10_True_True_0_all.tsv \
	--padSize 210 \
	--mapDir "dvg_influenza_analysis/data/*_map.tsv" \
	--allRuns dvg_influenza_analysis/data/all_runs.tsv

python3 dvg_influenza_analysis/scripts/plot_figure_s2.py \
	--dat dvg_influenza_analysis/output/parsed_dvgs_10_True_True_0_PB2_PB1_PA.tsv \
	--padSize 210 \
	--mapDir "dvg_influenza_analysis/data/*_map.tsv" \
	--allRuns dvg_influenza_analysis/data/all_runs.tsv       

python3 dvg_influenza_analysis/scripts/plot_figure_s3.py \
	--dat dvg_influenza_analysis/output/parsed_dvgs_10_True_True_0_PB2_PB1_PA.tsv

python3 dvg_influenza_analysis/scripts/plot_figure_s4.py \
	--allRuns dvg_influenza_analysis/data/all_runs.tsv \
	--dat dvg_influenza_analysis/output/parsed_dvgs_10_True_True_500_PB2_PB1_PA.tsv \
	--mapDir "dvg_influenza_analysis/data/*_map.tsv" \
	--padSize 210 \
	--minTspan 1

python3 dvg_influenza_analysis/scripts/plot_figure_s5.py \
	--allRuns dvg_influenza_analysis/data/all_runs.tsv \
	--dat dvg_influenza_analysis/output/parsed_dvgs_10_True_True_500_PB2_PB1_PA.tsv


