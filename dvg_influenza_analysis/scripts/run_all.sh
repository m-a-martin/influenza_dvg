# align and get reference maps
# to do automate seperation of reference sequence
for seg in pb2 pb1 pa ha np na m ns
do 
	mafft --auto data/${seg}_ref.fasta > data/${seg}_ref_aln.fasta
	python3 scripts/get_pos_map.py \
		--seqs data/${seg}_ref_aln.fasta \
		--name $(echo "$seg" | tr '[:lower:]' '[:upper:]')

done

# process data
# first tabulate from runtable and ensure downloaded data is correct 
python3 scripts/get_n_sample_stats.py \
	--dataDir "*_output/Virema/*.par" \
	--libs HK2 HK7 perth perth_2 HK1 HK8 HK6 vic_2 vic \
	--runTable data/SraRunTable.txt \
	--dateDat data/elife-35962-fig1-data1-v3.csv

# this combines data across runs
# output file has one row per specid
# also adds mapped start and stop locations
python3 scripts/process_dat.py \
	--runTable data/SraRunTable.txt \
	--dataDir "*_output/Virema/*.par" \
	--mapDir "data/*_map.tsv" \
	--dpDir "*_output/depth/*_dp.tsv" \
	--segmentAcc data/acc_segment.tsv \
	--titerTable data/elife-35962-fig1-data1-v3.csv

# we want to filter the above data
# filtering: 
# 1. remove any DVGs not present in both technical replicates
# 3. remove any DVGs present in plasmid controls
python3 scripts/filter_dat.py \
	--dat output/parsed_dvgs.tsv \
	--minRelReads 0 \
	--onlyShared \
	--filterPlasmid

# 2. remove any DVGs with 5E-3 rel supporting reads
python3 scripts/filter_dat.py \
	--dat output/parsed_dvgs.tsv \
	--minRelReads 5E-3 \
	--onlyShared \
	--filterPlasmid

# alternative filtering where we
# just look at DVGs in PB2, PB1, PA segments
python3 scripts/filter_dat.py \
	--dat output/parsed_dvgs.tsv \
	--minRelReads 5E-3 \
	--onlyShared \
	--segments PB2 PB1 PA \
	--filterPlasmid

# alternative filtering where we additionally just look at DVGs 
# with >500nt deleted and only in PB2, PB1, PA segments
python3 scripts/filter_dat.py \
	--dat output/parsed_dvgs.tsv \
	--minRelReads 5E-3 \
	--onlyShared \
	--segments PB2 PB1 PA \
	--minDel 500 \
	--filterPlasmid

# compare filtered data to raw data
python3 scripts/compare_filtered_unfiltered_n.py \
	--filteredDat 'output/parsed_dvgs_0.005_True_True_0_all.tsv' \
	--rawDat  'output/parsed_dvgs.tsv' \
	--minRelReads 5E-3 \
	--padSize  210 \
	--onlyShared


python3 scripts/null_premature_stop.py \
	--dat output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv \
	--fastaDir "*_output/consensus/*.fasta" \
	--padSize 210 \
	--mapDir "data/*_map.tsv"


# figures
python3 scripts/plot_figure_1.py \
	--dat output/parsed_dvgs_0.005_True_True_0_all.tsv \
	--mapDir "data/*_map.tsv" \
	--padSize 210 \
	--repSpecid HS1530 \
	--allRuns data/all_runs.tsv 

python3 scripts/plot_figure_2.py \
	--dat output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv \
	--padSize 210 \
	--mapDir "data/*_map.tsv" \
	--allRuns data/all_runs.tsv \
	--repEnrollID 50319

python3 scripts/plot_figure_3.py \
	--dat output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv \
	--mapDir  "data/*_map.tsv" \
	--padSize 210 \
	--allRuns data/all_runs.tsv

python3 scripts/plot_figure_4.py \
	--dat output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv \
	--allRuns data/all_runs.tsv

python3 scripts/plot_figure_5.py \
	--dat output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv \
	--allRuns data/all_runs.tsv \
	--pairDat data/elife-35962-fig3-data1-v3.csv \
	--doubleDat data/transmission_pairs.csv

python3 scripts/plot_figure_s1.py \
	--dataDir '*_output/Virema/*.par' \
	--allRuns data/all_runs.tsv \
	--dpDir "*_output/depth/*_dp.tsv"

python3 scripts/plot_figure_s2.py \
	--dat output/parsed_dvgs_0.005_True_True_0_all.tsv \
	--allRuns data/all_runs.tsv

python3 scripts/plot_figure_s3.py \
	--dat output/parsed_dvgs_0.005_True_True_0_all.tsv \
	--padSize 210

python3 scripts/plot_figure_s4.py \
	--dat output/parsed_dvgs_0.005_True_True_0_all.tsv \
	--padSize 210 \
	--fastaDir "*_output/consensus/*.fasta" \
	--allRuns data/all_runs.tsv \
	--nullPremature output/null_premature_stop.tsv


python3 scripts/plot_figure_s5.py \
	--dat output/parsed_dvgs_0.005_True_True_0_PB2_PB1_PA.tsv \
	--allRuns data/all_runs.tsv \
	--mapDir "data/*_map.tsv" \
	--padSize 210


python3 scripts/plot_figure_s6.py \
	--dat output/parsed_dvgs_0.005_True_True_0_PB2_PB1_PA.tsv 


python3 scripts/plot_figure_s7.py \
	--allRuns data/all_runs.tsv \
	--dat output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv \
	--mapDir "data/*_map.tsv" \
	--padSize 210 \
	--minTspan 1


python3 scripts/plot_figure_s8.py \
	--allRuns data/all_runs.tsv \
	--dat output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv 


python3 scripts/plot_figure_s9.py \
	--allRuns data/all_runs.tsv \
	--dat output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv  \
	--mapDir "data/*_map.tsv" \
	--padSize 210


python3 scripts/plot_figure_s10.py \
	--dat output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv 


python3 scripts/plot_figure_s11.py \
	--dat output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv \
	--fastaDir "*_output/consensus/*.fasta" \
	--nullPremature output/null_premature_stop.tsv \
	--padSize 210



