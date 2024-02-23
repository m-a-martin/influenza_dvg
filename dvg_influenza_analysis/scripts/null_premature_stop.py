import pandas as pd
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
try: 
	from utils import import_seqs, find_premature_stop
except:
	from scripts.utils import import_seqs, import_fasta, find_premature_stop




def find_null_premature_stop(in_dat, cds_dict, ref_dict, min_del=500, rng=np.random.default_rng(seed=111)):
	if in_dat.Segment.drop_duplicates().shape[0] > 1:
		raise Exception('multiple input genome segments')
	if in_dat.Specid.drop_duplicates().shape[0] > 1:
		raise Exception('multiple input specid')
	segment = in_dat.Segment.values[0]
	ref_seq_len = len(ref_dict[segment])
	# generate random starting locations in the first half of the genome
	starts = rng.integers(low=in_dat.start_low.iloc[0], 
		high=in_dat.start_high.iloc[0], size=in_dat.shape[0])
	stops = rng.integers(low=in_dat.stop_low.iloc[0], 
		high=in_dat.stop_high.iloc[0], 
		size=starts.size)
	random_dvg = pd.DataFrame({
		'Specid': in_dat.Specid.iloc[0],
		'Run_x': in_dat.Run_x.iloc[0],
		'Run_y': in_dat.Run_y.iloc[0],
		'Segment': segment, 
		'Start': starts, 
		'Stop':stops, 
		'Rel_support': in_dat['Rel_support'].iloc[rng.integers(in_dat.shape[0], size=in_dat.shape[0])]})
	random_dvg['premature'] = random_dvg.\
		apply(lambda k: 
			find_premature_stop(k.Start, k.Stop, cds_dict[k.Segment], ref_dict[k.Segment + '_' + k.Run_x]),
			axis=1)
	if in_dat.iloc[0].Run_y != '':
		random_dvg['premature'] = random_dvg['premature'] |  (random_dvg.\
			apply(lambda k: 
				find_premature_stop(k.Start, k.Stop, cds_dict[k.Segment], ref_dict[k.Segment + '_' + k.Run_y]),
				axis=1))
	return(random_dvg)


def run():
	parser = argparse.ArgumentParser()
	# input files
	# NEE DTO ADD ALL RUNS AND LIB
	parser.add_argument('--allRuns', default=None)
	parser.add_argument('--dat', default=None)
	parser.add_argument('--fastaDir', default=None)
	parser.add_argument('--refDir', default='ref/*/*')
	#parser.add_argument('--mapDir', default=None)
	parser.add_argument('--padSize', default=None, type=int)
	#parser.add_argument('--repEnrollID', default=None, type=int)
	args = parser.parse_args()
	#args.allRuns =  "data/all_runs.tsv"
	#args.dat = 'output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv'
	#args.mapDir = "data/*_map.tsv"
	#args.fastaDir = "*_output/consensus/*.fasta"
	#args.padSize = 210


	all_runs = pd.read_csv(args.allRuns, sep='\t')
	all_runs = all_runs[(all_runs['type'] == 'clinical')]

	# add days post symptom onset
	day_dict = {i['SPECID']:i['days_post_onset'] for 
			idx,i in all_runs.iterrows()}

	dat = pd.read_csv(args.dat, sep='\t')	

	dat = dat.\
		assign(
			days_post_onset = lambda k: k['Specid'].map(day_dict))


	# read in reference data
	ref_dict = {}
	cds_dict = {}
	fasta_files = glob.glob(args.fastaDir)
	for idx, fasta in enumerate(fasta_files):
		seqs_names, seqs = import_seqs(fasta)
		ref_dict.update({i[0].split(' ')[0] + '_' + fasta_files[idx].split('/')[-1].split('_')[0]: \
				np.array(list(i[1][args.padSize:-args.padSize].upper())) for i in zip(seqs_names, seqs)})


	for ref in glob.glob(args.refDir + '.gb'):
		with open(ref) as fp:
			gb = fp.readlines()
		name = [i.strip().split(' ')[-1] for idx, i in enumerate(gb) if i.strip()[:7] == 'VERSION'][0]
		cds_idx = [idx for idx, i in enumerate(gb) if i.strip()[:3] == 'CDS']
		cds = [gb[idx].strip().split(' ')[-1].replace('join(','').replace(')','') for 
			idx in cds_idx]
		cds_names = [gb[idx+1].strip().split('=')[1].replace('"', '') for idx in cds_idx]
		cds_pos = [np.hstack([np.arange(int(i.split('..')[0])-1, int(i.split('..')[1])) for 
			i in j.split(',')]) for j in cds]
		cds_dict[name] = {cds_names[idx]: cds_pos[idx] for idx, i in enumerate(cds_pos)}


	for ref in glob.glob(args.refDir + '.fasta'):
		if 'pad' not in ref:
			seqs_names, seqs = import_seqs(ref)
			ref_dict.update({i[0].split(' ')[0]: \
					np.array(list(i[1].upper())) for i in zip(seqs_names, seqs)})


	sample_dat = dat[['Specid', 'Run_x', 'Run_y', 'Segment', 'segment', 'Start', 'Stop', 'Rel_support']].\
		assign(Run_y = lambda k: np.where(k.Run_y.isna(), '', k.Run_y),
			Start = lambda k: k.Start - args.padSize, Stop = lambda k: k.Stop - args.padSize,
			start_low = lambda k: k.groupby('Segment')['Start'].transform(min),
			start_high = lambda k: k.groupby('Segment')['Start'].transform(max),
			stop_low = lambda k: k.groupby('Segment')['Stop'].transform(min),
			stop_high = lambda k: k.groupby('Segment')['Stop'].transform(max))

	rng=np.random.default_rng(seed=111)
	random_dvg = []
	for rep in np.arange(0,1000):
		print(rep)
		for g, g_dat in sample_dat.groupby(['Specid', 'Run_x', 'Run_y', 'Segment', 'segment']):
			random_dvg.append(find_null_premature_stop(g_dat, cds_dict, ref_dict, rng=rng).assign(rep = rep))


	random_dvg = pd.concat(random_dvg)


	random_dvg.to_csv('output/null_premature_stop.tsv', sep='\t', index=None)




if __name__ == "__main__":
    run()
