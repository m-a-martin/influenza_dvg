import pandas as pd
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
try: 
	from utils import import_seqs, find_premature_stop, plot_style, format_map_dict
except:
	from scripts.utils import import_seqs, import_fasta, find_premature_stop, plot_style, format_map_dict


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat', default=None)
	parser.add_argument('--fastaDir', default=None)
	parser.add_argument('--refDir', default='ref/*/*')
	#parser.add_argument('--mapDir', default=None)
	parser.add_argument('--padSize', default=None, type=int)
	parser.add_argument('--mapDir', default=None)
	#parser.add_argument('--repEnrollID', default=None, type=int)
	args = parser.parse_args()
	#args.allRuns =  "data/all_runs.tsv"
	#args.dat = 'output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv'
	#args.fastaDir = "*_output/consensus/*.fasta"
	#args.padSize = 210
	#args.mapDir = "data/*_map.tsv"
	
	dat = pd.read_csv(args.dat, sep='\t')	

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

	
	starts = {idx: i.Start.values-args.padSize for idx, i in dat.groupby('Segment')}
	stops = {idx: i.Stop.values-args.padSize for idx, i in dat.groupby('Segment')}

	reps = 1000

	random_dvg = dat[['Specid', 'Run_x', 'Run_y', 'Segment', 'segment', 'Rel_support']].\
			iloc[np.tile(np.arange(0,dat.shape[0]), reps), :].\
			assign(rep = np.repeat(np.arange(0,reps), dat.shape[0]))


	rng=np.random.default_rng(seed=111)
	# low value in rng.integers is inclusive, high value is exclusive
	random_dvg = random_dvg.\
		groupby('Segment', group_keys=False).\
			apply(lambda g: 
				g.assign(
					n = lambda k: int(k.shape[0]/reps))).\
			groupby('Segment', group_keys=False).\
			apply(lambda g: 
				g.assign(
					start_idx = lambda k: rng.integers(k.n, size=k.shape[0]),
					stop_idx = lambda k: rng.integers(k.n, size=k.shape[0]))).\
		assign(
			rstart = lambda k: k.apply(lambda g: starts[g.Segment][g.start_idx], axis=1) +
				rng.integers(-10,11,size=k.shape[0]),
			rstop = lambda k: k.apply(lambda g: stops[g.Segment][g.stop_idx], axis=1) +
				rng.integers(-10,11,size=k.shape[0]))

	# rejection sampling if < min number of nucleotides deleted
	where_small = np.where(random_dvg.rstop - random_dvg.rstart < 500)[0]
	while where_small.shape[0] > 0:
		# low value in rng.integers is inclusive, high value is exclusive
		random_dvg.iloc[where_small,:] = random_dvg.iloc[where_small,:].\
			groupby('Segment', group_keys=False).\
			apply(lambda g: 
				g.assign(
					start_idx = lambda k: rng.integers(k.n, size=k.shape[0]),
					stop_idx = lambda k: rng.integers(k.n, size=k.shape[0]))).\
			assign(
				rstart = lambda k: k.apply(lambda g: starts[g.Segment][g.start_idx], axis=1) +
					rng.integers(-10,11,size=k.shape[0]),
				rstop = lambda k: k.apply(lambda g: stops[g.Segment][g.stop_idx], axis=1) +
					rng.integers(-10,11,size=k.shape[0]))
		where_small = np.where(random_dvg.rstop - random_dvg.rstart < 500)[0]

	# which have null premature stop
	random_dvg['premature'] = random_dvg.\
		apply(lambda k: 
			find_premature_stop(k.rstart, k.rstop, cds_dict[k.Segment], ref_dict[k.Segment + '_' + k.Run_x]),
			axis=1)
	# map start and stop locations to universal coordinates
	map_dict = format_map_dict(glob.glob(args.mapDir))
	random_dvg['mapped_rstop'] = \
		random_dvg.apply(lambda k: map_dict[k['Segment']][k['rstop']], axis=1)
	random_dvg['mapped_rstart'] = \
		random_dvg.apply(lambda k: map_dict[k['Segment']][k['rstart']], axis=1)

	segment_order = {'PB2': 0, 'PB1': 1, 'PA': 2}
	segment_dat = \
			[pd.read_csv(i, sep='\t', index_col=0) for i in glob.glob(args.mapDir)]
	segment_len = {i.index.name: max(i.index)-(args.padSize*2) for i in segment_dat}

	bw=10
	plot_style()
	fig, axs = plt.subplots(2,3,figsize=(6.4*3,4.8*2), 
		constrained_layout=True)
	for seg, idx in segment_order.items():
		bins = np.arange(0, segment_len[seg]+bw, bw)
		seg_dat = dat.query("segment == @seg")
		seg_rdat = random_dvg.query("segment == @seg")
		_ = axs[0,idx].hist(seg_dat.mapped_start-args.padSize, bins=bins, 
			facecolor='steelblue', edgecolor='#333333', alpha=0.75, density=True,
			label='start', zorder=3)
		_ = axs[0,idx].hist(seg_dat.mapped_stop-args.padSize, bins=bins,
			facecolor='indianred', edgecolor='#333333', alpha=0.75, density=True,
			label='stop', zorder=4)
		_ = axs[0,idx].set_xlabel('')
		_ = axs[0,idx].set_ylabel('')
		_ = axs[0,idx].set_yticks([])
		_ = axs[0,idx].set_title(seg + ', empirical')
		_ = axs[1,idx].hist(seg_rdat.mapped_rstart, bins=bins, 
			facecolor='steelblue', edgecolor='#333333', alpha=0.75, density=True,
			zorder=3)
		_ = axs[1,idx].hist(seg_rdat.mapped_rstop, bins=bins,
			facecolor='indianred', edgecolor='#333333', alpha=0.75, density=True,
			zorder=4)
		_ = axs[1,idx].set_xlabel('')
		_ = axs[1,idx].set_ylabel('')
		_ = axs[1,idx].set_yticks([])
		_ = axs[1,idx].set_title(seg + ', simulated')
		


	axs[0,0].legend()
	fig.supxlabel('position (nt)', size=18)
	fig.supylabel('density', size=18)
	fig.savefig('figures/simulated_empirical_juncs.pdf')
	plt.close()
	
	random_dvg.to_csv('output/null_premature_stop.tsv', sep='\t', index=None)


if __name__ == "__main__":
    run()
