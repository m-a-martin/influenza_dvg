import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import string
import glob
from scipy.stats import kruskal
from collections import defaultdict
try: 
	from utils import plot_style, jitter_boxplot, import_seqs
except:
	from scripts.utils import plot_style, jitter_boxplot, import_seqs


def consensus_fuzz(del_start, del_stop, ref_seq):
		# start and stop are 1-indexed
		# adjust start and stop
		del_start -= 1
		del_stop -= 1
		# start and stop represent the last non-deleted base on the 3' end
		# and the first non-deleted base on the 5' end
		# arange doesn't include end point so don't need +1 
		# assumes junctions all pushed towards end of reference
		# first, undeleted portion of the beginning of the genome
		seg1 = ref_seq[:del_start+1]
		# second, deleted nucleotides 
		seg2 = ref_seq[del_start+1:del_stop]
		go = True
		fuzz = 0
		while go == True:
			if (seg1[-(fuzz+1):] == seg2[-(fuzz+1):]).all():
				fuzz += 1
			else:
				go = False
		return(fuzz)
		


def run():
	parser = argparse.ArgumentParser()
	# input files
	# NEE DTO ADD ALL RUNS AND LIB
	parser.add_argument('--allRuns', default=None)
	parser.add_argument('--dat', default=None)
	parser.add_argument('--fastaDir', default=None)
	parser.add_argument('--nullPremature', default=None)
	parser.add_argument('--padSize', default=0, type=int)
	args = parser.parse_args()
	#args.dat = 'output/parsed_dvgs_0.005_True_True_0_all.tsv'
	#args.fastaDir = "*_output/consensus/*.fasta"
	#args.nullPremature = "output/null_premature_stop.tsv"
	#args.padSize = 210
	# read in reference data
	ref_dict = defaultdict(lambda: np.nan)
	fasta_files = glob.glob(args.fastaDir)
	for idx, fasta in enumerate(fasta_files):
		seqs_names, seqs = import_seqs(fasta)
		ref_dict.update({i[0].split(' ')[0] + '_' + fasta_files[idx].split('/')[-1].split('_')[0]: \
				np.array(list(i[1].upper())) for i in zip(seqs_names, seqs)})

	#read in data and remove pads
	dat = pd.read_csv(args.dat, sep='\t')

	# start and stop are inclusive of pad
	# reference sequences include pad
	dat['fuzz1'] = \
		dat.apply(lambda g: consensus_fuzz(g.Start, g.Stop, ref_dict[g.Segment + '_' + g.Run_x]), axis=1)
	dat['fuzz2'] = \
		dat.apply(lambda g: 
			consensus_fuzz(g.Start, g.Stop, ref_dict[g.Segment + '_' + g.Run_x]) if 
				type(g.Run_y) == str else np.nan, axis=1)
	dat['fuzz'] = dat[['fuzz1', 'fuzz2']].max(axis=1)

	# get null
	null_dvg = pd.read_csv(args.nullPremature, sep='\t')
	null_dvg['fuzz1'] = \
		null_dvg.apply(lambda g: consensus_fuzz(g.rstart, g.rstop, ref_dict[g.Segment + '_' + g.Run_x]), axis=1)

	null_dvg['fuzz2'] = \
		null_dvg.apply(lambda g: 
			consensus_fuzz(g.rstart, g.rstop, ref_dict[g.Segment + '_' + g.Run_x]) if 
				type(g.Run_y) == str else np.nan, axis=1)

	null_dvg['fuzz'] = null_dvg[['fuzz1', 'fuzz2']].max(axis=1)
	null_quantiles = null_dvg.groupby(['rep', 'fuzz']).size().reset_index().\
		groupby('fuzz').\
		apply(lambda k: np.quantile(k[0], [0, 0.025, 0.5, 0.975,1]))
	null_quantiles = pd.DataFrame(null_quantiles.to_list(), 
		index=null_quantiles.index, 
		columns=['q0', 'q0.025', 'q0.5', 'q0.975','q1'])


	fuzz_sum = dat.groupby('fuzz').size()
	ufuzz = dat.groupby(['segment', 'mapped_start', 'mapped_stop', 'fuzz']).size().reset_index().\
		rename(columns={0:'count'})
	ufuzz_sum = ufuzz.groupby('fuzz').size()

	pol_fuzz_sum = dat.query('(segment == "PB2" | segment == "PA" | segment == "PB1") & (Stop - Start >= 500)').\
		groupby('fuzz').size()
	pol_ufuzz = dat.query('(segment == "PB2" | segment == "PA" | segment == "PB1") & (Stop - Start >= 500)').\
		groupby(['segment', 'mapped_start', 'mapped_stop', 'fuzz']).size().reset_index().\
		rename(columns={0:'count'})
	pol_ufuzz_sum = pol_ufuzz.groupby('fuzz').size()
	
	r_ufuzz = null_dvg.\
		groupby(['rep', 'segment', 'mapped_rstart', 'mapped_rstop', 'fuzz']).size().reset_index().\
		rename(columns={0:'count'})
	r_ufuzz_sum = r_ufuzz.groupby('fuzz').\
		apply(lambda k: np.quantile(k['count'], [0, 0.025, 0.5, 0.975,1]))

	plot_style()


	common = dat.\
		groupby(['segment', 'mapped_start', 'mapped_stop', 'fuzz']).\
		size().\
		reset_index().\
		sort_values(by=0, ascending=False).iloc[0,:].\
		rename({0:'n'})
	n = dat[['segment', 'mapped_start', 'mapped_stop', 'fuzz']].drop_duplicates().shape[0]
	output = []
	output.append(f'there are a total of {n} unique DVG species')
	output.append(f'the most common DVG is {common.segment} {common.mapped_start - args.padSize}_{common.mapped_stop-args.padSize}')
	output.append(f'which has a fuzz value of {common.fuzz} and a count of {common.n}')
	others = dat.query('fuzz >= @common.fuzz').\
		groupby(['segment', 'mapped_start', 'mapped_stop', 'fuzz']).\
		size().\
		reset_index().\
		sort_values(by=0).\
		rename(columns={0:'n'})
	output.append('other DVGs with a fuzz >= most common DVG fuzz:')
	for idx, i in others.iterrows():
		output.append(f'{i.segment} {i.mapped_start - args.padSize}_{i.mapped_stop - args.padSize}, fuzz={i.fuzz}, n={i.n}')

	with open('figures/final/figure_s4.txt', 'w') as fp:
		for i in output:
			fp.write(i+'\n')


	rng = np.random.default_rng(seed=111)
	fig, axs = plt.subplots(2,2,figsize=(6.4*2, 4.8*2), constrained_layout=True)
	axs[0,0].bar(fuzz_sum.index, fuzz_sum, facecolor='#eaeaea', edgecolor='#333333', zorder=3, label='empirical')
	axs[0,0].set_ylabel('DVG species', size=18)
	axs[0,0].set_title('all segments', size=18)
	axs[0,0].set_xlim(-0.75, fuzz_sum.index.max()+0.75)
	axs[1,0].scatter(ufuzz.fuzz + rng.uniform(low = -0.25, high=0.25, size=ufuzz.shape[0]), 
		ufuzz['count'], alpha=0.5, facecolor='none', edgecolor='#333333')


	axs[1,0].set_ylabel('samples', size=18)
	axs[1,0].set_xlim(-0.75, fuzz_sum.index.max()+0.75)
	axs[1,0].set_xlabel('repeat length', size=18)

	axs[0,1].bar(pol_fuzz_sum.index, pol_fuzz_sum, facecolor='#eaeaea', edgecolor='#333333', zorder=3, label='empirical')
	axs[0,1].scatter(null_quantiles.index, null_quantiles['q0.5'], color='#333333', zorder=5, label='null')
	axs[0,1].set_xlim(-0.75, fuzz_sum.index.max()+0.75)
	for idx, i in null_quantiles.iterrows():
		axs[0,1].plot([idx,idx],i[['q0.025', 'q0.975']].values, color='#333333', zorder=4)



	axs[0,1].legend()
	axs[0,1].set_ylabel('\n', size=18)

	axs[0,1].set_title(r'polymerase, $\geq$500nt deletion', size=18)
	axs[1,1].scatter(pol_ufuzz.fuzz + rng.uniform(low = -0.25, high=0.25, size=pol_ufuzz.shape[0]), 
		pol_ufuzz['count'], alpha=0.5, facecolor='none', edgecolor='#333333')


	axs[1,1].set_ylabel('\n', size=18)
	axs[1,1].set_xlabel('repeat length', size=18)

	axs[1,1].set_xlim(-0.75, fuzz_sum.index.max()+0.75)

	x_pos = [-0.15, -0.17, -0.13, -0.1]
	for ax_idx, ax in enumerate(axs.flatten()):
		ax.text(x_pos[ax_idx], 1.0, 
			string.ascii_uppercase[ax_idx], color='#333333', 
			transform=ax.transAxes, size=16, fontweight='bold')

	fig.suptitle('Figure S4')
	fig.savefig('figures/final/figure_s4.pdf')
	plt.close()




if __name__ == "__main__":
    run()

