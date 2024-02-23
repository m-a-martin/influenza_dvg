import pandas as pd
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import string
try: 
	from utils import plot_style,  find_premature_stop
except:
	from scripts.utils import plot_style, find_premature_stop


def run():
	parser = argparse.ArgumentParser()
	# input files
	# NEE DTO ADD ALL RUNS AND LIB
	parser.add_argument('--dat', default=None)
	parser.add_argument('--fastaDir', default=None)
	parser.add_argument('--refDir', default='ref/*/*')
	parser.add_argument('--padSize', default=None, type=int)
	parser.add_argument('--nullPremature', default=None)
	#parser.add_argument('--repEnrollID', default=None, type=int)
	args = parser.parse_args()
	#args.dat = 'output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv'
	#args.mapDir = "data/*_map.tsv"
	#args.fastaDir = "*_output/consensus/*.fasta"
	#args.nullPremature = "output/null_premature_stop.tsv"
	#args.padSize = 210

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

	null_premature_stop = pd.read_csv(args.nullPremature, sep='\t')
	sum_null = null_premature_stop.assign(premature_support = lambda k: k.premature * k.Rel_support).\
		groupby(['Specid', 'Segment', 'rep']).agg({'premature_support': sum}).\
		groupby(['Specid', 'Segment']).agg({'premature_support': 
			[lambda k: k.quantile(0.025), lambda k: k.quantile(0.5), lambda k: k.quantile(0.975)]}).\
		reset_index().\
		rename(columns={'<lambda_0>': 'q2.5',
			'<lambda_1>': 'q5',
			'<lambda_2>': 'q97.5'})
	sum_null.columns = ['Specid', 'Segment', 'q2.5', 'q5', 'q97.5']

	# We have already removed the pad from the sample-wise reference genomes, 
	# so need to remove it here from the DVG coordinates
	dat['premature_stop'] = dat.apply(lambda k: find_premature_stop(k.Start - args.padSize, k.Stop-args.padSize, 
		cds_dict[k.Segment], ref_dict[k.Segment]), axis=1)

	premature_stop_dat = dat.assign(rel_premature_stop = lambda k: k.Rel_support * k.premature_stop).\
		groupby(['Specid', 'Segment', 'segment']).\
		agg({
			'Rel_support': sum,
			'premature_stop': sum,
			'rel_premature_stop': sum}).reset_index().\
		assign(p_premature_obs = lambda k: k.rel_premature_stop / k.Rel_support).\
		merge(sum_null, on=['Specid', 'Segment'], how='left')

	print(dat)
	output = []
	for seg in ['PB2', 'PB1', 'PA']:
		t = mannwhitneyu(
			dat.query('segment == @seg & premature_stop == True').Rel_support,
			dat.query('segment == @seg & premature_stop == False').Rel_support)
		output.append('mann whitney u test ')
		output[-1] += f'comparing relative DVG read support of {seg} DVGs with '
		output[-1] += f'premature stop codons to those without '
		output[-1] += f': statistic={t.statistic}; pvalue={t.pvalue}' 
	
	with open('figures/final/figure_s10.txt', 'w') as fp:
		for line in output:
			fp.write(line + '\n')

	segments = ['PB2', 'PB1', 'PA']
	plot_style()
	
	rng = np.random.default_rng(seed=111)
	max_val1 = premature_stop_dat[['q97.5', 'rel_premature_stop']].max().max()*1.025
	min_val1 = -premature_stop_dat[['q97.5', 'rel_premature_stop']].max().max()*0.025
	
	max_val2 = dat.Rel_support.max()*1.025
	min_val2 = -dat.Rel_support.max()*0.025

	fig, axs = plt.subplots(2,3,figsize=(6.4*3, 4.8*2), constrained_layout=True)
	for idx, seg in enumerate(segments):
		# top row
		# total relative support of DVGs with premature stops
		# expected v. observed
		axs[0,idx].axline((min_val1,min_val1), (max_val1, max_val1), ls='--', color='indianred')
		for jdx, j in premature_stop_dat.query('segment == @seg').iterrows():
			axs[0,idx].plot([j['q2.5'], j['q97.5']], [j.rel_premature_stop, j.rel_premature_stop], color='#c1c1c1', zorder=2)
		axs[0,idx].scatter(
			premature_stop_dat.query('segment == @seg')['q5'],
			premature_stop_dat.query('segment == @seg').rel_premature_stop,
			facecolor='#eaeaea', edgecolor='#333333', alpha=0.5, zorder=3)
		axs[0,idx].set_xlabel('expected relative support', size=16)
		if idx == 0:
			axs[0,idx].set_ylabel('observed relative support', size=16)
		else: 
			axs[0,idx].set_ylabel('', size=16)
		axs[0,idx].set_title(seg)
		axs[0,idx].grid(zorder=0, color='#eaeaea')
		#axs[0,idx].set_xlim(min_val1, max_val1)
		#axs[0,idx].set_ylim(min_val1, max_val1)
		#axs[0,idx].set_yscale('log')
		#axs[0,idx].set_xscale('log')
		# bottom row
		# support of individual DVGs w/ and w/o premature stops
		axs[1,idx].scatter(dat.query('segment == @seg').premature_stop + \
				rng.uniform(low = -0.25/2, high=0.25/2, size=dat.query('segment == @seg').shape[0]), 
			dat.query('segment == @seg').Rel_support, 
			edgecolor='#333333', facecolor='none',
			 marker='o', zorder=4, alpha=0.5)
		axs[1,idx].set_xlabel('premature stop', size=16)
		if idx == 0:
			axs[1,idx].set_ylabel('observed relative support', size=16)
		else:
			axs[1,idx].set_ylabel(' ', size=16)
		axs[1,idx].set_xticks([0,1])
		axs[1,idx].set_xticklabels(['False', 'True'])
		axs[1,idx].set_title('\n')
		axs[1,idx].grid(zorder=0, axis='y', color='#eaeaea')
		axs[1,idx].set_ylim(min_val2, max_val2)
		axs[1,idx].set_xlim(-0.5, 1.5)

	x_pos = [-0.12, -0.075, -0.075, -0.1, -0.075, -0.075]
	for ax_idx, ax in enumerate(axs.flatten()):
		ax.text(x_pos[ax_idx], 1.0, 
			string.ascii_uppercase[ax_idx], color='#333333', 
			transform=ax.transAxes, size=16, fontweight='bold')

	fig.suptitle('Figure S10')
	fig.savefig('figures/final/figure_s10.pdf')
	plt.close()




if __name__ == "__main__":
    run()
'''
for i,g in dat.groupby('Specid'):
	rng = np.random.default_rng(seed=111)
	fig, axs = plt.subplots(1,3,figsize=(6.4*3, 4.8), constrained_layout=True)
	axs[0].scatter(g.query('segment == "PB2"').premature_stop + \
			rng.uniform(low = -0.25/2, high=0.25/2, size=g.query('segment == "PB2"').shape[0]), 
		g.query('segment == "PB2"').Rel_support, 
		edgecolor='#333333', facecolor='none',
		 marker='o', zorder=4, alpha=0.5)
	axs[0].set_yscale('log')
	axs[0].set_ylabel('relative support')
	axs[0].set_xlabel('premature stop')
	axs[0].set_xticks([0,1])
	axs[0].set_xticklabels(['False', 'True'])
	axs[0].set_title('PB2')
	axs[1].scatter(g.query('segment == "PB1"').premature_stop + \
			rng.uniform(low = -0.25/2, high=0.25/2, size=g.query('segment == "PB1"').shape[0]), 
		g.query('segment == "PB1"').Rel_support, 
		edgecolor='#333333', facecolor='none',
		 marker='o', zorder=4, alpha=0.5)
	axs[1].set_yscale('log')
	axs[1].set_ylabel('relative support')
	axs[1].set_xlabel('premature stop')
	axs[1].set_xticks([0,1])
	axs[1].set_xticklabels(['False', 'True'])
	axs[1].set_title('PB1')
	axs[2].scatter(g.query('segment == "PA"').premature_stop + \
			rng.uniform(low = -0.25/2, high=0.25/2, size=g.query('segment == "PA"').shape[0]), 
		g.query('segment == "PA"').Rel_support, 
		edgecolor='#333333', facecolor='none',
		 marker='o', zorder=4, alpha=0.5)
	axs[2].set_yscale('log')
	axs[2].set_ylabel('relative support')
	axs[2].set_xlabel('premature stop')
	axs[2].set_xticks([0,1])
	axs[2].set_xticklabels(['False', 'True'])
	axs[2].set_title('PA')
	fig.savefig(f'figures/stops/{i}.pdf')
	plt.close()
'''