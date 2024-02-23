import pandas as pd
import argparse
import glob
import numpy as np
import matplotlib.pyplot as plt
import string
from scipy.stats import mannwhitneyu
try: 
	from utils import plot_style, jitter_boxplot
except:
	from scripts.utils import plot_style, jitter_boxplot


def run():
	parser = argparse.ArgumentParser()
	# input files
	# NEE DTO ADD ALL RUNS AND LIB
	parser.add_argument('--dat', default=None)
	parser.add_argument('--allRuns', default=None)
	parser.add_argument('--mapDir', default=None)
	parser.add_argument('--padSize', default=None, type=int)
	args = parser.parse_args()
	#args.dat = 'output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv'
	#args.mapDir = "data/*_map.tsv"
	#args.padSize = 210
	#args.allRuns = "data/all_runs.tsv"
	all_runs = pd.read_csv(args.allRuns, sep='\t')
	# these maps include the padding
	segment_len = {
		pd.read_csv(i, sep='\t', index_col=0).index.name:
		pd.read_csv(i, sep='\t', index_col=0).index.max()-2*args.padSize for 
			i in glob.glob(args.mapDir)}

	# add days post symptom onset
	day_dict = {i['SPECID']:i['days_post_onset'] for 
			idx,i in all_runs.iterrows()}

	dat = pd.read_csv(args.dat, sep='\t')	

	# "HS" specids are household sampls
	# and "MH" specids are clinic samples
	# this can be verified in McCrone's data 
	# e.g. https://github.com/lauringlab/Host_level_IAV_evolution/blob/master/data/processed/secondary/meta_one.sequence.success.csv
	# therefore, when sampling date is the same, 
	# based on McCrone methods, 
	# household sample comes before clinic
	dat = dat.\
		assign(
			days_post_onset = lambda k: k['Specid'].map(day_dict),
			seg_len = lambda k: k['segment'].map(segment_len),
			dvg_len = lambda k: k['seg_len'] - (k.mapped_stop - k.mapped_start - 1),
			log10_rel_support = lambda k: np.log10(k.Rel_support))

	long_dat = dat.\
		merge(
			dat[['enrollid', 'Specid']].\
				drop_duplicates().\
				groupby('enrollid').\
				size().reset_index().\
				rename(columns={0:'size'}).\
				query('size > 1'),
			on='enrollid',
			how = 'inner').\
		assign(
			t0 = lambda k: 
					k.groupby('enrollid').days_post_onset.transform(lambda g: min(g)),
			t1 = lambda k: 
					k.groupby('enrollid').days_post_onset.transform(lambda g: max(g)),
			clinic = lambda k: k.Specid.str[0:2] != 'HS',
			within_host_order = lambda k:
				1*((k.days_post_onset == k.t1) & ((k.days_post_onset != k.t0) | (k.clinic == True)))).\
		pivot(index=['enrollid', 'segment', 'seg_len', 'dvg_len', 'mapped_start', 'mapped_stop', 't0', 't1'], 
			columns='within_host_order', values='Rel_support').reset_index().\
		fillna(0).\
		assign(
			init = lambda k: k[0] > 0,
			persistent = lambda k: (k[0] > 0) & (k[1] > 0))


	
	
	output = []
	output.append(f'there is sequence data taken 1 day apart for ')
	n = long_dat.query('(t1 - t0 == 1)').enrollid.drop_duplicates().shape[0]
	output[-1] += f'{n} individuals'
	for j, seg in enumerate(['PB2', 'PB1', 'PA']):	
		m_p = long_dat.query('segment == @seg & init & persistent & (t1 - t0 == 1)').dvg_len.mean()
		s_p = long_dat.query('segment == @seg & init & persistent & (t1 - t0 == 1)').dvg_len.std()
		m_n = long_dat.query('segment == @seg & init & ~persistent & (t1 - t0 == 1)').dvg_len.mean()
		s_n = long_dat.query('segment == @seg & init & ~persistent & (t1 - t0 == 1)').dvg_len.std()
		output.append(f'mean [sd] read support of persistent {seg} DVGs: ')
		output[-1] += f"{m_p} [{s_p}]"
		output.append(f'mean [sd] read support of non_persistent {seg} DVGs: ')
		output[-1] += f"{m_n} [{s_n}]"
		t = mannwhitneyu(long_dat.query('segment == @seg & init & ~persistent & (t1 - t0 == 1)').dvg_len,
				long_dat.query('segment == @seg & init & persistent & (t1 - t0 == 1)').dvg_len)
		output.append('mann whitney u test ')
		output[-1] += 'comparing relative DVG read support of DVGs <= 1000nt and > 1000nt on the'
		output[-1] += f' {seg} segment'
		output[-1] += f': statistic={t.statistic}; pvalue={t.pvalue}' 

	with open('figures/final/figure_3.txt', 'w') as fp:
		for line in output:
			fp.write(line + '\n')

	max_support = dat.Rel_support.max()*1.025
	min_support = -dat.Rel_support.max()*0.025

	max_count = 0
	plot_style()
	fig, axs = plt.subplots(3,3,figsize=(6.4*3, 4.8*2.5), constrained_layout = True)
	for j, seg in enumerate(['PB2', 'PB1', 'PA']):
		j_hist = axs[0,j].hist(dat.query('segment == @seg').dvg_len, facecolor='#eaeaea', edgecolor='#333333',
				bins = range(0, dat.query('segment == @seg').seg_len.values[0], 25),  zorder=3)
		max_count = max(max_count, j_hist[0].max())
		axs[0,j].set_xlim(0, dat.query('segment == @seg').seg_len.values[0])
		axs[0,j].set_title(seg)
		axs[0,j].set_ylabel('', size=16)
		axs[0,j].set_xlabel('', size=16)
		#axs[0,j].set_yticks([])
		_ = axs[0,j].axvline(dat.query('segment == @seg').seg_len.values[0] - 500, ls='--', color='#4d4d4d')
		_ = axs[0,j].grid(color='#eaeaea', zorder=1)
		axs[1,j].scatter(dat.query('segment == @seg').dvg_len,
			dat.query('segment == @seg').Rel_support, edgecolor='#333333', 
			alpha=0.25, facecolor='none',zorder=3)
		axs[1,j].set_xlim(0, dat.query('segment == @seg').seg_len.values[0])
		axs[1,j].set_ylim(min_support, max_support)
		axs[1,j].grid(color='#eaeaea', zorder=1)
		axs[1,j].set_ylabel('')
		axs[2,j] = jitter_boxplot(ax=axs[2,j], 
			d=[long_dat.query('segment == @seg & init & ~persistent & (t1 - t0 == 1)').dvg_len,
				long_dat.query('segment == @seg & init & persistent & (t1 - t0 == 1)').dvg_len],
			i=[0,1],
			j=0.25, vert=False, zorder=3)
		axs[2,j].set_xlim(0, dat.query('segment == @seg').seg_len.values[0])
		axs[2,j].set_xlabel('DVG length (nt)', size=18)
		axs[2,j].grid(color='#eaeaea', axis='x', zorder=1)
		axs[2,j].set_yticks([0,1])
		axs[2,j].set_yticklabels(['',''])

	for ax in axs[0,:]:
		ax.set_ylim(0, max_count*1.025)

	x_pos = [-0.15, -0.125, -0.125, -0.09, -0.075, -0.075, -0.275, -0.05, -0.05]
	for ax_idx, ax in enumerate(axs.flatten()):
		ax.text(x_pos[ax_idx], 1.0, 
			string.ascii_uppercase[ax_idx], color='#333333', 
			transform=ax.transAxes, size=16, fontweight='bold')
	axs[0,0].set_ylabel('DVG species', size=16)
	axs[1,0].set_ylabel(r'relative support', size=18)
	axs[2,0].set_yticklabels(['not persistent', 'persistent'])
	fig.savefig('figures/final/figure_3.pdf')
	plt.close()

	

if __name__ == "__main__":
    run()
