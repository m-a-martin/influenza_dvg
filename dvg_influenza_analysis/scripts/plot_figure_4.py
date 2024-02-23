import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import string
from scipy.stats import kruskal
try: 
	from utils import plot_style, jitter_boxplot
except:
	from scripts.utils import plot_style, jitter_boxplot


def calc_h(dat):
	if dat.shape[0] == 1:
		return(0.0)
	else:
		pi = dat['Rel_support']/dat['Rel_support'].sum()
		return(-1*(pi * np.log(pi)).sum())


def calc_j(dat):
	if dat.shape[0] == 1:
		return(1)
	else:
		pi = dat['Rel_support']/dat['Rel_support'].sum()
		h = -1*(pi * np.log(pi)).sum()
		return(h / np.log(pi.shape[0]))


def run():
	parser = argparse.ArgumentParser()
	# input files
	# NEE DTO ADD ALL RUNS AND LIB
	parser.add_argument('--allRuns', default=None)
	parser.add_argument('--dat', default=None)
	args = parser.parse_args()
	#args.dat = 'output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv'
	#args.allRuns = "data/all_runs.tsv"
	#args.repEnrollID = '50538'
	
	# create an ordering of the specids, by sample date arbitrary
	# first, subset run table to just samples we have
	all_runs = pd.read_csv(args.allRuns, sep='\t')
	all_runs = all_runs[(all_runs['type'] == 'clinical')]
	
	day_dict = {i['SPECID']:i['days_post_onset'] for 
		idx,i in all_runs.iterrows()}

	dat = pd.read_csv(args.dat, sep='\t')	

	uniq_enrollids = \
		np.unique(dat[['Specid', 'enrollid']].drop_duplicates()['enrollid'].values, return_counts=True)
	keep_enrollids = uniq_enrollids[0][uniq_enrollids[1] > 1]
	
	# map days post symptom onset onto dat
	dat['days_post_onset'] = dat['Specid'].map(day_dict)

	# groupby and calculate shannon's diversity 
	h_dict = {idx: i[0].values for idx, i in 
		dat.groupby(['Specid', 'enrollid', 'days_post_onset']).\
			apply(calc_h).reset_index().groupby('days_post_onset')}

	# and evenness 
	j_dict = {idx: i[0].values for idx, i in 
		dat.groupby(['Specid', 'enrollid', 'days_post_onset']).\
			apply(calc_j).reset_index().groupby('days_post_onset')}

	output = []
	h_t = kruskal(*list(h_dict.values()))
	j_t = kruskal(*list(j_dict.values()))
	output.append('Kruskal-Wallis H-test  ')
	output[-1] += 'comparing diversity (H) of DVGs at each DPI '
	output[-1] += f': statistic={h_t.statistic}; pvalue={h_t.pvalue}' 
	output.append('Kruskal-Wallis H-test  ')
	output[-1] += 'comparing diversity (J) of DVGs at each DPI '
	output[-1] += f': statistic={j_t.statistic}; pvalue={j_t.pvalue}' 

	with open('figures/final/figure_4.txt', 'w') as fp:
		for line in output:
			fp.write(line + '\n')
	# get only DVGs present at both time points
	persistent_dat = dat.\
			merge(
				dat[['enrollid', 'Specid']].\
					drop_duplicates().\
					groupby('enrollid').\
					size().reset_index().\
					rename(columns={0:'size'}).\
					query('size > 1'),
				on='enrollid',
				how = 'inner').\
		assign(within_host_order=lambda k: 1*(k.Specid.str[0:2] == 'MH')).\
		groupby(['enrollid', 'Segment', 'segment', 'mapped_start', 'mapped_stop']).\
		filter(lambda k: (k.days_post_onset.unique().shape[0] > 1))

	persistent_h = persistent_dat.groupby(['Specid', 'enrollid', 'days_post_onset']).apply(calc_h).reset_index().\
		pivot(index='enrollid', columns='days_post_onset', values=0)

	seg_persistent_h = persistent_dat.\
		groupby(['Specid', 'enrollid', 'segment', 'days_post_onset', 'within_host_order']).\
		apply(calc_h).reset_index().\
		pivot(index=['enrollid','segment'], columns='days_post_onset', values=0).reset_index()

	persistent_j = persistent_dat.groupby(['Specid', 'enrollid', 'days_post_onset']).apply(calc_j).reset_index().\
		pivot(index='enrollid', columns='days_post_onset', values=0)

	seg_persistent_j = persistent_dat.\
		groupby(['Specid', 'enrollid', 'segment', 'days_post_onset', 'within_host_order']).\
		apply(calc_j).reset_index().\
		pivot(index=['enrollid','segment'], columns='days_post_onset', values=0).reset_index()


	plot_style()
	#fig = plt.figure(figsize=(6.4*3, 4.8*2)) 
	#gs = gridspec.GridSpec(2, 6)

	#axs = [plt.subplot(gs[0,0:2]),
	#	plt.subplot(gs[0,2:4]),
	#	plt.subplot(gs[0,4:6])]
	fig, axs_arr = plt.subplots(2,2, figsize=(6.4*2, 4.8*2), constrained_layout=True)
	axs = axs_arr.flatten()
	#axs[0].boxplot(h_dict.values(), 
	#		positions=[int(i) for i in list(j_dict.keys())],
	#		flierprops=dict(markeredgecolor='#333333'),
	#		medianprops=dict(color='steelblue'))
	axs[0] = jitter_boxplot(ax=axs[0], 
			d=h_dict.values(),
			i=[int(i) for i in list(h_dict.keys())],
			j=0.25, vert=True, zorder=3)
	axs[0].set_xlabel('\n', size=18)
	axs[0].set_ylabel('$H$ (diversity)', size=18)
	axs[0].grid(color='#eaeaea', zorder=0, axis='y')
	#axs[1].boxplot(j_dict.values(), 
	#		positions=[int(i) for i in list(j_dict.keys())],
	#		flierprops=dict(markeredgecolor='#333333'),
	#		medianprops=dict(color='steelblue'))
	axs[1] = jitter_boxplot(ax=axs[1], 
			d=j_dict.values(),
			i=[int(i) for i in list(j_dict.keys())],
			j=0.25, vert=True, zorder=3)
	axs[1].set_xlabel('\n', size=18)
	axs[1].set_ylabel('$J$ (evenness)', size=18)
	axs[1].grid(color='#eaeaea', zorder=0, axis='y')

	
	for idx, row in persistent_h.iterrows():
		axs[2].plot(row.dropna().index, row.dropna(), color='#333333', lw=0.5, zorder=3)
		axs[2].scatter(row.dropna().index, row.dropna(), facecolor='#eaeaea', edgecolor='#333333', zorder=4)


	axs[2].set_xlabel('days post symptom onset', size=18)
	axs[2].set_ylabel('$H$ (diversity)', size=18)
	axs[2].grid(color='#eaeaea', zorder=0, axis='y')
	axs[2].set_xticks(list(j_dict.keys()))

	for idx, row in persistent_j.iterrows():
		axs[3].plot(row.dropna().index, row.dropna(), color='#333333', lw=0.5, zorder=3)
		axs[3].scatter(row.dropna().index, row.dropna(), facecolor='#eaeaea', edgecolor='#333333', zorder=4)


	axs[3].set_xlabel('days post symptom onset', size=18)
	axs[3].set_ylabel('$J$ (evenness)', size=18)
	axs[3].grid(color='#eaeaea', zorder=0, axis='y')
	axs[3].set_xticks(list(j_dict.keys()))
	x_pos = [-0.1, -0.15, -0.125, -0.125]
	for ax_idx, ax in enumerate(axs):
			ax.text(x_pos[ax_idx], 1.0, 
				string.ascii_uppercase[ax_idx], color='#333333', 
				transform=ax.transAxes, size=16, fontweight='bold')

	fig.savefig('figures/final/figure_4.pdf')
	plt.close()

if __name__ == "__main__":
    run()


