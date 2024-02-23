import pandas as pd
import argparse
import glob
import numpy as np
import string
from scipy.stats import kruskal
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
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
	bins = [750,1250]
	bin_labels = [f'[0,{bins[0]})',
		f'[{bins[0]},{bins[1]})',
		f'[{bins[1]},)']
	dat = dat.\
		assign(
			days_post_onset = lambda k: k['Specid'].map(day_dict),
			seg_len = lambda k: k['segment'].map(segment_len),
			dvg_len = lambda k: k['seg_len'] - (k.mapped_stop - k.mapped_start - 1),
			log10_rel_support = lambda k: np.log10(k.Rel_support),
			dvg_len_bin = lambda k: np.digitize(k.dvg_len, bins=bins))
	
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
		pivot(index=['enrollid', 'segment', 'seg_len', 'dvg_len', 'dvg_len_bin', 'mapped_start', 'mapped_stop', 't0', 't1'], 
			columns='within_host_order', values='Rel_support').reset_index().\
		fillna(0).\
		assign(
			tspan = lambda k: k.t1 - k.t0,
			init = lambda k: k[0] > 0,
			persistent = lambda k: (k[0] > 0) & (k[1] > 0))

	p_persistence = long_dat.query('init == True').\
		groupby(['tspan', 'dvg_len_bin']).\
		agg({'enrollid': len, 'persistent': sum}).\
		reset_index().\
		assign(p = lambda k: k.persistent/k.enrollid)
	p_persistence_arr = p_persistence[['dvg_len_bin', 'tspan', 'p']].\
		pivot(index='dvg_len_bin', columns='tspan', values='p')

		
	max_val = 1.025*dat.Rel_support.max()
	min_val = -.025*dat.Rel_support.max()

	output = []
	output.append(f'bin edges = {bins}')
	for seg in ['PB2', 'PB1', 'PA']:
		t = kruskal(*[i.Rel_support for idx,i in dat.query('segment == @seg').groupby('dvg_len_bin')])
		output.append('Kruskal-Wallis H-test  ')
		output[-1] += 'comparing relative DVG read support of DVGs in each bin on the'
		output[-1] += f' {seg} segment'
		output[-1] += f': statistic={t.statistic}; pvalue={t.pvalue}' 

	# fit a model for the longitudinal data
	# first with binned read support
	# then with continuous
	logit_in = long_dat.\
		query('init == True')\
		[['enrollid', 'segment', 'mapped_start', 'mapped_stop', 'dvg_len', 
			'dvg_len_bin', 'tspan', 0, 1, 'persistent']].\
		assign(log10_rel_0_support = lambda k: np.log10(k[0]),
			shared = lambda k: 1*k.persistent)
	
	binned_fit = smf.logit("shared ~ log10_rel_0_support + C(tspan) + C(dvg_len_bin)", 
			data=logit_in).fit()
	output.append(str(binned_fit.summary()))
	continuous_fit = smf.logit("shared ~ log10_rel_0_support + C(tspan) + dvg_len", 
			data=logit_in).fit()
	output.append(str(continuous_fit.summary()))
	
	with open('figures/final/figure_s8.txt', 'w') as fp:
		for line in output:
			fp.write(line + '\n')


	plot_style()
	cmap = mpl.colors.LinearSegmentedColormap.from_list("", ['white', '#333333'])
	output = []
	fig, axs_arr = plt.subplots(2,2,figsize=(6.4*2, 4.8*2))
	axs = axs_arr.flatten()
	for j, seg in enumerate(['PB2', 'PB1', 'PA']):
		axs[j] = jitter_boxplot(ax=axs[j], 
			d=[i.Rel_support for idx,i in dat.query('segment == @seg').groupby('dvg_len_bin')],
			i=dat.dvg_len_bin.drop_duplicates().values,
			j=0.25, vert=True, zorder=3, point_alpha=0.5)
		#axs[j].set_yscale('log')
		axs[j].set_xticks(dat.dvg_len_bin.drop_duplicates().values)
		axs[j].set_xticklabels(bin_labels)
		axs[j].set_ylim(min_val, max_val)
		axs[j].grid(color='#eaeaea', zorder=1)
		axs[j].set_ylabel('relative support', size=16)
		axs[j].set_xlabel('DVG length (nt)', size=16)
		axs[j].set_title(seg)
		axs[j].text(-0.15, 1, 
				string.ascii_uppercase[j], color='#333333', 
				transform=axs[j].transAxes, size=18, fontweight='bold')

	im = axs[3].pcolormesh(p_persistence_arr.T, cmap=cmap, edgecolor='#333333',
		vmin=0, vmax=1)
	axs[3].set_xticks(np.arange(0.5, p_persistence_arr.shape[0]+0.5),
		labels = bin_labels)
	axs[3].set_yticks(np.arange(0.5, p_persistence_arr.shape[1]+0.5),
		labels=p_persistence_arr.columns.astype(int))
	axs[3].set_xlabel('DVG length (nt)', size=16)
	axs[3].set_ylabel('days between samples', size=16)
	divider = make_axes_locatable(axs[3])
	cax = divider.append_axes('right', size='5%', pad=0.1)
	cb = fig.colorbar(im, cax=cax, orientation='vertical')
	cb.set_label('prop. persistent DVGs', size=16)
	axs[3].text(-0.1, 1, 
				string.ascii_uppercase[3], color='#333333', 
				transform=axs[3].transAxes, size=18, fontweight='bold')
	fig.suptitle('Figure S8')
	fig.tight_layout()
	fig.savefig('figures/final/figure_s8.pdf')
	plt.close()


if __name__ == "__main__":
    run()
