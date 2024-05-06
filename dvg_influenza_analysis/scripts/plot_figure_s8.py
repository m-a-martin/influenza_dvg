import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import statsmodels.formula.api as smf
try: 
	from utils import plot_style
except:
	from scripts.utils import plot_style


def run():
	parser = argparse.ArgumentParser()
	# input files
	# NEE DTO ADD ALL RUNS AND LIB
	parser.add_argument('--allRuns', default=None)
	parser.add_argument('--dat', default=None)
	args = parser.parse_args()
	#args.dat = 'output/parsed_dvgs_10_True_True_500_PB2_PB1_PA.tsv'
	#args.allRuns = "data/all_runs.tsv"

	# create an ordering of the specids, by sample date arbitrary
	# first, subset run table to just samples we have
	all_runs = pd.read_csv(args.allRuns, sep='\t')
	all_runs = all_runs[(all_runs['type'] == 'clinical')]

	day_dict = {i['SPECID']:i['days_post_onset'] for 
		idx,i in all_runs.iterrows()}

	dat = pd.read_csv(args.dat, sep='\t')	

	# map days post symptom onset onto dat
	dat['days_post_onset'] = dat['Specid'].map(day_dict)

	# get longitudinal data
	# only get rows from longitudinal data
	uniq_enrollids = \
		np.unique(dat[['Specid', 'enrollid']].drop_duplicates()['enrollid'].values, return_counts=True)

	keep_enrollids = uniq_enrollids[0][uniq_enrollids[1] > 1]
	# based on this file
	# we assume that "HS" specids are household sampls
	# and "MH" specids are clinic samples
	# therefore, when sampling date is the same, 
	# assume HS comes before MS
	dat['within_host_order'] = dat['Specid'].apply(lambda k: 0 if k[0:2] == 'HS' else 1)

	# get just longitudinal data
	long_dat = dat[dat['enrollid'].isin(keep_enrollids)]\
		[['enrollid', 'Specid', 'days_post_onset', 'within_host_order', 'segment', 'Start', 'Stop', 
			'Rel_support']].sort_values(by=['days_post_onset', 'within_host_order'])
	# add tspan for each row
	tspan = long_dat.groupby(['enrollid', 'days_post_onset'])['days_post_onset'].sum()
	tspan_dict = {idx:i for idx, i in 
		long_dat[['enrollid', 'within_host_order', 'days_post_onset']].drop_duplicates().pivot(
			index=['enrollid'], columns='within_host_order', 
			values='days_post_onset').apply(
			lambda k: k[1] - k[0], axis=1).iteritems()}
	long_dat['tspan'] = long_dat['enrollid'].map(tspan_dict)

	# pivot long dat
	# we want a dataframe 
	# that pivots DVGs based on their enrollid, subject, start, stop,
	# then support t0 and support t1
	grouped_dat = long_dat[['enrollid', 'within_host_order', 'tspan',
		'segment', 'Start', 'Stop', 'Rel_support']].pivot(
			index=['enrollid', 'tspan', 'segment', 'Start', 'Stop'], 
			columns='within_host_order', 
			values='Rel_support').fillna(0)
	grouped_t0_dat = grouped_dat[grouped_dat[0] > 0]
	grouped_t0_dat['shared'] = 1*(grouped_t0_dat[1] > 0)


	# multivariate regression
	grouped_t0_dat = grouped_t0_dat.reset_index().rename(
			columns={0: 'rel_0_support'}).assign(
			log10_rel_0_support = lambda k: np.log10(k.rel_0_support))
	t_fit = smf.logit("shared ~ log10_rel_0_support + C(tspan)", 
		data=grouped_t0_dat).fit()

	x_vals = np.linspace(grouped_t0_dat.log10_rel_0_support.min(), grouped_t0_dat.log10_rel_0_support.max(), num=50)
	t_span_predict = {tspan: t_fit.predict(pd.DataFrame(np.vstack([x_vals,
			np.repeat(tspan, 50)]).T, 
			columns=['log10_rel_0_support', 'tspan'])) for tspan in 
		np.sort(grouped_t0_dat['tspan'].unique())}

	output = []
	counts = grouped_t0_dat[['enrollid', 'tspan']].drop_duplicates().groupby('tspan').size()
	for row_idx, row in enumerate(counts):
		output.append(f'longitudinal data available for {row} individuals sampled {counts.index[row_idx]} day apart')
	output.append(str(t_fit.summary()))
	output.append(str(t_fit.predict(pd.DataFrame([[1E-3, 0], 
		[-3, 1],
		[-3, 2],
		[-3, 3],
		[-3, 4],
		[-3, 6],
		[-4, 3],
		[-3, 3],
		[-2, 3]], columns=['log10_rel_0_support', 'tspan']))))

	with open('figures/final/figure_s8.txt', 'w') as fp:
			for line in output:
				fp.write(line+'\n')

	n_cols =3
	n_rows = int(np.ceil(len(t_span_predict.keys())/n_cols))
	plot_style()
	fig, axs = plt.subplots(n_rows, n_cols, 
		figsize=(6.4*n_cols, 4.8*n_rows),
		constrained_layout=True)

	for idx, (key, val) in enumerate(t_span_predict.items()):
		key_dat = grouped_t0_dat[grouped_t0_dat['tspan'] == key]
		key_dat_dict = {t: t_dat['log10_rel_0_support'] for t, t_dat in key_dat.groupby('shared')}
		#axs.flatten()[idx].boxplot(key_dat_dict.values(), 
		#	positions=list(key_dat_dict.keys()), 
		#	vert=False, 
		#	flierprops=dict(markeredgecolor='#333333'),
		#			medianprops=dict(color='steelblue'),
		#			zorder=4)
		#axs.flatten()[idx] = jitter_boxplot(ax=axs.flatten()[idx],
		#	i=list(key_dat_dict.keys()),
		#	d=key_dat_dict.values(),
		#	j=0.25,
		#	vert=False)
		x = key_dat.log10_rel_0_support.values
		y = key_dat.shared.values
		#	np.repeat(1, list(grouped_t0_dat_tspan1_dict.values())[1].shape[0])])
		rng = np.random.default_rng(seed=111)
		y = y + rng.uniform(low = -0.25/2, high=0.25/2, size=y.shape[0])
		axs.flatten()[idx].scatter(x,y,edgecolor='#333333', facecolor="None", marker='o', zorder=4, alpha=0.5)
		axs.flatten()[idx].set_yticks([0,1])
		axs.flatten()[idx].set_yticklabels(['NP', 'P'])
		secax = axs.flatten()[idx].secondary_yaxis('right')
		secax.set_ylabel(r'$P($'+'persistent'+r'$)$', size=mpl.rcParams['axes.titlesize'])
		secax.set_yticks([0, 0.25, 0.5, 0.75, 1])
		axs.flatten()[idx].plot(x_vals, val, color='#333333', zorder=3)
		axs.flatten()[idx].set_xlabel('relative read support', size=mpl.rcParams['axes.titlesize'])
		axs.flatten()[idx].grid(color='#eaeaea', zorder=0)
		axs.flatten()[idx].set_title(f'\n{int(key)} day{"s" if int(key) != 1 else ""} between samples')


	for ax_idx, ax in enumerate(axs.flatten()):
				ax.text(-0.1, 1.0, 
					string.ascii_uppercase[ax_idx], color='#333333', 
					transform=ax.transAxes, size=16, fontweight='bold')
	fig.suptitle('Figure S8')
	fig.savefig('figures/final/figure_s8.pdf')
	plt.close()




if __name__ == "__main__":
    run()
