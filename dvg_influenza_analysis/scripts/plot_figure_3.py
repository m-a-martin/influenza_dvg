import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import string
from scipy.stats import kruskal
from scipy.stats import binomtest
from scipy.stats.contingency import chi2_contingency
import statsmodels.formula.api as smf
try: 
	from utils import plot_style, filter_dvg_dat, format_map_dict, import_fasta, jitter_boxplot
except:
	from scripts.utils import plot_style, filter_dvg_dat, format_map_dict, import_fasta, jitter_boxplot


def run():
	parser = argparse.ArgumentParser()
	# input files
	# NEE DTO ADD ALL RUNS AND LIB
	parser.add_argument('--allRuns', default=None)
	parser.add_argument('--dat', default=None)
	parser.add_argument('--mapDir', default=None)
	parser.add_argument('--padSize', default=None, type=int)
	parser.add_argument('--repEnrollID', default=None, type=int)
	args = parser.parse_args()
	
	#args.dat = 'output/parsed_dvgs_10_True_True_500_PB2_PB1_PA.tsv'
	#args.allRuns = "data/all_runs.tsv"
	#args.repEnrollID = '50319'
	# create an ordering of the specids, by sample date arbitrary
	# first, subset run table to just samples we have
	all_runs = pd.read_csv(args.allRuns, sep='\t')
	all_runs = all_runs[(all_runs['type'] == 'clinical')]

	day_dict = {i['SPECID']:i['days_post_onset'] for 
		idx,i in all_runs.iterrows()}

	dat = pd.read_csv(args.dat, sep='\t')	

	# map days post symptom onset onto dat
	dat['days_post_onset'] = dat['Specid'].map(day_dict)

	# groupby and count
	day_dvg_count = dat.groupby(['Specid', 'days_post_onset']).size().reindex(
		pd.MultiIndex.from_frame(all_runs[['SPECID', 'days_post_onset']].\
			drop_duplicates())).fillna(0.0).reset_index()

	day_dvg_count_dict = {int(g): g_dat[0].values for g, g_dat in 
		day_dvg_count.groupby(['days_post_onset'])}

	day_dvg_count_p = kruskal(*list(day_dvg_count_dict.values()))
	output = []
	output.append('Kruskal-Wallis H-test for independent samples')
	output[-1] += ' of the number of unique DVGs per sample'
	output[-1] += ' grouped by the number of days post symptom onset'
	output[-1] += f' p value = {day_dvg_count_p.pvalue}'
	output.append('median [sd] number of unique DVGs per sample')
	output[-1] += ' for samples taken 6 days post onset = '
	output[-1] += f'{np.median(day_dvg_count_dict[6])} [{np.std(day_dvg_count_dict[6])}]'
	output.append('median [sd] number of unique DVGs per sample')
	output[-1] += ' for samples taken <6 days post onset = '
	output[-1] += f'{np.median(np.hstack([day_dvg_count_dict[i] for i in day_dvg_count_dict.keys() if i < 6]))}'
	output[-1] += f'[{np.std(np.hstack([day_dvg_count_dict[i] for i in day_dvg_count_dict.keys() if i < 6]))}]'

	# groupby and sum supporting reads 
	day_dvg_support = dat.groupby(['Specid', 'days_post_onset'])['Rel_support'].sum().reindex(
		pd.MultiIndex.from_frame(all_runs[['SPECID', 'days_post_onset']].\
			drop_duplicates())).fillna(0.0).reset_index()
	day_dvg_support_dict = {int(g): g_dat['Rel_support'].values for g, g_dat in 
		day_dvg_support.groupby(['days_post_onset'])}

	day_dvg_support_p = kruskal(*list(day_dvg_support_dict.values()))

	output.append('Kruskal-Wallis H-test for independent samples')
	output[-1] += ' of the total relative read support per sample'
	output[-1] += ' grouped by the number of days post symptom onset'
	output[-1] += f' p value = {day_dvg_support_p.pvalue}'
	output.append('median [sd] number of total relative DVG reads per sample')
	output[-1] += ' for samples taken 6 days post onset = '
	output[-1] += f'{np.median(day_dvg_support_dict[6])} [{np.std(day_dvg_support_dict[6])}]'
	output.append('median [sd] number of total relative DVG reads per sample')
	output[-1] += ' for samples taken <6 days post onset = '
	output[-1] += f'{np.median(np.hstack([day_dvg_support_dict[i] for i in day_dvg_support_dict.keys() if i < 6]))}'
	output[-1] += f'[{np.std(np.hstack([day_dvg_support_dict[i] for i in day_dvg_support_dict.keys() if i < 6]))}]'

	# get longitudinal data
	# only get rows from longitudinal data
	uniq_enrollids = \
		np.unique(dat[['Specid', 'enrollid']].drop_duplicates()['enrollid'].values, return_counts=True)
	keep_enrollids = uniq_enrollids[0][uniq_enrollids[1] > 1]
	# "HS" specids are household sampls
	# and "MH" specids are clinic samples
	# this can be verified in McCrone's data 
	# e.g. https://github.com/lauringlab/Host_level_IAV_evolution/blob/master/data/processed/secondary/meta_one.sequence.success.csv
	# therefore, when sampling date is the same, 
	# based on McCrone methods, 
	# household sample comes before clinic
	dat = dat.merge(
		dat[['enrollid', 'Specid', 'days_post_onset']].\
			drop_duplicates().\
			assign(clinic = lambda k: k.Specid.str[0:2] != 'HS').\
			sort_values(by=['clinic', 'days_post_onset']).\
			assign(within_host_order = lambda g: g.groupby('enrollid').cumcount()),
		how='left', on=['enrollid', 'Specid', 'days_post_onset'])
	# get just longitudinal data
	long_dat = dat[dat['enrollid'].isin(keep_enrollids)]\
		[['enrollid', 'Specid', 'days_post_onset', 'within_host_order', 'segment', 'Start', 'Stop', 
			'Rel_support']].sort_values(by=['days_post_onset', 'within_host_order'])
	print('NUMBER OF LONGITUDINAL')
	print(long_dat.enrollid.drop_duplicates().shape)
	# add tspan for each row
	tspan = long_dat.groupby(['enrollid', 'days_post_onset'])['days_post_onset'].sum()
	tspan_dict = {idx:i for idx, i in 
		long_dat[['enrollid', 'within_host_order', 'days_post_onset']].drop_duplicates().pivot(
			index=['enrollid'], columns='within_host_order', 
			values='days_post_onset').apply(
			lambda k: k[1] - k[0], axis=1).items()}
	long_dat['tspan'] = long_dat['enrollid'].map(tspan_dict)
	rep_long_dat = long_dat[long_dat['enrollid'] == int(args.repEnrollID)]
	rep_long_dat = rep_long_dat.pivot(index=['enrollid', 'segment', 'Start', 'Stop'], 
			columns=['Specid', 'days_post_onset', 'within_host_order'], 
			values='Rel_support').fillna(0)

	# pivot long dat
	# we want a dataframe 
	# that pivots DVGs based on their enrollid, subject, start, stop,
	# then support t0 and support t1
	grouped_dat = long_dat[['enrollid', 'within_host_order', 'tspan',
		'segment', 'Start', 'Stop', 'Rel_support']].pivot(
			index=['enrollid', 'tspan', 'segment', 'Start', 'Stop'], 
			columns='within_host_order', 
			values='Rel_support').fillna(0).reset_index()
	# get just DVGs present in the first sample
	grouped_t0_dat = grouped_dat[grouped_dat[0] > 0]
	# indicator for whether that DVG is shared
	grouped_t0_dat = grouped_t0_dat.\
		assign(shared = lambda k: 1*(k[1] > 0))
	# finally, just those that are seperated by 1 day
	grouped_t0_dat_tspan1 = grouped_t0_dat.query('tspan == 1')
	logit_in = grouped_t0_dat_tspan1.\
		rename(
			columns={
				0: 'log10_rel_0_support',
				1: 'log10_rel_1_support'}).\
		assign(
			log10_rel_0_support = lambda k: np.log10(k['log10_rel_0_support']),
			log10_rel_1_support = lambda k: np.log10(k['log10_rel_1_support']))
	print(logit_in)
	print(logit_in.enrollid.unique().shape)
	t_fit = smf.logit("shared ~ log10_rel_0_support", 
		data=logit_in).fit()
	output.append(f'logistic regression p value using t_0 relative read support')
	output[-1] += ' as a predictor for persistent/non-persistent DVGs'
	output[-1] += ' for DVGs sampled 1 day appart'
	output[-1] += f' {t_fit.pvalues["log10_rel_0_support"]}'
	output.append(str(t_fit.summary()))
	# group by individual then calculate proportion that are shared between time points
	t_span_dict = {g: g_dat['shared'] for 
		g, g_dat in grouped_t0_dat.reset_index().groupby('tspan')}
	t_span = grouped_t0_dat.groupby('tspan')['shared'].apply(lambda x: 
			pd.DataFrame([len(x)-sum(x), sum(x)])).reset_index().pivot(
			index='tspan', columns='level_1', values=0)

	t_span_p = chi2_contingency(t_span)[1]
	output.append(f'chi2 test of independence p value comparing counts of persistent/non-persistent')
	output[-1] += ' DVGs as a function of time between samples = '
	output[-1] += f' {t_span_p}'
	# now get proportion and confidence interval for each day
	t_span_prop = t_span.apply(lambda k: 
		(binomtest(k=k[1], n=k[1]+k[0],	 p=0.5).proportion_estimate, 
			binomtest(k=k[1], n=k[1]+k[0], p=0.5).proportion_ci()), axis=1)

	# want to fit model to this data
	# x = intercept, probability of existing one generation, time points
	def calc_sse(params, times, dat, model):
		times = np.array(times)
		predict_dat = np.apply_along_axis(model, 0, times, p=params)
		delta = predict_dat - np.array(dat)
		sse = (delta**2).sum()
		return(sse)

	persistent_model = lambda t, p: p[0]*((1-p[1])**(t*4))
	from scipy.optimize import minimize
	persistent_model_fit = minimize(calc_sse, (0.5, 0.5), 
		(list(t_span_prop.index), [i[0] for i in t_span_prop.values], persistent_model))

	output.append(f'model fit to the persistent of DVGs as a funciton of time between samples')
	output.append(f'with the functional form r((1-p)^t) has paramters (r,p) of ')
	output.append(f'{persistent_model_fit.x}')
	# sort all DVGs by relative read support, get 99th percentile
	rel_support_99th = np.log10(dat.sort_values(
		by='Rel_support').iloc[
			int(np.ceil(dat.shape[0]*0.99)),:]['Rel_support'])
	# what is the predicted probability of this DVG perssiting?
	prob_99th_persist = \
		t_fit.predict(pd.DataFrame([rel_support_99th], columns=['log10_rel_0_support']))

	output.append(f'the 99th percentile of observed relative read support is {rel_support_99th}')
	output[-1] += f' and it would be predicted to have a {prob_99th_persist} probability of persisting'


	with open('dvg_influenza_analysis/figures/final/figure_3.txt', 'w') as fp:
		for line in output:
			fp.write(line+'\n')

	plot_style()
	fig = plt.figure(figsize=(6.4*2.5, 4.8*2.5), constrained_layout=True)
	spec = fig.add_gridspec(ncols=6, nrows=2)
	axs = [fig.add_subplot(spec[0, :3]), fig.add_subplot(spec[0,3:]), 
		fig.add_subplot(spec[1,:2]), fig.add_subplot(spec[1,2:4]), 
		fig.add_subplot(spec[1,4:])]
	axs[0] = jitter_boxplot(ax=axs[0], 
		d=list(day_dvg_support_dict.values()),
		i=list(day_dvg_support_dict.keys()),
		j=0.5)
	axs[0].grid(axis='y', color='#eaeaea', zorder=1)
	axs[0].set_xlabel('days post symptom onset', size=mpl.rcParams['axes.titlesize'])
	axs[0].set_ylabel('total rel. DVG reads per sample', size=mpl.rcParams['axes.titlesize'])
	axs[0].set_yscale('log')
	axs[1] = jitter_boxplot(ax=axs[1], 
		d=list(day_dvg_count_dict.values()),
		i=list(day_dvg_count_dict.keys()),
		j=0.5)
	axs[1].grid(axis='y', color='#eaeaea', zorder=1)
	axs[1].set_xlabel('days post symptom onset', size=mpl.rcParams['axes.titlesize'])
	axs[1].set_ylabel('unique DVGs per sample', size=mpl.rcParams['axes.titlesize'])
	axs[1].set_yscale('log')
	axs[2].scatter([0]*rep_long_dat.shape[0],
		rep_long_dat.iloc[:,0], facecolor='#eaeaea', 
		edgecolor='#333333', zorder=4)
	axs[2].scatter([1]*rep_long_dat.shape[0],
		rep_long_dat.iloc[:,1], facecolor='#eaeaea', 
		edgecolor='#333333', zorder=4)
	for i,row in rep_long_dat.iterrows():
		axs[2].plot([i[-1] for i in row.index], 
			row, color='#333333', zorder=2)


	axs[2].set_xticks([0,1])
	axs[2].set_xticklabels([str(int(i[-2])) for i in rep_long_dat.columns])
	axs[2].set_xlabel('days post symptom onset', size=mpl.rcParams['axes.titlesize'])
	axs[2].set_ylabel('rel. read support', size=mpl.rcParams['axes.titlesize'])
	axs[2].grid(axis='y', color='#eaeaea')
	axs[2].set_xlim(-0.1, 1.1)
	axs[2].set_title(args.repEnrollID)
	
	grouped_t0_dat_tspan1_dict = \
		{1*g: g_dat[0] for g, g_dat in grouped_t0_dat_tspan1.groupby('shared')}
	#axs[4].boxplot([i.apply(np.log10) for i in grouped_t0_dat_dict.values()], 
	#	positions=list(grouped_t0_dat_dict.keys()), 
	#	vert=False, 
	#	flierprops=dict(markeredgecolor='#333333'),
	#			medianprops=dict(color='steelblue'), zorder=3)
	#axs[3] = jitter_boxplot(ax=axs[3], 
	#	d=[i.apply(np.log10) for i in grouped_t0_dat_tspan1_dict.values()],
	#	i=list(grouped_t0_dat_tspan1_dict.keys()),
	#	j=0.25,
	#	vert=False)
	x = pd.concat([i.apply(np.log10) for i in grouped_t0_dat_tspan1_dict.values()]).values
	y = np.hstack([np.repeat(0, list(grouped_t0_dat_tspan1_dict.values())[0].shape[0]),
		np.repeat(1, list(grouped_t0_dat_tspan1_dict.values())[1].shape[0])])
	rng = np.random.default_rng(seed=111)
	y = y + rng.uniform(low = -0.25/2, high=0.25/2, size=x.shape[0])
	axs[3].scatter(x,y,edgecolor='#333333', facecolor="None", marker='o', zorder=4, alpha=0.5)
	#axs[4].set_xscale('log')
	x_vals = pd.concat([logit_in.log10_rel_0_support, logit_in.log10_rel_1_support]).\
		reset_index().\
		rename(columns={
					0: 'x'}).\
		query('x > -inf').x
	x_vals = np.linspace(x_vals.min(), x_vals.max() , 1000)
	pred_y = t_fit.predict(pd.DataFrame(
		x_vals, 
		columns=['log10_rel_0_support']))
	axs[3].plot(x_vals, pred_y, color='#333333')
	print(t_fit.summary())
	print(t_fit.predict(pd.DataFrame(
		[-4.0, -2.0], 
		columns=['log10_rel_0_support'])))

	secax = axs[3].secondary_yaxis('right')
	secax.set_ylabel(r'$P($'+'persistent'+r'$)$', size=mpl.rcParams['axes.titlesize'])
	secax.set_yticks([0, 0.25, 0.5, 0.75, 1])
	axs[3].set_yticks([0,1])
	axs[3].set_yticklabels(['non persistent', 'persistent'])
	axs[3].grid(color='#eaeaea')
	axs[3].set_xlabel(r'$log_{10}$(relative read support)', size=mpl.rcParams['axes.titlesize'])
	print(t_span_prop)
	for idx, i in t_span_prop.items():
		_ = axs[4].scatter([idx], i[0], facecolor='#eaeaea', edgecolor='#333333', zorder=5)
		_ = axs[4].plot([idx, idx], i[1], color='#333333', zorder=4)
		_ = axs[4].plot([idx-0.1, idx+0.1], [i[1][1], i[1][1]], color='#333333', zorder=4)
		_ = axs[4].plot([idx-0.1, idx+0.1], [i[1][0], i[1][0]], color='#333333', zorder=4)

	x_vals = np.linspace(t_span_prop.index.min(), t_span_prop.index.max())
	y_vals = np.apply_along_axis(persistent_model, 0, x_vals, persistent_model_fit.x)
	axs[4].plot(x_vals, y_vals, ls='--', color='#4d4d4d')
	axs[4].grid(axis='y', color='#eaeaea', zorder=1)
	axs[4].set_xlabel('days between samples', size=mpl.rcParams['axes.titlesize'])
	axs[4].set_ylabel('proportion of persistent $t_0$ DVGs', size=mpl.rcParams['axes.titlesize'])
	axs[4].set_xticks([int(i) for i in t_span_prop.index])
	print(t_span_prop)


	x_vals=[-0.125, -0.125, -0.325, -0.4, -0.25]
	for ax_idx, ax in enumerate(axs):
			ax.text(x_vals[ax_idx], 1.0, 
				string.ascii_uppercase[ax_idx], color='#333333', 
				transform=ax.transAxes, size=16, fontweight='bold')

	fig.savefig('dvg_influenza_analysis/figures/final/figure_3.pdf')
	plt.close()




if __name__ == "__main__":
    run()


