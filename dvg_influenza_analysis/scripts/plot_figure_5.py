import glob 
import string
import pandas as pd
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
try: 
	from utils import plot_style
except:
	from scripts.utils import plot_style


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--allRuns', default=None)
	parser.add_argument('--dat', default=None)
	parser.add_argument('--pairDat', default=None)
	parser.add_argument('--doubleDat', default=None)
	args = parser.parse_args()
	#args.dat = 'output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv'
	#args.allRuns = "data/all_runs.tsv"
	#args.pairDat = 'data/elife-35962-fig3-data1-v3.csv'
	#args.doubleDat = 'data/transmission_pairs.csv'

	output = []
	pair_dat = pd.read_csv(args.pairDat)
	output.append(f'there are {pair_dat.Household.sum()} household pairs in the mccrone data\n')
	cutoff = pair_dat.\
		query('Household == False').\
		sort_values(by='L1_norm').\
		assign(
			n = lambda k: k.shape[0], 
			idx = lambda k: np.arange(k.shape[0])).\
		query('idx <= n*0.05')['L1_norm'].max()

	pair_dat = pair_dat.\
		assign(valid_household = lambda k: 
			k.Household & (k.L1_norm <= cutoff))

	output.append(f'of these, {pair_dat.valid_household.sum()} have a genetic distance below the L1-norm cutoff')

	# map specid to enrollid 
	# this serves two purposes
	# 1) allows us to look at all DVGs across longitudinal samples from a single house
	# 2) remove H1N1 pairs
	all_runs = pd.read_csv(args.allRuns, sep='\t').\
		query('type == "clinical"')
	enrollid_dict = {i['SPECID']: i['Isolate'] for idx, i in all_runs.iterrows()}

	pair_dat = pair_dat.\
		assign(
			enrollid_1 = lambda k: k.SPECID1.map(enrollid_dict),
			enrollid_2 = lambda k: k.SPECID2.map(enrollid_dict)).\
		query('~enrollid_1.isnull() & ~enrollid_2.isnull()')

	output.append(f'of these, {pair_dat.valid_household.sum()} are from H3N2 infections')

	# drop non-valid household pairs and add type column
	pair_dat = pair_dat.query('~Household | (Household & valid_household)').\
		assign(type = lambda k: (k.Household & k.valid_household).map({True: 'household', False: 'community'}))

	double_dat = pd.read_csv(args.doubleDat).\
		query('(snv_qualified_pair == True) & (valid == True) & (quality_distance == True)').\
		query('(ENROLLID1 in @enrollid_dict.values()) & (ENROLLID2 in @enrollid_dict.values())')\
		[['ENROLLID1', 'ENROLLID2', 'double']].\
		rename(columns={'ENROLLID1': 'enrollid_1', 'ENROLLID2': 'enrollid_2'})

	if double_dat.shape[0] != 39:
		raise Exception('double dat is wrong length')

	pair_dat = pair_dat.merge(double_dat, how='left', on=['enrollid_1', 'enrollid_2']).\
		assign(double = lambda k: np.where(pd.isnull(k.double), False, k.double))

	# read in dvg dat 
	dat = pd.read_csv(args.dat, sep='\t').assign(enrollid = lambda k: k.enrollid.astype(str))
	# which individuals have longitudinal data
	long = dat[['enrollid', 'Specid']].drop_duplicates().\
		groupby('enrollid').filter(lambda k: k.shape[0] > 1)['enrollid'].\
		drop_duplicates()
		
	# in the case of longitudinal samples, get maximum read support across the two
	dat = dat[['enrollid', 'Segment', 'mapped_start', 'mapped_stop', 'Rel_support']].\
		groupby(['enrollid', 'Segment', 'mapped_start', 'mapped_stop'])\
		['Rel_support'].\
		apply(max).reset_index()
	
	trans_pairs = pair_dat.query('type == "household"').\
		assign(enrollid_1_long = lambda k: np.isin(k.enrollid_1, long.values),
			enrollid_2_long = lambda k: np.isin(k.enrollid_2, long.values))

	
	#### TABULATE HOW MANY PAIRS HAVE LONGITUDINAL DATA ###
	long_pairs = pair_dat.\
		query('(type == "household") & ((enrollid_1 in @long.values) | (enrollid_2 in @long.values))')
	long_pair_ind = \
		pd.concat([pair_dat.query('(type == "household") & ((enrollid_1 in @long.values))').\
				enrollid_1.reset_index(drop=True),
			pair_dat.query('(type == "household") & ((enrollid_2 in @long.values))').\
				enrollid_2.reset_index(drop=True)]).drop_duplicates()

	output.append(f'a total of {long_pair_ind.shape[0]} individuals across')
	output[-1] += f' {long_pairs.shape[0]} transmission pairs have longitudinal data'
	
	#### TABULATE HOW MANY PAIRS HAVE AMBIGUOUS ORDERING OF TRANSMISSION ####
	output.append(f'the order of transmission is ambiguous for ')
	output[-1] += f'''{pair_dat.query('type == "household"').double.sum()} transmission pairs'''
	
	all_pair_dat = []
	np.random.seed(seed=111)
	for pair_idx, pair in pair_dat.iterrows():
		# when random community pair
		# shuffle pair so it's random
		# no true donor or recipient
		pair_enrollids = pair[['enrollid_1', 'enrollid_2']].values
		if pair.type == 'community':
			np.random.shuffle(pair_enrollids)
		pair_order = {i:idx for idx, i in enumerate(pair_enrollids)}
		p_dat = dat.query('enrollid in @pair_enrollids')
		p_dat = p_dat.assign(
			pair = pair_idx,
			pair_order = lambda k: k.enrollid.map(pair_order))
		p_dat = p_dat.pivot(index=['pair', 'Segment', 'mapped_start', 'mapped_stop'], 
			columns='pair_order', values='Rel_support').reset_index().fillna(0).\
			rename(columns={0: 'rel_support_0', 1: 'rel_support_1'})
		# if not dvgs in one of the samples
		p_dat = p_dat.assign(
			rel_support_0 = 
				lambda k: k.rel_support_0 if 'rel_support_0' in k.columns else 0,
			rel_support_1 = 
				lambda k: k.rel_support_1 if 'rel_support_1' in k.columns else 0,
			subject_0 = pair_enrollids[0],
			subject_1 = pair_enrollids[1],
			type = pair.type,
			double = pair.double,
			cat = lambda k: 
				((k.rel_support_0 > 0).astype(str) + (k.rel_support_1 > 0).astype(str)).map(
					{'TrueTrue': 'shared',
					 'TrueFalse': 'donor_only',
					 'FalseTrue': 'recipient_only'}))
		all_pair_dat.append(p_dat)

	all_pair_dat = pd.concat(all_pair_dat)

	n_shared = all_pair_dat.query('cat == "shared" & type == "household"').shape[0]
	n_shared_pairs = all_pair_dat.query('cat == "shared" & type == "household"').pair.drop_duplicates().shape[0]

	all_pair_counted = \
		all_pair_dat.groupby(['type', 'pair']).\
			apply(lambda k: (k.cat == 'shared').sum()).\
			reset_index().\
			rename(columns={0:'num_shared'})

	all_pair_counted_plot = \
		all_pair_counted.\
			groupby(['type', 'num_shared']).\
			size().reset_index().\
			rename(columns={0:'count'}).\
			assign(
				num_shared_x = lambda k: np.where(k['num_shared'] < 10, k['num_shared'], 10),
				num_shared_label = lambda k: np.where(k['num_shared'] < 10, k['num_shared'].astype(str), "10+")).\
			groupby(['type', 'num_shared_x', 'num_shared_label'])['count'].sum().\
			reset_index().\
			rename(columns={0:'count'}).\
			assign(
				total = lambda k: k.groupby('type')['count'].transform('sum'),
				density = lambda k: k['count'] / k['total'])

	# shared between donor/recipient
	# need to duplicate data for pairs with ambiguous order of transmission 
	pair_shared_dat = \
		pd.concat([
			all_pair_dat.query("type == 'household' & cat == 'shared'")['rel_support_0'],
			all_pair_dat.query("type == 'household' & cat == 'shared' & double == True")['rel_support_1']])

	# present only in donor
	# need to duplicate data for pairs with ambiguous order of transmission 
	pair_donor_only_dat = \
		pd.concat([
			all_pair_dat.query("type == 'household' & cat =='donor_only'").\
				query('rel_support_0 > 0').rel_support_0,
			all_pair_dat.query("type == 'household' & cat =='recipient_only' & double == True").\
				query('rel_support_1 > 0').rel_support_1])


	#### TABULATE NUMBER OF SHARED BETWEEN PAIRS / NOT PAIRS ####
	# for pairs with ambiguous ordering, we consider them only *once* in this analysis
	#  mann whitney u 
	t = mannwhitneyu(all_pair_counted.query('type == "household"').num_shared,
		all_pair_counted.query('type == "community"').num_shared)
	output.append(f'among community samples there are an average of ')
	output[-1] += str(all_pair_counted.query("type == 'community'")['num_shared'].mean())
	output[-1] += ' ['
	output[-1] += str(all_pair_counted.query("type == 'community'")['num_shared'].std())
	output[-1] += '] shared DVGs'

	output.append(f'among household samples there are an average of ')
	output[-1] += str(all_pair_counted.query("type == 'household'")['num_shared'].mean())
	output[-1] += ' ['
	output[-1] += str(all_pair_counted.query("type == 'household'")['num_shared'].std())
	output[-1] += '] shared DVGs'

	output.append(f'Mann-Whitney U test comparing the number of shared DVGs among community and household samples ')
	output[-1] += f'statistic = {t.statistic}, p-value = {t.pvalue}'

	output.append(f'In total {n_shared} DVGs from {n_shared_pairs} transmission pairs are shared between donor and recipient')

	#### TABULATE READ SUPPORT OF PAIR SHARED/UNSHARED ####
	output.append(f'shared DVGs have mean [sd] donor relative read support of {pair_shared_dat.mean()} [{pair_shared_dat.std()}]')
	output.append(f'non-shared DVGs have mean [sd] donor relative read support of {pair_donor_only_dat.mean()} [{pair_donor_only_dat.std()}]')
	

	all_pair_dat.to_csv('pair_dat_test.tsv', sep='\t')
	
	mwu = mannwhitneyu(pair_shared_dat, pair_donor_only_dat)

	output.append(f'mann-whitney u test p-value comparing the relative support of shared and non-shared pair DVGs = {mwu.pvalue}')


	with open('figures/final/figure_5.txt', 'w') as fp:
		for line in output:
			fp.write(line+'\n')

	# now plot
	plot_style()
	fig, axs = plt.subplots(1, 3, figsize=(6.4*2.4, 4.8), constrained_layout=True)
	axs[0].bar(
		all_pair_counted_plot.query('type == "community"')['num_shared_x']-0.2,
		all_pair_counted_plot.query('type == "community"')['density'],
		width=0.4,
		edgecolor='#333333',
		facecolor='#4d4d4d',
		label='unlinked pair', zorder=3)
	axs[0].bar(
		all_pair_counted_plot.query('type == "household"')['num_shared_x']+0.2,
		all_pair_counted_plot.query('type == "household"')['density'],
		width=0.4,
		edgecolor='#333333',
		facecolor='#eaeaea',
		label='transmission pair', zorder=3)
	axs[0].set_xticks(
		ticks=all_pair_counted_plot[['num_shared_x', 'num_shared_label']].drop_duplicates()['num_shared_x'],
		labels=all_pair_counted_plot[['num_shared_x', 'num_shared_label']].drop_duplicates()['num_shared_label'])
	axs[0].set_xlabel('number of shared DVGs', size=mpl.rcParams['xtick.labelsize'])
	axs[0].set_ylabel('density', size=mpl.rcParams['xtick.labelsize'])
	axs[0].legend()
	axs[0].grid(color='#eaeaea', zorder=0, axis='y')

	xmin = np.floor(min(np.log10(pair_shared_dat).min(),
		np.log10(pair_donor_only_dat).min()))
	xmax = np.ceil(max(np.log10(pair_shared_dat).max(),
		np.log10(pair_donor_only_dat).max()))

	bins = np.arange(xmin, xmax, 0.2)

	# duplicate data when order is unknown
	axs[1].hist(np.log10(pair_donor_only_dat),
		bins = bins,
		facecolor='#eaeaea', edgecolor='#333333', density=True, label='non-shared', zorder=3)
	axs[1].grid(color='#eaeaea', zorder=0)
	axs[1].set_title('non-shared DVGs')
	axs[1].set_xlabel(r'$log_{10}$(relative read support)', size=mpl.rcParams['axes.titlesize'])
	axs[2].hist(np.log10(pair_shared_dat),
		bins = bins,
		facecolor='#eaeaea', edgecolor='#333333', density=True, label='non-shared', zorder=3)
	axs[2].grid(color='#eaeaea', zorder=0)
	axs[2].set_title('shared DVGs')
	axs[2].set_xlabel(r'$log_{10}$(relative read support)', size=mpl.rcParams['axes.titlesize'])
	
	x_pos = [-0.14, -0.13, -0.12]
	for ax_idx, ax in enumerate(axs.flatten()):
		ax.text(x_pos[ax_idx], 1.0, 
			string.ascii_uppercase[ax_idx], color='#333333', 
			transform=ax.transAxes, size=16, fontweight='bold')

	fig.savefig('figures/final/figure_5.pdf')
	plt.close()

	with open('figures/final/figure_5.txt', 'w') as fp:
		for line in output:
			fp.write(line + '\n')


if __name__ == "__main__":
    run()