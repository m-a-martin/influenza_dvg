import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import string
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
	parser.add_argument('--mapDir', default=None)
	parser.add_argument('--minTspan', default=-np.inf, type=float)
	parser.add_argument('--padSize', default=None, type=int)
	args = parser.parse_args()
	segment_order = {'PB2': 0, 'PB1': 1, 'PA': 2, 'HA': 3, 
		'NP': 4, 'NA': 5, 'M': 6, 'NS': 7}
	#args.dat = 'output/parsed_dvgs_10_True_True_500_PB2_PB1_PA.tsv'
	#args.allRuns = "data/all_runs.tsv"
	#args.repEnrollID = '50538'
	# create an ordering of the specids, by sample date arbitrary
	# first, subset run table to just samples we have
	all_runs = pd.read_csv(args.allRuns, sep='\t')
	all_runs = all_runs[(all_runs['type'] == 'clinical')]
	all_runs['Collection_date'] = pd.to_datetime(all_runs['Collection_date'])
	all_runs = all_runs.sort_values(by='Collection_date')

	day_dict = {i['SPECID']:i['days_post_onset'] for 
		idx,i in all_runs.iterrows()}

	dat = pd.read_csv(args.dat, sep='\t')	

	# map days post symptom onset onto dat
	dat['days_post_onset'] = dat['Specid'].map(day_dict)

	# get longitudinal data
	# only get rows from longitudinal data
	keep_enrollids = dat[['Specid', 'enrollid', 'days_post_onset']].\
		groupby(['enrollid']).\
		filter(lambda k: 
			len(set(k.Specid)) > 1).\
		groupby(['enrollid']).\
		filter(lambda k:  
			(max(k.days_post_onset) - min(k.days_post_onset)) >= args.minTspan)\
		['enrollid'].drop_duplicates()
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
			lambda k: k[1] - k[0], axis=1).items()}
	long_dat['tspan'] = long_dat['enrollid'].map(tspan_dict)

	long_specid = all_runs[all_runs['Isolate'].astype(int).isin(
		long_dat['enrollid'])]['Isolate'].drop_duplicates().reset_index()

	specid_sort_dict = {i['Isolate']: idx for idx, i in long_specid.iterrows()}

	n_per_page = 5
	specid_page = {key: 
		np.floor(value/n_per_page).astype(int) for key, value in specid_sort_dict.items()}
	specid_page_order = {key: value%n_per_page for key, value in specid_sort_dict.items()}

	within_host_order_dict = {idx: i.astype(int).astype(str) for idx, i in long_dat[['enrollid', 
		'within_host_order', 'days_post_onset']].drop_duplicates().pivot(index='enrollid', 
			columns='within_host_order', values='days_post_onset').iterrows()}
	# pivot the long dat
	long_dat = long_dat.pivot(index=['enrollid', 'segment', 'Start', 'Stop'],
		columns=['within_host_order'],
		values='Rel_support').fillna(0).reset_index()
	output = []
	output.append(f'longitudinal data with a tspan of at least {args.minTspan} is available for ')
	output[-1] += f'{long_dat.enrollid.drop_duplicates().shape[0]} individuals'

	with open('figures/final/figure_s7.txt', 'w') as fp:
		for line in output:
			fp.write(line+'\n')

	plot_style()
	with PdfPages('figures/final/figure_s7.pdf') as pdf:
		for page in set(specid_page.values()):
			page_dat = long_dat[long_dat['enrollid'].astype(str).map(specid_page) == page]
			#fig, axs = plt.subplots(n_per_page, n_segs,
			#		figsize=(6.4*1.5, 4.8*4), constrained_layout=True)
			fig = plt.figure(figsize=(6.4*2, 4.8*3), constrained_layout=True)
			fig.suptitle(f'Figure S7 page {page+1}')
			subfigs = fig.subfigures(nrows=n_per_page, ncols=1)
			for page_order in range(n_per_page):
				if (page)*n_per_page + page_order < long_specid.shape[0]:
					group_dat = page_dat[page_dat['enrollid'].astype(str).\
						map(specid_page_order) == page_order]
					specid = [key for key, val in specid_page_order.items() if 
						val == page_order and specid_page[key] == page][0]
					# now also group by segment
					subfig = subfigs[page_order]
					# blank line for spacing
					subfig.suptitle(f'\n{specid}')
					axs = subfig.subplots(nrows=1, ncols=3)
					for row_idx, row in group_dat.iterrows():
						_ = axs[segment_order[row['segment']]].scatter([0,1], 
							[row[0], row[1]],
							facecolor='#eaeaea', edgecolor='#333333', zorder=4)
						_ = axs[segment_order[row['segment']]].plot([0,1], 
							[row[0], row[1]],
							color='#333333', zorder=2)
					for ax_idx, ax in enumerate(axs):
						_ = ax.set_xlabel('days post symptom onset', size=16)
						_ = ax.set_xticks([0,1])
						_ = ax.set_xticklabels(within_host_order_dict[int(specid)].values)
						_ = ax.set_title(f'{[key for key, val in segment_order.items() if val == ax_idx][0]}')
						#_ = ax.set_yticks([])
						ax.grid(color='#eaeaea',axis='y')
						ax.set_ylabel('rel. read support', size=16)
			pdf.savefig(fig)
			plt.close()




if __name__ == "__main__":
    run()
