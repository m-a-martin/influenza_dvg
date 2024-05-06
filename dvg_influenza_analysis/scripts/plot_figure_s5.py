import glob 
import pandas as pd
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
try: 
	from utils import plot_style
except:
	from scripts.utils import plot_style


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat', default=None)
	parser.add_argument('--allRuns', default=None)
	parser.add_argument('--mapDir', default=None)
	parser.add_argument('--padSize', default=None, type=int)
	args = parser.parse_args()
	segment_order = {'PB2': 0, 'PB1': 1, 'PA': 2, 'HA': 3, 
		'NP': 4, 'NA': 5, 'M': 6, 'NS': 7}
	#args.padSize = 210
	#args.dat = 'output/parsed_dvgs_10_True_0_PB2_PB1_PA.tsv' 
	#args.dat = 'output/parsed_dvgs_10_True_0_all.tsv' 

	#args.allRuns = 'data/all_runs.tsv'
	#args.libs = ['HK', 'perth', 'vic']
	#args.mapDir = "data/*_map.tsv"
	#args.padSize = 210

	# also get the maximum size for each segment
	segment_dat = \
		[pd.read_csv(i, sep='\t', index_col=0) for i in glob.glob(args.mapDir)]
	segment_len = {i.index.name: max(i.index)-(args.padSize*2) for i in segment_dat}

	# read in and process data
	dat = pd.read_csv(args.dat, sep='\t')
	dat['segment'] = dat['segment'].fillna("NA")
	# adjust positions to account for pad
	dat['mapped_start'] -= args.padSize
	dat['mapped_stop'] -= args.padSize

	all_runs = pd.read_csv(args.allRuns, sep='\t')
	all_runs['lib'] = all_runs['lib'].apply(lambda k: ''.join([i for i in k.split('_')[0] if not i.isdigit()]))
	all_runs = all_runs[(all_runs['type'] == 'clinical')]
	all_runs['Collection_date'] = pd.to_datetime(all_runs['Collection_date'])
	all_runs = all_runs[['SPECID', 'Collection_date']].drop_duplicates().\
		sort_values(by="Collection_date").reset_index(drop=True)
	specid_sort_dict = {i['SPECID']: idx for idx, i in all_runs.iterrows()}
	n_per_page = 10
	specid_page = {key: 
		np.floor(value/n_per_page).astype(int) for key, value in specid_sort_dict.items()}
	specid_page_order = {key: value%n_per_page for key, value in specid_sort_dict.items()}
	# need to sort data by specid
	
	plot_style()
	n_segs = dat['segment'].unique().shape[0]
	# group by specid and for each plot junction
	with PdfPages('figures/final/figure_s5.pdf') as pdf:
		for page in set(specid_page.values()):
			page_dat = dat[dat['Specid'].map(specid_page) == page]
			#fig, axs = plt.subplots(n_per_page, n_segs,
			#		figsize=(6.4*1.5, 4.8*4), constrained_layout=True)
			fig = plt.figure(figsize=(6.4*1.75, 4.8*4), constrained_layout=True)
			fig.suptitle(f'Figure S5 page {page+1}')
			subfigs = fig.subfigures(nrows=n_per_page, ncols=1)
			for page_order in range(n_per_page):
				if (page)*n_per_page + page_order < all_runs.shape[0]:
					group_dat = page_dat[page_dat['Specid'].map(specid_page_order) == page_order]
					specid = [key for key, val in specid_page_order.items() if 
						val == page_order and specid_page[key] == page][0]
					# now also group by segment
					subfig = subfigs[page_order]
					# blank line for spacing
					subfig.suptitle(f'\n{specid}')
					axs = subfig.subplots(nrows=1, ncols=3)
					for row_idx, row in group_dat.iterrows():
						_ = axs[segment_order[row['segment']]].plot(
							[row['mapped_start'], row['mapped_stop']], 
							[1,0], color='#333333', lw=0.5)
					for ax_idx, ax in enumerate(axs):
						_ = [ax.spines[i].set_visible(False) for i in ['left', 'right']]
						_ = ax.set_xlabel('position (nt)', size=mpl.rcParams['axes.titlesize'])
						_ = ax.set_title(f'{[key for key, val in segment_order.items() if val == ax_idx][0]}')
						_ = ax.set_yticks([])
			pdf.savefig(fig)
			plt.close()

					

if __name__ == "__main__":
    run()







