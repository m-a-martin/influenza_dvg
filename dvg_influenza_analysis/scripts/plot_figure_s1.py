import matplotlib.pyplot
import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import glob
import string
try: 
	from utils import plot_style, jitter_boxplot
except:
	from scripts.utils import plot_style, jitter_boxplot
	

def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dataDir', default=None)
	parser.add_argument('--dpDir', default=None)
	parser.add_argument('--cutoff', default=5E-3, type=int)
	parser.add_argument('--segmentAcc', default='data/acc_segment.tsv')
	parser.add_argument('--allRuns', default=None)
	args = parser.parse_args()
	#args.allRuns = 'data/all_runs.tsv'
	#args.dataDir = '*_output/Virema/*.par'
	#args.dpDir = '*_output/depth/*_dp.tsv'
	
	all_runs = pd.read_csv(args.allRuns, sep='\t')
	plasmid_runs = all_runs[all_runs['type'] == 'plasmid']
	segment_acc = pd.read_csv(args.segmentAcc, sep='\t', header=None)
	segment_acc_dict = {i[0]: i[1] for i in segment_acc.values}
	# todo infer from file??
	segment_order = {'PB2': 0, 'PB1': 1, 'PA': 2, 'HA': 3, 
		'NP': 4, 'NA': 5, 'M': 6, 'NS': 7}
	#args.runTable = 'data/SraRunTable.txt'
	#args.dataDir = 'hk*_output/Virema/*.par'
	plasmid_dict = {i['lib']: i['Run'] for idx, i in all_runs.iterrows() if i['type'] == 'plasmid'}
	dir_files = glob.glob(args.dataDir)
	# create mapping of id to actually directory
	id_dir_dict = {i.split('/')[-1].split('_')[0]: i for i in dir_files}
	plasmid_dat = {i: 
		pd.read_csv(path, sep='\t').assign(run = i) for i,path in id_dir_dict.items() if i in plasmid_dict.values()}

	all_plasmid_dat = pd.concat(list(plasmid_dat.values()))
	all_plasmid_dat['segment'] = all_plasmid_dat['Segment'].map(segment_acc_dict)
	
	# add depth
	wt_dp_dict = {i.split('/')[-1].split('_')[0]: 
		{j[0]:j[1] for idx,j in pd.read_csv(i, sep='\t', header=None).\
					rename(columns={0: 'seg', 1: 'pos', 2: 'dp'}).\
					assign(mid_pos = lambda g: 
						g.groupby('seg').pos.transform(lambda k: int(k.max()/2))).\
					query('pos == mid_pos')[['seg', 'mid_pos']].iterrows()} 
			for i in glob.glob(args.dpDir)}
	all_plasmid_dat['Rel_support'] = all_plasmid_dat['Total_support']/\
		all_plasmid_dat.apply(lambda k: wt_dp_dict[k.run][k.Segment], axis=1)

	# number of DVGs in each seqment for each sequenced plasmid
	segment_dvg_count = all_plasmid_dat.groupby(['segment', 'run']).size().reindex(
		pd.MultiIndex.from_product(
			[segment_order.keys(), plasmid_runs['Run']])).fillna(0.0).reset_index()
	segment_dvg_count_dict = {idx: i[0] for idx, i in segment_dvg_count.groupby('level_0')}

	# stats output file
	output = []
	output.append('segment dvg count')
	output.extend([f'{key}: {np.median(val)} [{np.std(val)}]' for 
		key, val in segment_dvg_count_dict.items()])
	output.append(f'all: {np.median(np.hstack(list(segment_dvg_count_dict.values())))}')
	output[-1] += f' [{np.std(np.hstack(list(segment_dvg_count_dict.values())))}]'

	# group by segment and get read support for all DVGs observed in that segment
	# no need to fill in missing rows because this is 
	# dvg support given dvg observed
	segment_read_support_per_dvg = {idx: 
		i['Rel_support'] for idx, i in all_plasmid_dat.groupby('segment')}

	output.append('segment read support per dvg')
	output.extend([f'{key}: {np.median(val.values)} [{np.std(val.values)}' for 
		key, val in segment_read_support_per_dvg.items()])
	output.append(f'all: {np.median(pd.concat(segment_read_support_per_dvg.values()).values)}')
	output[-1] += f' [{np.std(pd.concat(segment_read_support_per_dvg.values()).values)}]'
	output.append(f'{(all_plasmid_dat["Rel_support"] < args.cutoff).sum()/all_plasmid_dat.shape[0]} of DVGs supported by <{args.cutoff} rel reads')
	with open('figures/final/figure_s1.txt', 'w') as fp:
		for line in output:
			fp.write(line + '\n')
	plot_style()

	fig, axs = plt.subplots(1,2, figsize=(6.4*2, 4.8), constrained_layout=True)
	axs[0] = jitter_boxplot(ax=axs[0], 
		d=list(segment_dvg_count_dict.values()),
		i=[segment_order[i] for i in segment_dvg_count_dict.keys()])
	axs[0].grid(axis='y', color='#eaeaea', zorder=0)
	axs[0].set_xticks(range(8))
	axs[0].set_xticklabels(segment_order.keys())
	axs[0].set_xlabel('segment', size=14)
	axs[0].set_ylabel('unique DVG species', size=mpl.rcParams['axes.titlesize'])
	axs[1] = jitter_boxplot(ax=axs[1], 
		d=list(segment_read_support_per_dvg.values()),
		i=[segment_order[i] for i in segment_read_support_per_dvg.keys()])
	#ax.boxplot(all_plasmid_dat['segment'].map(segment_order), 
		#all_plasmid_dat['Total_support'])
		#facecolor='steelblue', edgecolor='#333333',
		#alpha=0.75, zorder=2, s=150)
	axs[1].grid(axis='y', color='#eaeaea', zorder=0)
	axs[1].set_xticks(range(8))
	axs[1].set_xticklabels(segment_order.keys())
	axs[1].set_xlabel('segment', size=mpl.rcParams['axes.titlesize'])
	axs[1].set_ylabel('rel. supporting reads per DVG', size=mpl.rcParams['axes.titlesize'])
	axs[1].axhline(args.cutoff, color='#4d4d4d', ls='--', zorder=2)
	for ax_idx, ax in enumerate(axs):
		ax.text(-0.125, 1.0, 
			string.ascii_uppercase[ax_idx], color='#333333', 
			transform=ax.transAxes, size=16, fontweight='bold')

	fig.suptitle('Figure S1')
	plot_style()
	fig.savefig(f'figures/final/figure_s1.pdf')
	plt.close()


if __name__ == "__main__":
    run()