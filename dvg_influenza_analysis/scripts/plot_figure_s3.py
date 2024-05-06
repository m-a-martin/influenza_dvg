import argparse
import pandas as pd
import matplotlib.pyplot as plt
from collections import Counter
import string
import matplotlib.patches as patches
try: 
	from utils import plot_style
except:
	from scripts.utils import plot_style


def run():
	parser = argparse.ArgumentParser()
	# input files
	# NEE DTO ADD ALL RUNS AND LIB
	#parser.add_argument('--allRuns', default=None)
	parser.add_argument('--dat', default=None)
	#parser.add_argument('--mapDir', default=None)
	parser.add_argument('--padSize', default=None, type=int)
	args = parser.parse_args()
	#args.dat = 'output/parsed_dvgs_10_True_True_0_all.tsv'
	# todo infer from file??
	#args.padSize = 210
	#args.mapDir = "data/*_map.tsv"
	#args.readBreakDir = '*/Virema/*_read.txt'
	#args.seq = 'data/KJ942684.1.fasta'
	#args.gene='data/KJ942684.1.gff3'
	#args.highlightDVG = ['NS', '316', '545']
	#args.allRuns = 'data/all_runs.tsv'

	dat = pd.read_csv(args.dat, sep='\t')	
	# adjust start/stop locations
	# this helps us define the highlighted DVG easier
	dat['mapped_start'] -= args.padSize
	dat['mapped_stop'] -= args.padSize
	# group by segment, start, stop, and count number of occurences
	shared_dvg_counts = dat.groupby(['mapped_start', 'mapped_stop', 'segment']).size().reset_index()
	# now count
	shared_dvg_counted = Counter(shared_dvg_counts[0])
	# get most prominent DVG
	output = []
	output.append('AMONG ALL SEGMENTS:\n')
	output.append(f'the DVG present in the most samples is ')
	output[-1] += f'{shared_dvg_counts.sort_values(by=0, ascending=False).iloc[0].segment} '
	output[-1] += f'{shared_dvg_counts.sort_values(by=0, ascending=False).iloc[0].mapped_start} '
	output[-1] += f'{shared_dvg_counts.sort_values(by=0, ascending=False).iloc[0].mapped_stop}'
	output.append(f'which is present in ')
	output[-1] += f'{shared_dvg_counts.sort_values(by=0, ascending=False).iloc[0][0]} out of '
	output[-1] += f'{dat.Specid.unique().shape[0]} samples'
	output.append(f'the next most prevalent DVG is present in ')
	output[-1] += f'{shared_dvg_counts.sort_values(by=0, ascending=False).iloc[1][0]} samples'
	output[-1] += f'{shared_dvg_counts.sort_values(by=0, ascending=False).iloc[1].segment} '
	output[-1] += f'{shared_dvg_counts.sort_values(by=0, ascending=False).iloc[1].mapped_start} '
	output[-1] += f'{shared_dvg_counts.sort_values(by=0, ascending=False).iloc[1].mapped_stop}'
	output.append(f'{(shared_dvg_counts[0] == 1).sum()} out of ')
	output.append(f'{shared_dvg_counts.shape[0]} ')
	output[-1] += f'DVGs are present in only a single sample'
	output.append(f'{(shared_dvg_counts[0] == 1).sum()*100/dat.shape[0]}% of DVGs are observed in only 1 sample')
	output.append(f'DVGs present in mean [sd] {shared_dvg_counts[0].mean()} [{shared_dvg_counts[0].std()}] samples')


	# adjust start/stop locations
	# this helps us define the highlighted DVG easier
	dat['mapped_start'] -= args.padSize
	dat['mapped_stop'] -= args.padSize
	# group by segment, start, stop, and count number of occurences
	pol_shared_dvg_counts = dat.query('(segment == "PB2" | segment == "PB1" | segment == "PA") & (Stop - Start >= 500)').\
		groupby(['mapped_start', 'mapped_stop', 'segment']).size().reset_index()
	# now count
	pol_shared_dvg_counted = Counter(pol_shared_dvg_counts[0])

	output.append('AMONG POL SEGMENTS:\n')
	output.append(f'the DVG present in the most samples is ')
	output[-1] += f'{pol_shared_dvg_counts.sort_values(by=0, ascending=False).iloc[0].segment} '
	output[-1] += f'{pol_shared_dvg_counts.sort_values(by=0, ascending=False).iloc[0].mapped_start} '
	output[-1] += f'{pol_shared_dvg_counts.sort_values(by=0, ascending=False).iloc[0].mapped_stop}'
	output.append(f'which is present in ')
	output[-1] += f'{pol_shared_dvg_counts.sort_values(by=0, ascending=False).iloc[0][0]} out of '
	n = dat.query('segment == "PB2" | segment == "PB1" | segment == "PA"').Specid.unique().shape[0]
	output[-1] += f'{n} samples'
	output.append(f'the next most prevalent DVG is present in ')
	output[-1] += f'{pol_shared_dvg_counts.sort_values(by=0, ascending=False).iloc[1][0]} samples'
	output[-1] += f'{pol_shared_dvg_counts.sort_values(by=0, ascending=False).iloc[1].segment} '
	output[-1] += f'{pol_shared_dvg_counts.sort_values(by=0, ascending=False).iloc[1].mapped_start} '
	output[-1] += f'{pol_shared_dvg_counts.sort_values(by=0, ascending=False).iloc[1].mapped_stop}'
	output.append(f'{(pol_shared_dvg_counts[0] == 1).sum()} out of ')
	output.append(f'{pol_shared_dvg_counts.shape[0]} ')
	output[-1] += f'DVGs are present in only a single sample'
	output.append(f'{(pol_shared_dvg_counts[0] == 1).sum()*100/dat.shape[0]}% of DVGs are observed in only 1 sample')
	output.append(f'DVGs present in mean [sd] {pol_shared_dvg_counts[0].mean()} [{pol_shared_dvg_counts[0].std()}] samples')


	with open('figures/final/figure_s3.txt', 'w') as fp:
		for line in output:
			fp.write(line+'\n')


	plot_style()
	#fig, axs = plt.subplots(2,2,figsize=(6.4*2, 4.8*2), constrained_layout=True)
	fig, axs = plt.subplots(1,2, figsize=(6.4*2, 4.8), constrained_layout=True)
	axs[0].bar(shared_dvg_counted.keys(), 
		shared_dvg_counted.values(), facecolor='#eaeaea', edgecolor='#333333', zorder=2)
	axs[0].set_xlabel('number of samples', size=16)
	axs[0].set_ylabel('number of unique DVG species', size=16)
	axs[0].grid(axis='both', color='#eaeaea', zorder=0)
	axs[0].set_yscale('log')
	axs[0].set_title('all segments')
	axs[1].bar(pol_shared_dvg_counted.keys(), 
		pol_shared_dvg_counted.values(), facecolor='#eaeaea', edgecolor='#333333', zorder=2)
	axs[1].set_xlabel('number of samples', size=16)
	axs[1].set_ylabel('number of unique DVG species', size=16)
	axs[1].grid(axis='both', color='#eaeaea', zorder=0)
	axs[1].set_yscale('log')
	axs[1].set_title(r'PB2, PB1, PA, $\geq$500nt deletion')

	x_pos = [-0.13, -0.13]
	for ax_idx, ax in enumerate(axs.flatten()):
		ax.text(x_pos[ax_idx], 1.0, 
			string.ascii_uppercase[ax_idx], color='#333333', 
			transform=ax.transAxes, size=16, fontweight='bold')

	fig.suptitle('Figure S3')
	fig.savefig('figures/final/figure_s3.pdf')
	plt.close()





if __name__ == "__main__":
    run()

