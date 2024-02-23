import matplotlib.pyplot
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import string
from scipy.stats import mannwhitneyu
from mpl_toolkits.axes_grid1 import make_axes_locatable
try: 
	from utils import plot_style, format_map_dict, jitter_boxplot
except:
	from scripts.utils import plot_style, format_map_dict, jitter_boxplot


def get_coding_region(f, map_dict):
		gene_dat = pd.read_csv(f, sep='\t', comment='#', header=None)
		# need to map all coordinates universal coordinate system
		gene_dat[3] = gene_dat[3].map(map_dict[gene_dat.iloc[0,0]])
		gene_dat[4] = gene_dat[4].map(map_dict[gene_dat.iloc[0,0]])
		# subset to just genes
		gene_dat = gene_dat[gene_dat[2] == 'gene']
		gene_cords = np.unique(np.hstack([np.arange(i[3], i[4]+1) for i in gene_dat.values]))
		return(gene_cords)


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat', default=None)
	parser.add_argument('--allRuns', default=None)
	parser.add_argument('--mapDir', default=None)
	parser.add_argument('--padSize', default=None, type=int)
	parser.add_argument('--repSpecid')
	#parser.add_argument('--pb2Gene')
	#parser.add_argument('--pb1Gene')
	#parser.add_argument('--paGene')
	args = parser.parse_args()
	#args.dat = 'output/parsed_dvgs_0.005_True_True_0_all.tsv'
	# todo infer from file??
	#args.padSize = 210
	#args.repSpecid = 'HS1530'
	#args.mapDir = "data/*_map.tsv"
	#args.allRuns = "data/all_runs.tsv"

	segment_dat = \
		[pd.read_csv(i, sep='\t', index_col=0) for i in glob.glob(args.mapDir)]
	segment_len = {i.index.name: max(i.index)-(args.padSize*2) for i in segment_dat}
	map_dict = format_map_dict(glob.glob(args.mapDir))
	
	segment_order = {'PB2': 0, 'PB1': 1, 'PA': 2, 'HA': 3, 
		'NP': 4, 'NA': 5, 'M': 6, 'NS': 7}
	dat = pd.read_csv(args.dat, sep='\t')
	dat['segment'] = dat['segment'].fillna('NA')

	
	# create an ordering of the specids, by sample date arbitrary
	# first, subset run table to just samples we have
	all_runs = pd.read_csv(args.allRuns, sep='\t')
	all_runs['lib'] = all_runs['lib'].apply(lambda k: ''.join([i for i in k.split('_')[0] if not i.isdigit()]))
	all_runs = all_runs[(all_runs['type'] == 'clinical')]
	all_runs['Collection_date'] = pd.to_datetime(all_runs['Collection_date'])
	all_runs = all_runs[['SPECID', 'Collection_date']].drop_duplicates().\
		sort_values(by="Collection_date").reset_index(drop=True)
	specid_sort_dict = {i['SPECID']: idx for idx, i in all_runs.iterrows()}
	
	# number of relative dvg reads per sample per segment
	# need to reset index to account for samples with no dvg reads
	sample_seg_rel_reads_df = \
		dat.groupby(['segment', 'Specid'])['Rel_support'].sum().reindex(
			pd.MultiIndex.from_product(
			[segment_order.keys(), all_runs['SPECID'].unique()], names=['segment', 'Specid'])).fillna(
		0.0).reset_index()
	sample_seg_rel_reads_df['x'] = \
		sample_seg_rel_reads_df['segment'].map(segment_order)
	sample_seg_rel_reads_df['y'] = \
		sample_seg_rel_reads_df['Specid'].map(specid_sort_dict)
	dat_arr = np.zeros((sample_seg_rel_reads_df['y'].max()+1, sample_seg_rel_reads_df['x'].max()+1))
	dat_arr[:] = np.nan
	dat_arr[sample_seg_rel_reads_df['y'], sample_seg_rel_reads_df['x']] = \
		np.log10(sample_seg_rel_reads_df['Rel_support'])
	# unique dvgs per segment per sample
	sample_seg_n_dvgs = dat.groupby(['segment', 'Specid']).size().reindex(
		pd.MultiIndex.from_product(
			[segment_order.keys(), all_runs['SPECID'].unique()], 
			names=['segment', 'Specid'])).fillna(
		0.0).reset_index()
	sample_seg_n_dvgs_dict = {idx: i[0] for idx, i in 
		sample_seg_n_dvgs[['segment', 0]].groupby('segment')}
	# how many specids do we observe dvgs in?
	output = []
	output.append(f'DVGs observed in {(sample_seg_n_dvgs[["Specid", 0]].groupby("Specid").sum() > 0).sum().values[0]}')
	output[-1] += f'/{sample_seg_n_dvgs["Specid"].unique().shape[0]} SPECIDS'

	# total relative reads per sample per segment
	sample_seg_rel_reads_dict = {idx: 
		i['Rel_support'].values for idx, i in sample_seg_rel_reads_df.groupby('segment')}


	# statistical test, DVG read support in PB2, PB1, PA, NS v. HA, NP, NA, M. 

	t = mannwhitneyu(sample_seg_rel_reads_df[
			sample_seg_rel_reads_df['segment'].isin(['PB2', 'PB1', 'PA', 'NS'])]['Rel_support'],
		sample_seg_rel_reads_df[
			~sample_seg_rel_reads_df['segment'].isin(['PB2', 'PB1', 'PA', 'NS'])]['Rel_support'])

	output.append('mann whitney u test ')
	output[-1] += 'comparing relative DVG read support on PB2, PB1, PA, NS v. HA, NP. NA, M'
	output[-1] += f': statistic={t.statistic}; pvalue={t.pvalue}' 
	
	# statistical test, unique DVG support in PB2, PB1, PA v. HA, NP, NA, M, NS
	t2 = mannwhitneyu(sample_seg_n_dvgs[
			sample_seg_n_dvgs['segment'].isin(['PB2', 'PB1', 'PA'])][0],
		sample_seg_n_dvgs[
			~sample_seg_n_dvgs['segment'].isin(['PB2', 'PB1', 'PA'])][0])

	output.append('Mann-Whtiney U test ')
	output[-1] += 'comparing number of unique DVGs per sample on PB2, PB1, PA v. HA, NP, NA, M, NS'
	output[-1] += f': statistic={t2.statistic}; pvalue={t2.pvalue}' 
	

	with open('figures/final/figure_1.txt', 'w') as fp:
		for line in output:
			fp.write(line + '\n')

	
	# finally, plot



	plot_style()
	fig = plt.figure(figsize=(6.4*3, 4.8*2.5), constrained_layout=True)
	# heatmap
	cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", ["white", "#333333"])
	spec = fig.add_gridspec(ncols=10, nrows=10)
	ax0 = fig.add_subplot(spec[:, :3])
	im = ax0.imshow(dat_arr, aspect='auto', cmap=cmap)
	cax = make_axes_locatable(ax0).append_axes('top', size='5%', pad=0.75)
	cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
	cbar.ax.set_xlabel(r'log$_{10}$(rel. read support)', size=16)
	ax0.set_xticks(list(segment_order.values()))
	ax0.set_xticklabels(list(segment_order.keys()))
	ax0.set_xlabel('segment', fontsize=20)
	ax0.set_yticks([])
	ax0.set_ylabel('sample', fontsize=20)
	
	

	# unique DVGs per sample per segment
	ax2 = fig.add_subplot(spec[:5, 4:7])
	ax2 = jitter_boxplot(ax=ax2, 
		i=[segment_order[i] for i in sample_seg_rel_reads_dict.keys()],
		d=sample_seg_rel_reads_dict.values())
	ax2.set_xticks(list(segment_order.values()))
	ax2.set_xticklabels(segment_order.keys())
	ax2.set_xlabel('segment', size=mpl.rcParams['axes.titlesize'])
	ax2.set_ylabel('total rel. DVG reads per sample', size=mpl.rcParams['axes.titlesize'])
	ax2.grid(axis='y', color='#eaeaea')

	# total relative dvg read support per sample per segment
	# unique DVGs per sample per segment
	ax3 = fig.add_subplot(spec[:5, 7:10])
	ax3 = jitter_boxplot(ax=ax3, 
		i=[segment_order[i] for i in sample_seg_n_dvgs_dict.keys()],
		d=sample_seg_n_dvgs_dict.values())
	ax3.set_xticks(list(segment_order.values()))
	ax3.set_xticklabels(segment_order.keys())
	ax3.set_xlabel('segment', size=mpl.rcParams['axes.titlesize'])
	ax3.set_ylabel('unique DVG species per sample', size=mpl.rcParams['axes.titlesize'])
	ax3.grid(axis='y', color='#eaeaea')


	
	# representative junction
	# pb2
	# need to read in and parse gff3
	ax4 = fig.add_subplot(spec[5:, 4:6])
	for idx, row in dat[(dat['Specid'] == args.repSpecid) & (dat['segment'] == 'PB2')].iterrows():
		_ = ax4.plot([row['mapped_start'] - args.padSize, row['mapped_stop'] - args.padSize], 
			[1,0], color='#333333', lw=0.5)

	
	_ = [ax4.spines[i].set_visible(False) for i in ['left', 'right']]
	ax4.set_ylim(0,1)
	ax4.set_xlabel('position (nt)', size=mpl.rcParams['axes.titlesize'])
	ax4.set_xlim(0, segment_len['PB2'])
	ax4.set_title('PB2')
	ax4.set_yticks([])


	# pb1
	ax5 = fig.add_subplot(spec[5:, 6:8])
	for idx, row in dat[(dat['Specid'] == args.repSpecid) & (dat['segment'] == 'PB1')].iterrows():
		_ = ax5.plot([row['mapped_start'] - args.padSize, row['mapped_stop'] - args.padSize], 
			[1,0], color='#333333', lw=0.5)

	_ = [ax5.spines[i].set_visible(False) for i in ['left', 'right']]
	ax5.set_ylim(0,1)
	ax5.set_xlabel('position (nt)', size=mpl.rcParams['axes.titlesize'])
	ax5.set_xlim(0, segment_len['PB1'])
	ax5.set_title(f'{args.repSpecid}\nPB1')
	ax5.set_yticks([])
	# pa
	ax6 = fig.add_subplot(spec[5:, 8:10])
	for idx, row in dat[(dat['Specid'] == args.repSpecid) & (dat['segment'] == 'PA')].iterrows():
		_ = ax6.plot([row['mapped_start'] - args.padSize, row['mapped_stop'] - args.padSize], 
			[1,0], color='#333333', lw=0.5)

	_ = [ax6.spines[i].set_visible(False) for i in ['left', 'right']]
	ax6.set_ylim(0,1)
	ax6.set_xlabel('position (nt)', size=mpl.rcParams['axes.titlesize'])
	ax6.set_xlim(0, segment_len['PA'])
	ax6.set_title(f'PA')
	ax6.set_yticks([])

	x_locs = [-0.07, -0.15, -0.175, -0.125, -0.125, -0.125]
	y_locs = [1.0, 1.0, 1.0, 1.05, 1.05, 1.05]
	for ax_idx, ax in enumerate([cbar.ax, ax2, ax3, ax4, ax5, ax6]):
			ax.text(x_locs[ax_idx], y_locs[ax_idx], 
				string.ascii_uppercase[ax_idx], color='#333333', 
				transform=ax.transAxes, size=18, fontweight='bold')

	fig.savefig('figures/final/figure_1.pdf')
	plt.close()


if __name__ == "__main__":
    run()








