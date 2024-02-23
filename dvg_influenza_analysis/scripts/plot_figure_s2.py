import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf
import matplotlib as mpl
import string
try: 
	from utils import plot_style, jitter_boxplot
except:
	from scripts.utils import plot_style, jitter_boxplot


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat', default=None)
	parser.add_argument('--allRuns', default=None)
	args = parser.parse_args()
	#args.dat = 'output/parsed_dvgs_0.005_True_True_0_all.tsv'
	#args.allRuns = "data/all_runs.tsv"

	segment_order = {'PB2': 0, 'PB1': 1, 'PA': 2, 'HA': 3, 
		'NP': 4, 'NA': 5, 'M': 6, 'NS': 7}

	# get index to add 0s	
	all_runs = pd.read_csv(args.allRuns, sep='\t')
	all_runs = all_runs[(all_runs['type'] == 'clinical')]
	all_seg_runs_index = pd.MultiIndex.from_product([all_runs.SPECID.drop_duplicates(),
			pd.Series(segment_order.keys(), name='segment')])
	
	dat = pd.read_csv(args.dat, sep='\t')
	dat['segment'] = dat['segment'].fillna('NA')

	# get genomes per uL dict
	# implicitly assumes all samples are in our DVG dat which we know to be true
	genomes_dict = {i.Specid: i.genomes_per_ul for 
		idx,i in dat[['Specid', 'genomes_per_ul']].drop_duplicates().iterrows()}
	
	# group by specid and segment and get total DVG support
	# then add 0s through multi index
	seg_rel_support = dat.groupby(['Specid', 'segment']).agg({'Rel_support': sum}).\
		reindex(all_seg_runs_index).fillna(0.0).reset_index().\
		assign(
			genomes_per_ul = lambda k: k.SPECID.map(genomes_dict),
			log10_genomes_per_ul = lambda k: np.log10(k.genomes_per_ul))

	# check that all specids have genomes per ul
	if seg_rel_support.genomes_per_ul.min() <= 0:
		raise Exception('missing genomes per ul data')

	lin_regress = seg_rel_support.groupby('segment').\
		apply(lambda k: smf.ols(formula='Rel_support ~ log10_genomes_per_ul', 
			data=k).fit()).reset_index().\
		assign(order = lambda k: k.segment.map(segment_order)).sort_values(by='order')

	output = []
	for idx,i in lin_regress.iterrows():
		output.append(i.segment)
		output.append(str(i[0].summary()))

	with open('figures/final/figure_s2.txt', 'w') as fp:
		for line in output:
			fp.write(line+'\n')

	plot_style()
	fig, axs = plt.subplots(4,2, figsize=(6.4*1.5, 4.8*3), constrained_layout=True)
	for seg, idx in segment_order.items():
		ax = axs.flatten()[idx]
		seg_dat = seg_rel_support.query('segment == @seg')
		ax.scatter(seg_dat.genomes_per_ul, seg_dat.Rel_support, facecolor='none', edgecolor='#333333', alpha=0.5)
		ax.plot(np.linspace(seg_dat.genomes_per_ul.min(), seg_dat.genomes_per_ul.max()),
			lin_regress.query('segment == @seg')[0].values[0].predict(
				pd.DataFrame(np.log10(np.linspace(seg_dat.genomes_per_ul.min(), seg_dat.genomes_per_ul.max())), 
					columns=['log10_genomes_per_ul'])), color='steelblue')
		ax.set_xscale('log')
		ax.set_xlabel(r'genomes/$\mu$L', size=14)
		ax.set_ylabel('relative read support', size=14)
		ax.tick_params(axis='both', labelsize=12)
		ax.set_title('\n' + seg, size=18)
		ax.text(-0.125, 1.0, 
			string.ascii_uppercase[idx], color='#333333', 
			transform=ax.transAxes, size=18, fontweight='bold')

	fig.suptitle('Figure S2')
	fig.savefig('figures/final/figure_s2.pdf')
	plt.close()





if __name__ == "__main__":
    run()

