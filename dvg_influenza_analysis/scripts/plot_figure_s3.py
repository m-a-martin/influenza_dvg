import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
try: 
	from utils import plot_style
except:
	from scripts.utils import plot_style

def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat', default=None)
	args = parser.parse_args()
	segment_order = {'PB2': 0, 'PB1': 1, 'PA': 2, 'HA': 3, 
		'NP': 4, 'NA': 5, 'M': 6, 'NS': 7}
	# read in and process data
	dat = pd.read_csv(args.dat, sep='\t')
	dat['segment'] = dat['segment'].fillna("NA")
	# count number of deleted bases
	dat['n_deleted'] = dat['mapped_stop'] - dat['mapped_start']

	bw = 50
	n_segs = dat['segment'].unique().shape[0]
	plot_style()
	fig, axs = plt.subplots(1, 3,
	figsize=(6.4*2, 4.8*0.85), constrained_layout=True)
	# now let's sum up across everything and do a cumulative plot
	for segment, segment_dat in dat.groupby('segment'):
		bins = np.arange(0, segment_dat['n_deleted'].max()+bw, bw)
		_ = axs.flatten()[segment_order[segment]].hist(segment_dat['n_deleted'], 
			bins=bins, facecolor='#eaeaea', edgecolor='#333333', zorder=2)
		axs.flatten()[segment_order[segment]].axvline(500, ls='--', color='#4d4d4d')
		axs.flatten()[segment_order[segment]].grid(color='#eaeaea', zorder=0)
		axs.flatten()[segment_order[segment]].set_title(segment)

	fig.suptitle('S3 Fig')
	fig.supxlabel('deleted nucleotides', size=16)
	fig.supylabel('DVG count', size=16)
	fig.savefig(f'figures/final/figure_s3.pdf')
	plt.close('all')


if __name__ == "__main__":
    run()





