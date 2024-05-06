import matplotlib.pyplot
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import string
try: 
	from utils import plot_style
except:
	from scripts.utils import plot_style
	

def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat', default=None)
	args = parser.parse_args()
	#args.allRuns =  "data/all_runs.tsv"
	#args.dat = 'output/parsed_dvgs_0.005_True_True_500_PB2_PB1_PA.tsv'
	#args.mapDir = "data/*_map.tsv"
	#args.fastaDir = "*_output/consensus/*.fasta"
	#args.nullPremature = "output/null_premature_stop.tsv"
	#args.padSize = 210
	dat = pd.read_csv(args.dat, sep='\t')	

	rng = np.random.default_rng(seed=111)
	# plot modulus
	dat = dat.assign(dvg_len_mod = lambda k: (k.Stop - k.Start - 1)%3)
	q_dat = {}
	# resample and calculate proportion
	for seg, seg_dat in dat.groupby('segment'):
		# resample mod
		samples = rng.choice(seg_dat.dvg_len_mod, 
				size=(seg_dat.dvg_len_mod.shape[0], 10000),
				replace=True, 
				p=seg_dat.Rel_support/seg_dat.Rel_support.sum())
		counted = np.apply_along_axis(lambda x: 
			np.bincount(x, minlength=3), axis=0, 
				arr=samples)
		proportioned = counted / counted.sum(axis=0)
		q = np.quantile(proportioned, [0, 0.025, 0.5, 0.975,1], axis=1)
		q_dat[seg] = q

	max_val = dat.Rel_support.max()*1.025
	min_val = -dat.Rel_support.max()*0.025

	output = []
	for seg in ['PB2', 'PB1', 'PA']:
		t = mannwhitneyu(
			dat.query('segment == @seg & dvg_len_mod == 0').Rel_support,
			dat.query('segment == @seg & dvg_len_mod > 0').Rel_support)
		output.append('mann whitney u test ')
		output[-1] += f'comparing relative DVG read support of {seg} DVGs with '
		output[-1] += f'mod == 0 to those with mod > 0 '
		output[-1] += f': statistic={t.statistic}; pvalue={t.pvalue}' 
	output.append(str(pd.concat([pd.DataFrame(value, columns=['mod = 0', 'mod = 1', 'mod = 2'], 
		index=['q0', 'q2.5', 'q50', 'q97.5', 'q1']).\
		assign(segment = key) for key, value in q_dat.items()])))
	with open('figures/final/figure_s10.txt', 'w') as fp:
		for line in output:
			fp.write(line + '\n')

		
	plot_style()
	segments = ['PB2', 'PB1', 'PA']
	fig, axs = plt.subplots(2,3,figsize=(6.4*3, 4.8*2), constrained_layout=True)
	for idx, seg in enumerate(segments):
		for col_idx, col in enumerate(q_dat[seg].T):
			axs[0,idx].plot([col_idx, col_idx], [col[1], col[3]], color='#333333', zorder=3)
		
		axs[0,idx].scatter([0,1,2], q_dat[seg][2,:], zorder=4, 
			facecolor='#eaeaea', edgecolor='#333333')
		axs[0,idx].set_xticks([0,1,2])
		if idx == 0:
			axs[0,idx].set_ylabel('Proportion of relative DVG reads', size=16)
		else:
			axs[0,idx].set_ylabel(' ', size=16)
		axs[0,idx].set_xlim(-0.5, 2.5)
		axs[0,idx].grid(axis='y', color='#eaeaea', zorder=0)
		axs[0,idx].set_title(seg)
		axs[0,idx].set_ylim(-0.025, 0.5)
		axs[1,idx].scatter(
			dat.query('segment == @seg').dvg_len_mod + rng.normal(0,0.1,size=dat.query('segment == @seg').shape[0]),
			dat.query('segment == @seg').Rel_support,
			facecolor='none', edgecolor='#333333', alpha=0.25, zorder=3)
		axs[1,idx].grid(axis='y', color='#eaeaea', zorder=0)
		#axs[1,idx].set_yscale("log")
		axs[1,idx].set_xticks([0,1,2])
		axs[1,idx].set_xlim(-0.5, 2.5)
		axs[1,idx].set_ylim(min_val, max_val)
		axs[1,idx].set_xlabel('deleted nucleotides % 3',size=16)
		if idx == 0:
			axs[1,idx].set_ylabel(r'relative read support', size=16)
		else:
			axs[1,idx].set_ylabel(r'', size=16)

	x_pos = [-0.13, -0.12, -0.12, -0.1, -0.08, -0.08]
	for ax_idx, ax in enumerate(axs.flatten()):
		ax.text(x_pos[ax_idx], 1.0, 
				string.ascii_uppercase[ax_idx], color='#333333', 
				transform=ax.transAxes, size=16, fontweight='bold')

	fig.suptitle('Figure S10')
	fig.savefig("figures/final/figure_s10.pdf")
	plt.close()


if __name__ == "__main__":
    run()
