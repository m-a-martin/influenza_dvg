import matplotlib.pyplot
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import glob
import string
try: 
	from utils import filter_dvg_dat
except:
	from scripts.utils import filter_dvg_dat




def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--rawDat', default=None)
	parser.add_argument('--filteredDat', default=None)
	parser.add_argument('--minRelReads', default=0, type=float)
	parser.add_argument('--padSize', default=None, type=int)
	parser.add_argument('--onlyShared', dest='only_shared', 
		default=False, action='store_true')
	parser.add_argument('--filterPlasmid', dest='filter_plasmid', 
		default=False, action='store_true')
	parser.add_argument('--minDel', default=0, type=int)
	args = parser.parse_args()
	#args.filteredDat = 'output/parsed_dvgs_0.005_True_True_0_all.tsv'
	#args.rawDat = 'output/parsed_dvgs.tsv'
	#args.minRelReads = 5E-3
	#args.padSize = 210
	#args.only_shared = True
	
	
	raw_dat = pd.read_csv(args.rawDat, sep='\t').\
		assign(mapped_start = lambda k: k.mapped_start - args.padSize, 
			mapped_stop = lambda k: k.mapped_stop - args.padSize)
	raw_dat['segment'] = raw_dat['segment'].fillna('NA')
	raw_dat = filter_dvg_dat(raw_dat, 
		only_shared=args.only_shared, 
		min_rel_reads=args.minRelReads, 
		filter_plasmid=args.filter_plasmid,
		min_del=args.minDel)
	filtered_dat = pd.read_csv(args.filteredDat, sep='\t').\
		assign(mapped_start = lambda k: k.mapped_start - args.padSize, 
			mapped_stop = lambda k: k.mapped_stop - args.padSize)
	filtered_dat['segment'] = filtered_dat['segment'].fillna('NA')
	
	# count number of observed DVGs
	raw_n = raw_dat.groupby(['segment', 'mapped_start', 'mapped_stop']).size().reset_index()
	filtered_n = filtered_dat.groupby(['segment', 'mapped_start', 'mapped_stop']).size().reset_index()

	n_combined = raw_n.merge(filtered_n, how='outer', on=['segment', 'mapped_stop', 'mapped_start'])

	of = 'output/' + \
		'_'.join([
			args.filteredDat.split('/')[-1].replace('.tsv',''), 
			args.rawDat.split('/')[-1].replace('.tsv',''),
			str(args.minRelReads),
			str(args.only_shared),
			str(args.filter_plasmid), 
			str(args.minDel)]) + '.tsv'

	n_combined.sort_values(by=['0_x', '0_y'], ascending=False).to_csv(of,
			sep='\t', index=None)




if __name__ == "__main__":
    run()

	