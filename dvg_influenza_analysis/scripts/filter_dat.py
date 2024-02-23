import pandas as pd
import numpy as np
import argparse
try: 
	from utils import filter_dvg_dat
except:
	from scripts.utils import filter_dvg_dat


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat', default=None)
	parser.add_argument('--segments', default=None, nargs='+')
	parser.add_argument('--minDel', default=0, type=int)
	parser.add_argument('--minRelReads', default=0, type=float)
	parser.add_argument('--onlyShared', dest='only_shared', 
		default=False, action='store_true')
	parser.add_argument('--filterPlasmid', dest='filter_plasmid', 
		default=False, action='store_true')
	args = parser.parse_args()
	# read in and process data
	dat = pd.read_csv(args.dat, sep='\t', low_memory=False)
	# filter data
	filtered_dat = \
		filter_dvg_dat(dat, min_rel_reads=args.minRelReads, 
			only_shared=args.only_shared, min_del=args.minDel, filter_plasmid=args.filter_plasmid)
	# get only the segments we want
	if args.segments:
		filtered_dat = filtered_dat[filtered_dat['segment'].isin(args.segments)]

	if args.segments:
		filtered_dat.to_csv('.'.join(args.dat.split('.')[:-1])+ \
			f'_{args.minRelReads}_{args.only_shared}_{args.filter_plasmid}_{args.minDel}_{"_".join(args.segments)}.tsv',
			sep='\t', index=None)
	else:
		filtered_dat.to_csv('.'.join(args.dat.split('.')[:-1])+ \
			f'_{args.minRelReads}_{args.only_shared}_{args.filter_plasmid}_{args.minDel}_all.tsv',
			sep='\t', index=None)



if __name__ == "__main__":
    run()




