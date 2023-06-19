import pandas as pd
import numpy as np
import argparse
try:
	from utils import import_fasta
except:
	from scripts.utils import import_fasta


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--seqs', default=None)
	parser.add_argument('--seqNameSep', default=' ')
	parser.add_argument('--name', default='')
	args = parser.parse_args()
	#args.seqs = 'data/pb2_ref_aln.fasta'
	s_arr, s_names = import_fasta(args.seqs)
	# format names
	s_names = np.array([i.split(args.seqNameSep)[0] for i in s_names])
	out = np.full(s_arr.shape, np.nan)
	out[np.where(s_arr != 45)] = \
		np.hstack([np.arange(i) for i in (s_arr != 45).sum(axis=1)])
	# adding 1 because DVG coordinates are 1-indexed 
	pd.DataFrame(out+1, index=s_names, columns=np.arange(out.shape[1])+1).T.rename_axis(args.name).to_csv('.'.join(args.seqs.split('.')[:-1])+'_map.tsv', sep='\t')




if __name__ == "__main__":
    run()



