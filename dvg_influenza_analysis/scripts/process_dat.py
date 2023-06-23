import argparse
import pandas as pd
import glob
from collections import defaultdict
from utils import format_map_dict
import numpy as np


def get_replicate_dat(run_table, dat):
	replicate_dat = []
	for specid, specdat in run_table.groupby('SPECID'):
		# check that the run table is properly formatted
		# either 1 or 2 runs per specid
		if specdat.shape[0] == 0:
			raise Exception('no data available')
		elif specdat.shape[0] > 2:
			raise Exception('three technical replicates unexpected')
		# check if *any* of the runs for this specid are in the dat dictionary
		if specdat['Run'].isin(dat.keys()).any():
			# if some are missing
			if not specdat['Run'].isin(dat.keys()).all():
				raise Exception('incorrect number of runs are in data directory')
			# if there are two runs for this sapmle, then we merge the data to find their union
			if specdat.shape[0] == 2:
				i_dat = dat[specdat['Run'].iloc[0]].merge(dat[specdat['Run'].iloc[1]], 
					on=['Segment', 'Start', 'Stop'],
					how='outer').fillna(0)
				i_dat['Specid'] = specid
				i_dat['Run_x'] = specdat['Run'].iloc[0]
				i_dat['Run_y'] = specdat['Run'].iloc[1]
			# else if there is only one run for this sample, just take the data for that run
			elif specdat.shape[0] == 1:
				i_dat = dat[specdat['Run'].iloc[0]]
				i_dat['Specid'] = specid
				i_dat['Total_support_x'] = i_dat['Total_support']
				i_dat['Total_support_y'] = None
				i_dat['Rel_support_x'] = i_dat['Rel_support']
				i_dat['Rel_support_y'] = None
				i_dat['Run_x'] = specdat['Run'].iloc[0]
				i_dat['Run_y'] = None
			else:
				raise Exception('unexpected number of replicates')
			replicate_dat.append(i_dat[['Specid', 'Run_x', 'Run_y', 'Segment', 'Start', 'Stop', 
				'Total_support_x', 'Total_support_y', 'Rel_support_x', 'Rel_support_y']])
	return(pd.concat(replicate_dat))


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--runTable', default=None)
	parser.add_argument('--dataDir', default=None)
	parser.add_argument('--dpDir', default=None)
	parser.add_argument('--mapDir', default=None)
	parser.add_argument('--segmentAcc', default=None)
	args = parser.parse_args()
	#args.runTable = 'data/SraRunTable.txt'
	#args.dataDir = 'hk_2014_2015_output/Virema/*.par'
	#args.dpDir = 'hk_2014_2015_output/depth/*_dp.tsv'
	run_table = pd.read_csv(args.runTable)
	run_table['lib'] = run_table['Library Name'].str.replace('_A', '').\
		str.replace('_B', '').\
		str.replace('_pool', '').\
		str.replace('_mp', '').\
		str.split('_', expand=True).iloc[:,1:].\
		apply(lambda k: k[1] if k[2] is None else '_'.join(k), axis=1)
	# create a dictionary linking each library to the run ID of its plasma control
	plasmid_dict = {i['lib']: i['Run'] for
		k, i in run_table[run_table['Isolate'].\
			str.contains('plasma')].iterrows()}
	lib_dict = {i['Run']: i['lib'] for k, i in run_table.iterrows()}
	# map each specid to enrollid
	enrollid_dict = {i['SPECID']: i['Isolate'] for k, i in run_table.iterrows()}
	# read in all files in data directory
	dir_files = glob.glob(args.dataDir)
	# create mapping of id to actually directory
	id_dir_dict = {i.split('/')[-1].split('both')[0]: i for i in dir_files}
	dat = {i: 
		pd.read_csv(path, sep='\t') for i,path in id_dir_dict.items()}
	# add relative support coumn
	wt_dp_dat = {i.split('/')[-1].split('both')[0]: 
				pd.read_csv(i, sep='\t', header=None,
					names=['Segment', 'length', 'reads_mapped', 'reads_unmapped'])
			for i in glob.glob(args.dpDir)}
	for key, value in dat.items():
		val_seg_dvg_reads = value.groupby(['Segment'])['Total_support'].sum()
		val_seg_all_reads = wt_dp_dat[key].merge(val_seg_dvg_reads, on=['Segment'])
		all_reads_dict = {i[0]: i[2]+i[4] for i in val_seg_all_reads.values}
		value['Rel_support'] = value['Total_support']/value['Segment'].map(all_reads_dict)
	# for each library, 
	# get set of identified DIPs if there is data 
	plasmid_dat = {key: defaultdict(lambda: False) for key in plasmid_dict.keys()}
	for key, value in plasmid_dict.items():
		if value in dat.keys():
			plasmid_dat[key].update({(i[0], i[1], i[2]): True for i in dat[value].values})
	# mapping of each accession to its library
	# now read in data, merge if technical replicate, don't if not
	# only do this for non plasmid controls
	replicate_dat = \
		get_replicate_dat(run_table[~run_table['Isolate'].str.contains('plasma')], dat)
	print(replicate_dat)
	# for each row, add a bool for whether that dvg is in X plasmid and Y plasmid
	# tuple ID for each dvg so we can look them up in the plasmid dat dictionaries 
	dvgs = list(zip(replicate_dat['Segment'], replicate_dat['Start'], replicate_dat['Stop']))
	replicate_dat['plasmid_x'] = \
		[plasmid_dat[i[1]][i[0]] for i in 
			zip(dvgs, replicate_dat['Run_x'].map(lib_dict))]
	replicate_dat['plasmid_y'] = \
		[plasmid_dat[i[1]][i[0]] if not pd.isna(i[1]) else None for i in 
			zip(dvgs, replicate_dat['Run_y'].map(lib_dict))]
	if args.segmentAcc:
	# add a column with the actual segment string
		segment_acc = pd.read_csv(args.segmentAcc, sep='\t', header=None).fillna("NA")
		segment_acc_dict = {i[0]: i[1] for i in segment_acc.values}
		replicate_dat['segment'] = replicate_dat['Segment'].map(segment_acc_dict)
	# add mapped start and stop locations
	if args.mapDir:
		map_dict = format_map_dict(glob.glob(args.mapDir))
		replicate_dat['mapped_start'] = \
			replicate_dat.apply(lambda k: map_dict[k['Segment']][k['Start']], axis=1)
		replicate_dat['mapped_stop'] = \
			replicate_dat.apply(lambda k: map_dict[k['Segment']][k['Stop']], axis=1)
	replicate_dat['enrollid'] = replicate_dat['Specid'].map(enrollid_dict)
	replicate_dat.to_csv('dvg_influenza_analysis/output/parsed_dvgs.tsv', sep='\t', index=None)



if __name__ == "__main__":
    run()




