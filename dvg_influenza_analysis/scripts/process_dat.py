import argparse
import pandas as pd
import glob
from collections import defaultdict
try: 
	from utils import format_map_dict
except:
	from scripts.utils import format_map_dict
	

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
	parser.add_argument('--titerTable', default=None)
	parser.add_argument('--dataDir', default=None)
	parser.add_argument('--dpDir', default=None)
	parser.add_argument('--mapDir', default=None)
	parser.add_argument('--segmentAcc', default=None)
	args = parser.parse_args()
	#args.runTable = 'data/SraRunTable.txt'
	#args.dataDir = 'hk_2014_2015_output/Virema/*.par'
	#args.dpDir = 'hk_2014_2015_output/depth/*_dp.tsv'
	#args.titerTable = 'data/elife-35962-fig1-data1-v3.csv'
	#args.mapDir = "data/*_map.tsv"
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
	id_dir_dict = {i.split('/')[-1].split('_')[0]: i for i in dir_files}
	dat = {i: 
		pd.read_csv(path, sep='\t') for i,path in id_dir_dict.items()}
	# add relative support coumn
	wt_dp_dict = {i.split('/')[-1].split('_')[0]: 
		{j[0]:j[1] for idx,j in pd.read_csv(i, sep='\t', header=None).\
					rename(columns={0: 'seg', 1: 'pos', 2: 'dp'}).\
					assign(mid_pos = lambda g: 
						g.groupby('seg').pos.transform(lambda k: int(k.max()/2))).\
					query('pos == mid_pos')[['seg', 'mid_pos']].iterrows()} 
			for i in glob.glob(args.dpDir)}
	for key, value in dat.items():
		value['Rel_support'] = value['Total_support']/value['Segment'].map(wt_dp_dict[key])
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
	# add sample titer information
	if args.titerTable:
		titers = pd.read_csv(args.titerTable, index_col=0).\
			rename(columns={
				'genomes.per.ul': 'genomes_per_ul',
				'SPECID': 'Specid'})\
			[['Specid', 'genomes_per_ul']]
		replicate_dat = replicate_dat.merge(titers, on='Specid', how='left')
	replicate_dat.to_csv('output/parsed_dvgs.tsv', sep='\t', index=None)



if __name__ == "__main__":
    run()




