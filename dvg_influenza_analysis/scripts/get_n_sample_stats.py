import argparse
import pandas as pd
import glob

def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat', default=None)
	parser.add_argument('--runTable', default=None)
	parser.add_argument('--dateDat', default=None)
	parser.add_argument('--dataDir', default=None)
	parser.add_argument('--libs', default=None, nargs='+')
	args = parser.parse_args()
	#args.runTable = 'data/SraRunTable.txt'
	#args.dataDir = "*_output/Virema/*.par"
	#args.libs = ['HK', 'vic', 'perth']
	#args.dateDat =  'data/elife-35962-fig1-data1-v3.csv'
	#args.libs = ['HK2', 'HK7', 'perth', 'perth_2', 'HK1', 'HK8', 'HK6', 'vic_2', 'vic']

	date_dat = pd.read_csv(args.dateDat, index_col=0).rename(columns={'SPECID': 'Specid'})
	date_dict = {i[0]: i[1] for i in date_dat.values}
	#str.split('_', expand=True).iloc[:,1].apply(
	#	lambda k: ''.join(i for i in k if not i.isdigit()))
	run_table = pd.read_csv(args.runTable)
	specid_dict = {i['Run']: i['SPECID'] for idx, i in run_table.iterrows()}
	run_table['lib'] = run_table['Library Name'].str.replace('_A', '').\
		str.replace('_B', '').\
		str.replace('_pool', '').\
		str.replace('_mp', '').\
		str.split('_', expand=True).iloc[:,1:].\
		apply(lambda k: k[1] if k[2] is None else '_'.join(k), axis=1)
		#apply(lambda k: ''.join([i for i in k.split('_')[0] if not i.isdigit()]))
	run_table['type'] = run_table['Isolate'].apply(lambda k: 'clinical' if not 'plasma' in k  else 'plasmid')
	lib_runs = run_table[run_table['lib'].isin(args.libs)].drop_duplicates()
	lib_runs['days_post_onset'] = lib_runs['SPECID'].map(date_dict)
	lib_runs.sort_values(by='type').to_csv('data/all_runs.tsv', sep='\t', index=None)
	lib_runs_counts	= lib_runs[['Run', 'lib', 'type']].\
		drop_duplicates().groupby(['lib', 'type']).size().\
		reset_index()
	lib_runs_counts.to_csv('data/lib_runs_counts.tsv', sep='\t', index=None)
	lib_specids_counts = lib_runs[['SPECID', 'lib', 'type']].copy()
	lib_specids_counts['lib'] = lib_specids_counts['lib'].apply(
			lambda k: ''.join([i for i in k.split('_')[0] if not i.isdigit()]))
	lib_specids_counts = lib_specids_counts.drop_duplicates().groupby(['lib', 'type']).size().\
		reset_index()
	lib_specids_counts.to_csv('data/lib_specids_counts.tsv', sep='\t', index=None)
	# count number of samples
	enrollid_runs_counts = \
		lib_runs[['Isolate', 'lib', 'type']].copy()
	enrollid_runs_counts = enrollid_runs_counts[enrollid_runs_counts['type'] == 'clinical']
	enrollid_runs_counts['lib'] = enrollid_runs_counts['lib'].apply(
		lambda k: ''.join([i for i in k.split('_')[0] if not i.isdigit()]))
	enrollid_runs_counts = enrollid_runs_counts.drop_duplicates().groupby(
			['lib', 'type']).size().reset_index()
	enrollid_runs_counts.to_csv('data/enrollids_runs_counts.tsv', sep='\t', index=None)
	# count number of longitudinal samples
	n_longitudinal =  lib_runs[lib_runs['type'] == 'clinical'][['Isolate', 'SPECID', 'days_post_onset']].\
	groupby('Isolate').apply(lambda k: 
		k['SPECID'].unique().shape[0] > 1).sum()
	output = []
	output.append(f'there are {n_longitudinal} individuals with longitudinal samples')
	with open('data/stats.tsv', 'w') as fp:
		for line in output:
			fp.write(line+'\n')
	# from data directory
	plasmid_dict = {i['lib']: i['Run'] for
		k, i in run_table[run_table['Isolate'].\
			str.contains('plasma')].iterrows()}
	dir_files = glob.glob(args.dataDir)
	# create mapping of id to actually directory
	id_dir_runs = [i.split('/')[-1].split('both')[0] for i in dir_files]
	# confirm that the runs match the run table
	# ensures everythign ran through the pipeline
	if not set(id_dir_runs) == set(lib_runs['Run'].values):
		raise Exception('files and run table do not match')


if __name__ == "__main__":
    run()