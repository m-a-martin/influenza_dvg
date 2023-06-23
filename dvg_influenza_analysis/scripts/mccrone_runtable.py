import pandas as pd
import os
import datetime


### From:  https://doi.org/10.7554/eLife.35962
# Briefly, reads were aligned to the reference sequence 
# H3N2 2010–2011 and 2011–2012: GenBank CY121496-503
# mm: these are labelled as "perth"
# H3N2 2012–2013: GenBank KJ942680-8
# mm: these are labelled as "vic"
# H3N2 2014–2015: Genbank CY207731-8
# mm: based on dates, these must be labelled as "hk"
# H1N1 GenBank: CY121680-8
# mm: these are labelled as "cali09"
lib_season_dict = {'perth': ('2010-07-01', '2012-06-30'),
       'vic': ('2012-07-01', '2013-06-30'),
       'hk': ('2014-07-01', '2015-06-30')}


dat = pd.read_csv('SraRunTable.txt')
dat['Collection_date'] = pd.to_datetime(dat['Collection_date'])
dat['season'] = dat['Collection_date'].apply(lambda k: f'{k.year-1}_{k.year}' if k.month < 7 else f'{k.year}_{k.year+1}')
dat['subtype'] = dat['Library Name'].apply(lambda k: 'H3N2' if 'cali' not in k else 'H1N1')
dat['lib'] = dat['Library Name'].str.split('_', expand=True).iloc[:,1:].apply(lambda k: '_'.join([i for i in k if i]), axis=1).str.replace('A_', '').\
       str.replace('B_', '').str.replace('_2', '').\
       str.replace('pool_', '').str.replace('mp_', '').apply(lambda k: 'hk' if 'HK' in k else k)


for i, i_dat in dat.groupby('lib'):
       if i != 'cali09':
              min_date = i_dat['Collection_date'].min()
              max_date = i_dat['Collection_date'].max()
              assert min_date > datetime.datetime.strptime(lib_season_dict[i][0], '%Y-%m-%d')
              assert max_date < datetime.datetime.strptime(lib_season_dict[i][1], '%Y-%m-%d')
              d = f'{i}_{min_date.year}_{max_date.year}'
       else:
           d = i
       if not os.path.exists(f'data/{d}'):
           os.mkdir(f'data/{d}')
       with open(f'data/{d}/download_{d}.sh', 'w') as fp:
              fp.write("#!/bin/sh\n")
              for k in i_dat['Run'].values:
                     fp.write(f'fasterq-dump -S -O data/{d} ')
                     fp.write(k+'\n')




