import numpy as np
import pandas as pd


def jitter_boxplot(ax, d, i, j=0.5, vert=True, zorder=4):
    # for each item in D get quantiles
    q = [np.quantile(i, q=[0.25, 0.5, 0.75]) for i in d]
    # for each item get outlier thresholds
    o = [(i[0] - 1.5*(i[2] - i[0]), i[2] + 1.5*(i[2] - i[0])) for i in q]
    # get just outliers
    outliers = [i[(i < o[idx][0]) | (i > o[idx][1])] for 
        idx, i in enumerate(d)]
    # get outlier x positions
    outlier_x = np.hstack([np.repeat(i, outliers[idx].shape[0]) for idx, i in enumerate(i)])
    # now add jitter
    rng = np.random.default_rng(seed=111)
    outlier_x = outlier_x + rng.uniform(low = -j/2, high=j/2, size=outlier_x.shape[0])
    # now flatten outliers
    outliers = np.hstack(outliers)
    # first boxplot without fliers
    _ = ax.boxplot(d, positions=i,
                sym = '',
                medianprops=dict(color='steelblue'), vert=vert, zorder=zorder)
    # now add the fliers
    if vert:
        _ = ax.scatter(outlier_x, outliers, edgecolor='#333333', facecolor="None", marker='o', zorder=zorder)
    elif vert == False:
         _ = ax.scatter(outliers, outlier_x, edgecolor='#333333', facecolor="None", marker='o', zorder=zorder)
    return(ax)

def map_arr(a, d):
    u,inv = np.unique(a,return_inverse = True)
    a_map = np.array([d[x] for x in u])[inv].reshape(a.shape)
    return(a_map)   

    
def format_map_dict(map_files):
    map_dfs = [pd.read_csv(i, sep='\t', index_col=0) for i in map_files]
    map_dict = {}
    for seg_map in map_dfs:
        for name, name_map in seg_map.items():
            map_dict[name] = {int(i):idx for idx, i in name_map.items() if not pd.isna(i)} 
    return(map_dict)



def filter_dvg_dat(dat, only_shared=False, min_reads=0, min_del=0, filter_plasmid=False):
    # replace nan support values with 0
    # if we want only shared DVGs
    # note: this *only* applies to samples with technical rpelicates
    # otherwise we just take the identified reads
    replicated_dat = dat[(~dat['Run_y'].isnull()) & 
        (~dat['Run_x'].isnull())].copy()
    unreplicated_dat = dat[(dat['Run_y'].isnull()) | 
        (dat['Run_x'].isnull())].copy()
    if only_shared:
        replicated_dat = replicated_dat[(~replicated_dat['Total_support_x'].isnull())
            & (~replicated_dat['Total_support_y'].isnull())]
    # get average read support for eplicated data
    replicated_dat['Total_support'] = \
        replicated_dat[['Total_support_x', 'Total_support_y']].fillna(0).mean(axis=1)
    replicated_dat['Rel_support'] = \
        replicated_dat[['Rel_support_x', 'Rel_support_y']].fillna(0).mean(axis=1)
    # for unreplicated dat, we just take raw support values
    unreplicated_dat['Total_support'] = \
        unreplicated_dat[['Total_support_x', 'Total_support_y']].fillna(0).max(axis=1)
    unreplicated_dat['Rel_support'] = \
        unreplicated_dat[['Rel_support_x', 'Rel_support_y']].fillna(0).max(axis=1)
    # now combine and filter
    filtered_dat = pd.concat([replicated_dat, unreplicated_dat])
    filtered_dat = \
        filtered_dat[(filtered_dat['Total_support'] >= min_reads)
        & (filtered_dat['Stop'] - filtered_dat['Start'] >= min_del)]
    if filter_plasmid == True: 
        filtered_dat = filtered_dat[(filtered_dat['plasmid_x'] != True) & 
            (filtered_dat['plasmid_y'] != True)]
    return(filtered_dat)


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def format_seqs_arr(s, n_seqs):
    n_seqs = int(n_seqs)
    size = int(len(s)/n_seqs)
    seqs_arr = \
        np.frombuffer(s.lower().encode(), dtype=np.int8)
    seqs_arr = np.copy(seqs_arr)
    #[(97, 'A'), (114, 'R'), (119, 'W'), (109, 'M'), (100, 'D'), (104, 'H'), (118, 'V'), 
    #(110, 'N'), (99, 'C'), (121, 'Y'), (115, 'S'), (109, 'M'), (98, 'B'), (104, 'H'), 
    #(118, 'V'), (110, 'N'), (117, 'U'), (121, 'Y'), (119, 'W'), (107, 'K'), (98, 'B'), 
    #(100, 'D'), (104, 'H'), (110, 'N'), (103, 'G'), (114, 'R'), (115, 'S'), (107, 'K'), 
    #(98, 'B'), (100, 'D'), (118, 'V'), (110, 'N'), (116, 'T')]
    seqs_arr = \
        seqs_arr.reshape((n_seqs, int(seqs_arr.shape[0]/n_seqs)))
    return(seqs_arr)





def import_fasta(fasta_path):
    s_names = []
    all_s = ''
    fh = open(fasta_path, 'rt')
    with fh as fasta:
        for h,s in read_fasta(fasta):
            s_names.append(h)
            all_s += s
    fh.close()
    s_arr = format_seqs_arr(all_s, len(s_names))
    return(s_arr, np.array(s_names))




def plot_style(grey='#333333'):
    import matplotlib as mpl
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = 'Arial'
    mpl.rcParams['font.weight'] = 'light'
    mpl.rcParams['text.color'] = grey
    mpl.rcParams['axes.labelcolor'] = grey
    mpl.rcParams['xtick.color'] = '#707070'
    mpl.rcParams['ytick.color'] = '#707070'
    # Font sizes
    mpl.rcParams['figure.titlesize'] = 16
    mpl.rcParams['axes.titlesize'] = 16
    mpl.rcParams['xtick.labelsize'] = 16
    mpl.rcParams['ytick.labelsize'] = 16
    # Border colors
    mpl.rcParams['axes.edgecolor'] = grey
    mpl.rcParams['grid.color'] = '#eaeaea'
    #mpl.rcParams['grid.zorder'] = 0
    # Legend
    mpl.rcParams['legend.fontsize'] = 16
    mpl.rcParams['legend.frameon'] = False


