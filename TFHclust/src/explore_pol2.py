#%%
import pandas as pd
import numpy as np
import bioframe as bf
import matplotlib.pyplot as plt
import os
from importlib import reload
import detect_cluster
reload(detect_cluster)
#%%
ctrl_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/ctrl/ENCFF062SCH.bam"
obs_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/ENCFF343QJV.bam"
#%%
chromo = 'chr22'
bg_tbl = detect_cluster.get_chrom_read_tbl_from_bam(ctrl_bam,chromo,"/home/vipink/Documents/TFHclust/TFHclust/data/tmp").assign(start = lambda df_: df_.start.astype(int))
read_tbl = detect_cluster.get_chrom_read_tbl_from_bam(obs_bam,chromo,"/home/vipink/Documents/TFHclust/TFHclust/data/tmp").assign(start = lambda df_: df_.start.astype(int))
#%%
pol2_hdbscan_res = detect_cluster.perform_HDBScan_clustering(read_tbl,3,10)
full_graph = pol2_hdbscan_res.condensed_tree_.to_networkx()
tmp_cl_tbl = detect_cluster.collect_hdb_cluster_read(full_graph)
tmp_cl_data_tbl = detect_cluster.build_cluster_data_df(tmp_cl_tbl,read_tbl)
# %%
bg_coord_df = bg_tbl.loc[:,['chrom','start']].assign(start = lambda df_: df_.start.astype(int)).assign(end = lambda df_: df_.start + 1)
hdb_cluster_bg_rc_df = bf.count_overlaps(tmp_cl_data_tbl.assign(chrom = chromo).reset_index().loc[:,['chrom','start','end','HDB_cluster']],bg_coord_df)
tmp_cl_data_tbl = tmp_cl_data_tbl.merge(hdb_cluster_bg_rc_df.rename(columns={'count':'bg_count'}))
# %%
tmp_cl_data_tbl.plot.scatter(x='bg_count',y='rc',s=0.1,logx=True,logy=True,xlabel = "control read count", ylabel = 'Obs. read count')
# %%
data_tbl = (tmp_cl_data_tbl
            .loc[:,['HDB_cluster','rc','bg_count']]
            .assign(lrc = lambda df_: np.log10(df_.rc), lbg = lambda df_: np.log10(df_.bg_count))
            .assign(lbg2 = lambda df_: np.where(df_.bg_count.lt(1),0,df_.lbg))
            )
mod_res_tbl = detect_cluster.build_hic_zscore(data_tbl,['lbg2'],[10],[3],'lrc')
data_tbl = data_tbl.merge(mod_res_tbl)
# %%
data_tbl.sort_values('zscore').plot.scatter(x='bg_count',y='rc',c='zscore',s=0.1,cmap='coolwarm',logx=True,logy=True)
# %%
