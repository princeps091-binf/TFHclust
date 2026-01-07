#%%
import pandas as pd
import numpy as np
import bioframe as bf
import matplotlib.pyplot as plt
import os
from importlib import reload
import detect_cluster
reload(detect_cluster)
import integrate_replicate
reload(integrate_replicate)
#%%
ctrl_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/ctrl/ENCFF062SCH.bam"
rep1_obs_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/ENCFF343QJV.bam"
rep2_obs_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/ENCFF386OKF.bam"

#%%
chromo = 'chr22'
bg_tbl = detect_cluster.get_chrom_read_tbl_from_bam(ctrl_bam,chromo,"/home/vipink/Documents/TFHclust/TFHclust/data/tmp").assign(start = lambda df_: df_.start.astype(int))
rep1_read_tbl = detect_cluster.get_chrom_read_tbl_from_bam(rep1_obs_bam,chromo,"/home/vipink/Documents/TFHclust/TFHclust/data/tmp").assign(start = lambda df_: df_.start.astype(int))
rep2_read_tbl = detect_cluster.get_chrom_read_tbl_from_bam(rep2_obs_bam,chromo,"/home/vipink/Documents/TFHclust/TFHclust/data/tmp").assign(start = lambda df_: df_.start.astype(int))
# %%
# %%
rep1_clustering = detect_cluster.read_clustering(rep1_read_tbl)
rep1_clustering.HDBScan_clustering(3,10)
rep1_clustering.collect_hdb_cluster_read()
rep1_clustering.build_regression_tbl(bg_tbl)
rep1_clustering.compute_zscore(['lbg2'],[10],[3],'lrc')

rep2_clustering = detect_cluster.read_clustering(rep2_read_tbl)
rep2_clustering.HDBScan_clustering(3,10)
rep2_clustering.collect_hdb_cluster_read()
rep2_clustering.build_regression_tbl(bg_tbl)
rep2_clustering.compute_zscore(['lbg2'],[10],[3],'lrc')

# %%
rep1_clustering.summary_tbl.sort_values('zscore').plot.scatter(x='bg_count',y='rc',c='zscore',s=0.1,cmap='coolwarm',logx=True,logy=True)

# %%
rep2_clustering.summary_tbl.sort_values('zscore').plot.scatter(x='bg_count',y='rc',c='zscore',s=0.1,cmap='coolwarm',logx=True,logy=True)

# %%
merged_hic_clustering = integrate_replicate.merged_clustering(rep1_clustering,rep2_clustering)
merged_hic_clustering.match_rep_cluster()
merged_hic_clustering.evaluate_overlap_significance(10)
merged_hic_clustering.filter_significant_overlap_cluster(0.5,0.5)
# %%
merged_hic_clustering.inter_rep_match_tbl.sort_values('w')