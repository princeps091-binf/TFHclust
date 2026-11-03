# %%
import pandas as pd
import numpy as np
import bioframe as bf
import matplotlib.pyplot as plt
from matplotlib.cm import viridis
import os
from importlib import reload
import src.detect_cluster
reload(src.detect_cluster)
import src.integrate_replicate
reload(src.integrate_replicate)

import matplotlib.patches as pat  # Patches like pat.Polygon()
from matplotlib.collections import PolyCollection  # Collections of patches

# %%
ctrl_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/ctrl/ENCFF062SCH.bam"
rep1_obs_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/rep3_ENCFF343QJV.bam"
rep2_obs_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/rep1_ENCFF386OKF.bam"

# %%
chromo = 'chr19'
bg_tbl = src.detect_cluster.get_chrom_read_tbl_from_bam(ctrl_bam,chromo,"/home/vipink/Documents/TFHclust/TFHclust/data/tmp").assign(start = lambda df_: df_.start.astype(int))
rep1_read_tbl = src.detect_cluster.get_chrom_read_tbl_from_bam(rep1_obs_bam,chromo,"/home/vipink/Documents/TFHclust/TFHclust/data/tmp").assign(start = lambda df_: df_.start.astype(int))
rep2_read_tbl = src.detect_cluster.get_chrom_read_tbl_from_bam(rep2_obs_bam,chromo,"/home/vipink/Documents/TFHclust/TFHclust/data/tmp").assign(start = lambda df_: df_.start.astype(int))

# %%
rep1_clustering = src.detect_cluster.read_clustering(rep1_read_tbl)
rep1_clustering.HDBScan_clustering(3,10)
rep1_clustering.collect_hdb_cluster_read()
rep1_clustering.build_regression_tbl(bg_tbl)
rep1_clustering.compute_zscore(['lbg2'],[10],[3],'lrc')

rep2_clustering = src.detect_cluster.read_clustering(rep2_read_tbl)
rep2_clustering.HDBScan_clustering(3,10)
rep2_clustering.collect_hdb_cluster_read()
rep2_clustering.build_regression_tbl(bg_tbl)
rep2_clustering.compute_zscore(['lbg2'],[10],[3],'lrc')

# %%
rep1_clustering.summary_tbl.sort_values('zscore').plot.scatter(x='bg_count',y='rc',c='zscore',s=0.1,cmap='coolwarm',logx=True,logy=True)
plt.show()
# %%
rep2_clustering.summary_tbl.sort_values('zscore').plot.scatter(x='bg_count',y='rc',c='zscore',s=0.1,cmap='coolwarm',logx=True,logy=True)
plt.show()
# %%
merged_hic_clustering = src.integrate_replicate.merged_clustering(rep1_clustering,rep2_clustering)
merged_hic_clustering.match_rep_cluster()
merged_hic_clustering.evaluate_overlap_significance(10)
merged_hic_clustering.filter_significant_overlap_cluster(0.5,0.5)
# %%
merged_hic_clustering.inter_rep_match_tbl.plot.scatter(x='zscore_rep1',y='zscore_rep2')
plt.show()
# %%
# 1. Define your target interval
x_min, x_max = 1000000, 1250000

# 2. Filter data where the triangle falls within [x_min, x_max]
tmp_plt_df = merged_hic_clustering.inter_rep_match_tbl.sort_values('norm_lvl_rep1')
plt_df = tmp_plt_df[(tmp_plt_df.end >= x_min) & (tmp_plt_df.start <= x_max)]


triangles = []
for _,t in plt_df.iterrows():
    xmid = t.start + (t.w/2)  # Middle x-coord
    xleft = t.start
    xright = t.end
    y1 = 0  # y-coords
    y2 = (t.w/2)
    coordinates = [[xleft, y1], [xright, y1], [xmid, y2]]
    triangles.append(coordinates)  # Append to collection

z = (plt_df.norm_lvl_rep2 + plt_df.norm_lvl_rep1)/2
collec = PolyCollection(triangles, array=z, cmap=viridis)

fig, ax = plt.subplots(1, 1)
ax.add_collection(collec)  # Plot polygon collection
ax.autoscale_view()
plt.show()

# %%
chrom_peak_df = (pd.read_csv('/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/bed/rep1_ENCFF320MOY.bed.gz',sep='\t',usecols=[0,1,2,4,9],header=None)
 .rename(columns={0:'chrom',
                  1:'start',
                  2:'end',
                  4:'IDR',
                  9:'w'})
 .assign(start = lambda df: df.start.astype(int),
         end = lambda df: df.end.astype(int))).query('chrom == @chromo')
# %%
(bf.count_overlaps(chrom_peak_df,rep2_read_tbl.assign(end = lambda df: df.start + 1))
.rename(columns={'count':'obs'})
.merge(
bf.count_overlaps(chrom_peak_df,bg_tbl.assign(end = lambda df: df.start + 1))
.rename(columns={'count':'bg'})
)
.assign(FC = lambda df: df.obs/(df.bg + 1))
.FC.plot.kde()
)
plt.show()
# %%
(bf.count_overlaps(rep2_clustering.summary_tbl,chrom_peak_df)
 .assign(peak_in = lambda df: np.where(df.loc[:,'count'].lt(1),'out','in'))
 .groupby(['peak_in'])
 .pvalue.plot.kde(legend=True)
)
plt.show()
#%%
merged_hic_clustering.inter_rep_match_tbl.assign(lw=lambda df: np.log10(df.w)).lw.plot.kde()
#%%
bf.merge(merged_hic_clustering.inter_rep_match_tbl.loc[:,['chrom','start','end','pvalue_rep1','pvalue_rep2']]).sort_values('n_intervals')

#%%
#%%
target_domains_df = (merged_hic_clustering.inter_rep_match_tbl
                     .sort_values('norm_lvl_rep1')
                     .query('start > 18588000 and end <18591176'))
triangles = []
for _,t in target_domains_df.iterrows():
    xmid = t.start + (t.w/2)  # Middle x-coord
    xleft = t.start
    xright = t.end

    y1 = 0  # y-coords
    y2 = (t.w/2)

    coordinates = [[xleft, y1], [xright, y1], [xmid, y2]]

    print(coordinates) 
    triangles.append(coordinates)  # Append to collection

z = (target_domains_df.norm_lvl_rep2 + target_domains_df.norm_lvl_rep1)/2
collec = PolyCollection(triangles, array=z, cmap=viridis)

fig, ax = plt.subplots(1, 1)
ax.add_collection(collec)  # Plot polygon collection
ax.autoscale_view()
plt.show()

# %%
