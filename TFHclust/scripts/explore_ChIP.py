# %%
# run from project root (ease of module access)

import pandas as pd
import numpy as np
import bioframe as bf
import matplotlib.pyplot as plt
import os
from importlib import reload
import src.detect_cluster
reload(src.detect_cluster)
import src.integrate_replicate
reload(src.integrate_replicate)

import matplotlib.patches as pat  # Patches like pat.Polygon()
from matplotlib.collections import PolyCollection  # Collections of patches

# %%
chromo = 'chr19'

# %%

# 1. Define your replicate paths or raw data containers
replicate_configs = {
    'rep1': {'reads': "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/rep3_ENCFF343QJV.bam", 'control': "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/ctrl/ENCFF062SCH.bam"},
    'rep2': {'reads': "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/rep1_ENCFF386OKF.bam", 'control': "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/ctrl/ENCFF062SCH.bam"}
    }

# 2. Instantiate and run the ReadClustering pipeline for each replicate
clustering_objs_dict = {}

for rep_name, config in replicate_configs.items():
    tmp_read_tbl = src.detect_cluster.get_chrom_read_tbl_from_bam(config['reads'],chromo)
    # Initialize your custom ReadClustering class (assuming standard naming from earlier context)
    tmp_ctrl_tbl = src.detect_cluster.get_chrom_read_tbl_from_bam(config['control'],chromo)
    
    clusterer = src.detect_cluster.read_clustering(tmp_read_tbl)
    
    # Run your upstream processing steps to populate summary_tbl and read tables
    clusterer.HDBScan_clustering(3,10)
    clusterer.collect_hdb_cluster_read()
    clusterer.build_regression_tbl(tmp_ctrl_tbl)
    clusterer.compute_zscore(['lbg2'],[10],[3],'lrc')
    
    # Store the fitted object in the dictionary
    clustering_objs_dict[rep_name] = clusterer

# %%
multi_rep_obj = src.integrate_replicate.DirectedMultiRepGraphClustering(clustering_objs_dict)

multi_rep_obj.build_graph(0.5)
multi_rep_obj.parallel_compute_cvm_edge_weights(max_workers=8)
multi_rep_obj.extract_consensus_components(jaccard_thresh=0.5,cvm_pval_thresh=0.5)

# %%

x_min, x_max = 48000000,50000000
tmp_plt_df = multi_rep_obj.consensus_loci_tbl.assign(w = lambda df: df.end - df.start).sort_values('w',ascending=False)
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

# z = (plt_df.norm_lvl_rep2 + plt_df.norm_lvl_rep1)/2
collec = PolyCollection(triangles, alpha=0.3, facecolors='skyblue')

fig, ax = plt.subplots(1, 1)
ax.add_collection(collec)  # Plot polygon collection
ax.autoscale_view()
plt.show()



