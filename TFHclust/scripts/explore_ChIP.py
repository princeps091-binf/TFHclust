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

ctrl_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/ctrl/ENCFF062SCH.bam"
rep1_obs_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/rep3_ENCFF343QJV.bam"
rep2_obs_bam = "/home/vipink/Documents/TFHclust/TFHclust/data/POLR2A/rep1_ENCFF386OKF.bam"

# %%
chromo = 'chr22'

bg_tbl = src.detect_cluster.get_chrom_read_tbl_from_bam(ctrl_bam,chromo)
rep1_read_tbl = src.detect_cluster.get_chrom_read_tbl_from_bam(rep1_obs_bam,chromo)
rep2_read_tbl = src.detect_cluster.get_chrom_read_tbl_from_bam(rep2_obs_bam,chromo)


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
