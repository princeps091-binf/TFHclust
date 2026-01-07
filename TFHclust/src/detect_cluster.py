import pysam
import pandas as pd
import hdbscan
import networkx as nx
import statsmodels.api as sm
from scipy.stats import t
import subprocess

def get_chrom_read_tbl_from_bam(bam_file,tmp_chrom,tmp_folder):
    tmp_file = f"{tmp_folder}/tmp.bam"
    tmp_index_file = f"{tmp_folder}/tmp.bai"
    subprocess.run(['touch', tmp_file])
    subprocess.run(["samtools", "view", "-b", "-o",tmp_file, bam_file, tmp_chrom])
    subprocess.run(["samtools", "index", tmp_file, tmp_index_file])
    tmp_bamfile = pysam.AlignmentFile(tmp_file, "rb")
    obs_iter = tmp_bamfile.fetch(tmp_chrom)
    obs_dfs = []
    for ix in obs_iter:
        obs_dfs.append(pd.DataFrame({'ID':[ix.to_dict()['name']],'chrom':[ix.to_dict()['ref_name']],'start':[ix.to_dict()['ref_pos']],'flag':[ix.to_dict()['flag']]}))
    subprocess.run(["rm", tmp_file, tmp_index_file])
    return pd.concat(obs_dfs)    

def perform_HDBScan_clustering(read_tbl,min_cluster,processes):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster,
                        metric='euclidean',
                        core_dist_n_jobs=processes)
    clusterer.fit(read_tbl.loc[:,['start']])
    return clusterer

def collect_hdb_cluster_read(g):
    
    leaves = set([v for v, d in g.out_degree() if d == 0])
    HDB_clusters = [v for v, d in g.out_degree() if d > 0]

    cl_read_idx = [list(nx.descendants(g,i).intersection(leaves)) for i in HDB_clusters]
    cl_read_tbl = pd.DataFrame({"HDB_cluster":HDB_clusters,"read_id_set":cl_read_idx})
    return(cl_read_tbl)

def build_cluster_data_df(tmp_cl_tbl,read_tbl):
    long_cl_read_tbl = tmp_cl_tbl.explode('read_id_set')
    long_cl_read_tbl = (long_cl_read_tbl
                        .assign(start = lambda df_:read_tbl.reset_index().start.to_numpy()[df_.read_id_set.to_numpy(dtype=int)])
                        )

    tmp_chr_hdb_summary_tbl = (long_cl_read_tbl
                                    .groupby('HDB_cluster')
                                    .agg(start = ('start','min'),
                                        end=('start','max'),
                                        rc = ('start','count')
                                        )
                                    .assign(w= lambda df_:df_.end - df_.start)
                                    .reset_index()
                                    .merge(tmp_cl_tbl)
                                )
    return tmp_chr_hdb_summary_tbl   


def build_hic_zscore(data_tbl,vars,vars_df,vars_degree,target_var):
    """
        Computes statistical values (z-scores and p-values) for a given dataset using a smoothing model.

        This function applies a Generalized Additive Model (GAM) with B-splines to analyze the relationship
        between selected variables and a target variable. It returns the original dataset with two additional
        columns: 
        - `zscore`: A measure of how much each data point deviates from the average.
        - `pvalue`: A value indicating the significance of each data point.

        Parameters:
        -----------
        data_tbl : pandas.DataFrame
            The dataset containing the variables to be analyzed.
        
        vars : list of str
            The names of the columns to be used as predictors (independent variables).
        
        vars_df : int
            The number of degrees of freedom for the spline smoothing.
        
        vars_degree : int
            The degree of the spline function, which controls how smooth the curve is.
        
        target_var : str
            The column name of the target variable (dependent variable) being analyzed.

        Returns:
        --------
        pandas.DataFrame
            The original dataset with two additional columns:
            - `zscore`: The standardized residuals (how different each value is from the model's prediction).
            - `pvalue`: A probability score indicating the significance of each data point.
        
        Example:
        --------
        Suppose `data_tbl` contains genetic data, where `vars` represents different genetic markers, 
        and `target_var` is a measure of gene expression. This function helps determine which data 
        points significantly deviate from the expected pattern.
    """
    x_spline = data_tbl[vars].to_numpy(dtype=float)
    y = data_tbl[target_var].to_numpy(dtype=float)
    bs = sm.gam.BSplines(x_spline, df=vars_df, degree=vars_degree)

    chr_gam = sm.GLMGam(y,smoother=bs)
    chr_gam_res = chr_gam.fit()
    gam_infl = chr_gam_res.get_influence()
    bs_tranform_exog = bs.transform(data_tbl[vars].to_numpy())
    tmp_rng = chr_gam_res.get_distribution(exog=bs_tranform_exog)
    mod_pvalue = tmp_rng.sf(data_tbl[target_var].to_numpy())
    new_data_tbl = data_tbl.assign(zscore = gam_infl.resid_studentized,pvalue = mod_pvalue)
    return new_data_tbl
