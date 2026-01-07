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

class read_clustering:
    def __init__(self,read_df):
        self.read_tbl = read_df
        self.chrom = read_df.iloc[0,:].chrom
        self.full_graph = None
        self.cl_read_tbl = None
        self.regression_tbl = None
        self.res_tbl = None
        self.summary_tbl = None
    
    def HDBScan_clustering(self,min_size,njobs):
        clusterer = hdbscan.HDBSCAN(min_cluster_size=min_size,
                        metric='euclidean',
                        core_dist_n_jobs=njobs)
        clusterer.fit(self.read_tbl.loc[:,['start']])
        self.full_graph = clusterer.condensed_tree_.to_networkx()
    
    def collect_hdb_cluster_read(self):
    
        leaves = set([v for v, d in self.full_graph.out_degree() if d == 0])
        HDB_clusters = [v for v, d in self.full_graph.out_degree() if d > 0]

        cl_read_idx = [list(nx.descendants(self.full_graph,i).intersection(leaves)) for i in HDB_clusters]
        tmp_cl_read_tbl = pd.DataFrame({"HDB_cluster":HDB_clusters,"read_id_set":cl_read_idx})
        node_depth = nx.shortest_path_length(self.full_graph,source=tmp_cl_read_tbl.HDB_cluster.min())
        tmp_cl_read_tbl = (tmp_cl_read_tbl
                           .assign(lvl = lambda df_: [node_depth[i] for i in df_.HDB_cluster.to_list()])
                           .assign(norm_lvl = lambda df_: df_.lvl/df_.lvl.max()))

        long_cl_read_tbl = tmp_cl_read_tbl.explode('read_id_set')

        long_cl_read_tbl = (long_cl_read_tbl
                            .assign(start = lambda df_:self.read_tbl.reset_index().start.to_numpy()[df_.read_id_set.to_numpy(dtype=int)])
                            )

        self.cl_read_tbl = (long_cl_read_tbl
                                        .groupby('HDB_cluster')
                                        .agg(start = ('start','min'),
                                            end=('start','max'),
                                            rc = ('start','count')
                                            )
                                        .assign(w= lambda df_:df_.end - df_.start)
                                    ).reset_index().merge(tmp_cl_read_tbl)

    def  build_regression_tbl(self,ctrl_read_df):

        bg_coord_df = ctrl_read_df.loc[:,['chrom','start']].assign(start = lambda df_: df_.start.astype(int)).assign(end = lambda df_: df_.start + 1)
        hdb_cluster_bg_rc_df = bf.count_overlaps(self.cl_read_tbl.assign(chrom = self.chrom).loc[:,['chrom','start','end','HDB_cluster']],bg_coord_df)
        tmp_chr_hdb_summary_tbl = self.cl_read_tbl.merge(hdb_cluster_bg_rc_df.rename(columns={'count':'bg_count'}))
        
        self.regression_tbl = (tmp_chr_hdb_summary_tbl
                    .loc[:,['HDB_cluster','rc','bg_count']]
                    .assign(lrc = lambda df_: np.log10(df_.rc), lbg = lambda df_: np.log10(df_.bg_count))
                    .assign(lbg2 = lambda df_: np.where(df_.bg_count.lt(1),0,df_.lbg))
                    )
    def compute_zscore(self,explanatory_vars,vars_df,vard_degree,target_variable):
        mod_res_tbl = build_hic_zscore(self.regression_tbl,explanatory_vars,vars_df,vard_degree,target_variable)
        self.res_tbl = self.regression_tbl.merge(mod_res_tbl)
        self.summary_tbl =   (self.res_tbl
                                .merge(self.cl_read_tbl.assign(chrom = self.chrom).loc[:,['HDB_cluster','chrom','start','end','w','norm_lvl']])
                                .assign(FC = lambda df_: (df_.rc/df_.w)/((df_.bg_count+1)/df_.w))
                             )
