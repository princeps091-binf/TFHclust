# %%
import pandas as pd
import fast_hdbscan
import numpy as np
import bioframe as bf
import networkx as nx
import statsmodels.api as sm
import subprocess
import numba
import io

# %%
def get_chrom_read_tbl_from_bam(bam_file, tmp_chrom):
    # 1. On lance samtools view en extrayant uniquement les colonnes nécessaires :
    # QNAME (ID = %q), FLAG (%f), RNAME (chrom = %r), POS (start = %pos)
    # L'argument -@ 2 active le multi-threading pour accélérer la décompression
    cmd = ["samtools", "view", "-@", "2",  bam_file, tmp_chrom]
    
    # 2. Exécution du processus et capture du flux de sortie (stdout)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    
    # 3. Lecture directe du flux texte par Pandas (sans aucun fichier temporaire)
    # On spécifie les colonnes SAM standard : 0=ID, 1=FLAG, 2=CHROM, 3=START
    df = pd.read_csv(
        process.stdout, 
        sep='\t', 
        header=None, 
        usecols=[0, 1, 2, 3],
        names=['ID', 'flag', 'chrom', 'start'],
        dtype={'ID': str, 'flag': int, 'chrom': str, 'start': int}
    )
    
    # S'assurer que le processus s'est terminé correctement
    process.wait()
    
    # Réorganiser l'ordre des colonnes pour correspondre exactement à votre format d'origine
    return df[['ID', 'chrom', 'start', 'flag']]


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
    # Ensure degrees of freedom and degree parameters match the dimension of predictors if given as ints
    if isinstance(vars_df, int):
        vars_df = [vars_df] * x_spline.shape[1]
    if isinstance(vars_degree, int):
        vars_degree = [vars_degree] * x_spline.shape[1]

    bs = sm.gam.BSplines(x_spline, df=vars_df, degree=vars_degree)

    # 2. Configure Negative Binomial family (automatically uses log-link for counts)
    nb_family = sm.families.NegativeBinomial(link=sm.families.links.Log(), alpha=1.0)

    # 3. Fit GLMGam with Negative Binomial distribution
    # Note: GLMGam takes the endog (y) and optional linear exog, plus the smoother
    chr_gam = sm.GLMGam(endog=y, smoother=bs, family=nb_family)
    chr_gam_res = chr_gam.fit()

    # 4. Extract influence metrics and standardized residuals
    gam_infl = chr_gam_res.get_influence()
    
    # 5. Calculate empirical survival function (p-values) from the fitted Negative Binomial distribution
    # Transform exog to evaluate the fitted distribution per data point
    bs_transform_exog = bs.transform(x_spline)
    tmp_rng = chr_gam_res.get_distribution(exog=bs_transform_exog)
    mod_pvalue = tmp_rng.sf(y)

    # 6. Assign results back to the dataframe
    new_data_tbl = data_tbl.assign(
        zscore=gam_infl.resid_studentized,
        pvalue=mod_pvalue
    )
    
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
        self.cluster_read_coords_cache: dict[int, np.ndarray] = {}
    
    def HDBScan_clustering(self,min_cluster,processes):
        total_cores = numba.config.NUMBA_NUM_THREADS
        # 2. Ajuster selon le choix de l'utilisateur (ex: gestion du -1 classique)
        if processes == -1 or processes > total_cores:
            target_threads = total_cores
        elif processes <= 0:
            raise ValueError("processes must be greater than 0 or equal to -1")
        else:
            target_threads = processes
        # 3. Sauvegarder le nombre initial de threads pour pouvoir le restaurer après
        previous_threads = numba.get_num_threads()
        try:
            # 4. Assigner temporairement le nombre de cœurs choisi par l'utilisateur
            numba.set_num_threads(target_threads)
            clusterer = fast_hdbscan.HDBSCAN(min_cluster_size=min_cluster,
                            metric='euclidean')
            clusterer.fit(self.read_tbl.loc[:,['start']])
            self.full_graph = clusterer.condensed_tree_.to_networkx()
        finally:
            # 6. Toujours restaurer la configuration initiale (bonne pratique)
            numba.set_num_threads(previous_threads)

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
        # Build the optimized cluster-to-read-start coordinate cache for fast parallel task dispatch
        self.cluster_read_coords_cache = {
            int(cluster_id): group['start'].to_numpy(dtype=float)
            for cluster_id, group in long_cl_read_tbl.groupby('HDB_cluster')
        }
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
