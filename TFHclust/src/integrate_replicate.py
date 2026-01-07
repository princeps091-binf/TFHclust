import hdbscan
import networkx as nx
import numpy as np
import bioframe as bf
import pandas as pd
import multiprocessing
from functools import partial
import scipy.stats as stats

def match_cluster(rep1_cl_read_tbl,rep2_cl_read_tbl):

    rep_overlap_df = (bf.overlap(rep1_cl_read_tbl.loc[:,['chrom','start','end','HDB_cluster']],
                                 rep2_cl_read_tbl.loc[:,['chrom','start','end','HDB_cluster']],
                                 return_overlap=True,
                                 how='inner',suffixes=['_rep1','_rep2'])
    .assign(inter_w=lambda df_:(df_.overlap_end-df_.overlap_start),
            end_point=lambda df_:df_[["end_rep1","end_rep2"]].values.tolist(),
            start_points = lambda df_:df_[["start_rep1","start_rep2"]].values.tolist())
    .assign(jaccard = lambda df_:df_.inter_w/(df_.end_point.apply(max) - df_.start_points.apply(min)))
    .drop(['end_point','start_points'],axis=1)
    )

    rep1_max_jaccard_idx = (rep_overlap_df
            .groupby(['chrom_rep1','start_rep1','end_rep1','HDB_cluster_rep1'])
            .jaccard.idxmax()
            )

    rep2_max_jaccard_idx = (rep_overlap_df
            .groupby(['chrom_rep2','start_rep2','end_rep2','HDB_cluster_rep2'])
            .jaccard.idxmax()
            )

    rep1_matching_tbl = rep_overlap_df.iloc[rep1_max_jaccard_idx,:]
    rep2_matching_tbl = rep_overlap_df.iloc[rep2_max_jaccard_idx,:]

    return (
        rep1_matching_tbl,
        rep2_matching_tbl
        )

def get_CvM_pvalue(i,tbl):
    return stats.cramervonmises_2samp(tbl.iloc[i].rep1_start,tbl.iloc[i].rep2_start).pvalue

def get_matched_cluster_agreement_pvalue(rep_match_tbl,rep1_clustering,rep2_clustering,njobs):
    tmp_rep1_read_start_tbl = (rep_match_tbl.loc[:,['HDB_cluster_rep1']]
    .drop_duplicates()
    .merge(rep1_clustering.cl_read_tbl.loc[:,['HDB_cluster','read_id_set']],left_on = ['HDB_cluster_rep1'],right_on = ['HDB_cluster'])
    .drop('HDB_cluster',axis=1)
    .assign(DNA_start = lambda df_: df_.apply(lambda x:rep1_clustering.read_tbl.start.iloc[x.read_id_set].to_numpy(),axis=1))
    .drop('read_id_set',axis=1)
    .rename(columns={'DNA_start':"rep1_start"})
    )

    tmp_rep2_read_start_tbl = (rep_match_tbl.loc[:,['HDB_cluster_rep2']]
    .drop_duplicates()
    .merge(rep2_clustering.cl_read_tbl.loc[:,['HDB_cluster','read_id_set']],left_on = ['HDB_cluster_rep2'],right_on = ['HDB_cluster'])
    .drop('HDB_cluster',axis=1)
    .assign(DNA_start = lambda df_: df_.apply(lambda x:rep2_clustering.read_tbl.start.iloc[x.read_id_set].to_numpy(),axis=1))
    .drop('read_id_set',axis=1)
    .rename(columns={'DNA_start':"rep2_start"})
    )

    rep_match_read_coord_tbl = (rep_match_tbl
    .loc[:,['HDB_cluster_rep1','HDB_cluster_rep2','jaccard']]
    .merge(tmp_rep1_read_start_tbl)
    .merge(tmp_rep2_read_start_tbl)
    )

    with multiprocessing.Pool(processes=njobs) as pool:
            # Using map_async method to perform square operation on all numbers parallely
            read_agreement_pvalue = pool.map(partial(get_CvM_pvalue, tbl = rep_match_read_coord_tbl),
                                            range(rep_match_read_coord_tbl.shape[0])) 
    return rep_match_read_coord_tbl.assign(cvm_pvalue = read_agreement_pvalue).loc[:,['HDB_cluster_rep1','HDB_cluster_rep2','jaccard','cvm_pvalue']]   

       
class merged_clustering:
    def __init__(self,rep1_clustering,rep2_clustering):
        self.rep1 = rep1_clustering
        self.rep2 = rep2_clustering
        self.chrom = self.rep1.read_tbl.iloc[0,:].chrom
        self.rep1_match_tbl = None
        self.rep2_match_tbl = None
        self.robust_cluster_coord_tbl = None
        self.inter_rep_match_tbl = None
        
    def match_rep_cluster(self):
        self.rep1_match_tbl, self.rep2_match_tbl = match_cluster(self.rep1.summary_tbl,self.rep2.summary_tbl)        

    def evaluate_overlap_significance(self,njobs):
            self.rep1_match_pvalue_tbl = get_matched_cluster_agreement_pvalue(self.rep1_match_tbl,self.rep1,self.rep2,njobs)
            self.rep2_match_pvalue_tbl = get_matched_cluster_agreement_pvalue(self.rep2_match_tbl,self.rep1,self.rep2,njobs)


    def filter_significant_overlap_cluster(self,CvM_thresh,specificity_thresh):
        
        rep1_matched_cluster_ID_list = self.rep1_match_pvalue_tbl.query("cvm_pvalue > @CvM_thresh").HDB_cluster_rep1.drop_duplicates().to_list()
        rep2_matched_cluster_ID_list = self.rep2_match_pvalue_tbl.query("cvm_pvalue > @CvM_thresh").HDB_cluster_rep2.drop_duplicates().to_list()

        rep1_robust_cluster_tbl = (self.rep1
                                    .summary_tbl.loc[:,['HDB_cluster','chrom','start','end','zscore','pvalue']]
                                    .query("HDB_cluster in @rep1_matched_cluster_ID_list")
                                    .assign(rep = "rep1")
                                    .merge(self.rep1_match_pvalue_tbl,left_on='HDB_cluster',right_on="HDB_cluster_rep1")
                                    .drop('HDB_cluster',axis=1)
                                    .query("HDB_cluster_rep2 in @rep2_matched_cluster_ID_list")
                                    .merge(self.rep2.summary_tbl.loc[:,['HDB_cluster','pvalue']].rename(columns={'HDB_cluster':'HDB_cluster_rep2','pvalue':'pvalue_rep2'}))
                                    .query('pvalue < @specificity_thresh and pvalue_rep2 < @specificity_thresh')
                                    .loc[:,['chrom','start','end','zscore','pvalue','rep','HDB_cluster_rep1','HDB_cluster_rep2','jaccard','cvm_pvalue']]
                                )
        rep2_robust_cluster_tbl = (self.rep2
                                    .summary_tbl.loc[:,['HDB_cluster','chrom','start','end','zscore','pvalue']]
                                    .query("HDB_cluster in @rep2_matched_cluster_ID_list")
                                    .assign(rep = "rep2")
                                    .merge(self.rep2_match_pvalue_tbl,left_on='HDB_cluster',right_on="HDB_cluster_rep2")
                                    .drop('HDB_cluster',axis=1)
                                    .query("HDB_cluster_rep1 in @rep1_matched_cluster_ID_list")
                                    .merge(self.rep1.summary_tbl.loc[:,['HDB_cluster','pvalue']].rename(columns={'HDB_cluster':'HDB_cluster_rep1','pvalue':'pvalue_rep1'}))
                                    .query('pvalue < @specificity_thresh and pvalue_rep1 < @specificity_thresh')
                                    .loc[:,['chrom','start','end','zscore','pvalue','rep','HDB_cluster_rep2','HDB_cluster_rep1','jaccard','cvm_pvalue']]
                                )
        if (rep1_robust_cluster_tbl.shape[0] > 0 and rep2_robust_cluster_tbl.shape[0]):
            self.robust_cluster_match_tbl = bf.cluster(pd.concat([rep1_robust_cluster_tbl.drop('HDB_cluster_rep2',axis=1).rename(columns={'HDB_cluster_rep1':'HDB_cluster'}),
                                                                  rep2_robust_cluster_tbl.drop('HDB_cluster_rep1',axis=1).rename(columns={'HDB_cluster_rep2':'HDB_cluster'})])).assign(w = lambda df_: df_.cluster_end -df_.start).sort_values('w')

            self.robust_cluster_coord_tbl = (self.robust_cluster_match_tbl
                            .loc[:,['chrom','cluster_start','cluster_end','cluster']]
                            .drop_duplicates()
                            .rename(columns = {'cluster_start':'start','cluster_end':"end"})
                            .assign(w = lambda df_: df_.end - df_.start)
                            )
            self.inter_rep_match_tbl = (pd.concat([rep2_robust_cluster_tbl
                                        .loc[:,['HDB_cluster_rep1','HDB_cluster_rep2','cvm_pvalue','jaccard']]
                                        .merge(self.rep1.summary_tbl.loc[:,['HDB_cluster','chrom','start','end','zscore','pvalue','norm_lvl']]
                                               .rename(columns = {'start':'start_rep1','end':'end_rep1','pvalue':'pvalue_rep1','zscore':'zscore_rep1','norm_lvl':'norm_lvl_rep1'}),left_on = 'HDB_cluster_rep1',right_on='HDB_cluster').drop('HDB_cluster',axis=1)
                                        .merge(self.rep2.summary_tbl.loc[:,['HDB_cluster','start','end','zscore','pvalue','norm_lvl']]
                                               .rename(columns = {'start':'start_rep2','end':'end_rep2','pvalue':'pvalue_rep2','zscore':'zscore_rep2','norm_lvl':'norm_lvl_rep2'}),left_on = 'HDB_cluster_rep2',right_on='HDB_cluster').drop('HDB_cluster',axis=1)
                                        ,
                                        rep1_robust_cluster_tbl
                                        .loc[:,['HDB_cluster_rep1','HDB_cluster_rep2','cvm_pvalue','jaccard']]
                                        .merge(self.rep1.summary_tbl.loc[:,['HDB_cluster','chrom','start','end','zscore','pvalue','norm_lvl']]
                                               .rename(columns = {'start':'start_rep1','end':'end_rep1','pvalue':'pvalue_rep1','zscore':'zscore_rep1','norm_lvl':'norm_lvl_rep1'}),left_on = 'HDB_cluster_rep1',right_on='HDB_cluster').drop('HDB_cluster',axis=1)
                                        .merge(self.rep2.summary_tbl.loc[:,['HDB_cluster','start','end','zscore','pvalue','norm_lvl']]
                                               .rename(columns = {'start':'start_rep2','end':'end_rep2','pvalue':'pvalue_rep2','zscore':'zscore_rep2','norm_lvl':'norm_lvl_rep2'}),left_on = 'HDB_cluster_rep2',right_on='HDB_cluster').drop('HDB_cluster',axis=1)
                                        ]).drop_duplicates()
                                        .assign(start = lambda df_: df_.apply(lambda row: min(row.start_rep1,row.start_rep2),axis=1),
                                                end = lambda df_: df_.apply(lambda row: max(row.end_rep1,row.end_rep2),axis=1))
                                        .assign(w = lambda df_: df_.end -df_.start))
        else:
            print("no significant and replicable cluster")



    