import fast_hdbscan
import networkx as nx
import numpy as np
import bioframe as bf
import pandas as pd
import multiprocessing
from functools import partial
import scipy.stats as stats
import concurrent.futures
import os

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


#------------------------------

class DirectedMultiRepGraphClustering:
    """Builds a multi-replicate consensus graph using directed pairwise best-match edges, safely handling nested intervals."""

    def __init__(self, clustering_objs_dict: dict[str, object]) -> None:
        """
        clustering_objs_dict: Dict mapping replicate names to their fitted ReadClustering objects.
        e.g., {'rep1': obj1, 'rep2': obj2, 'rep3': obj3}
        """
        self.replicates = clustering_objs_dict
        self.graph = nx.DiGraph()  # Directed graph to capture directional best-match relationships
        self.consensus_loci_tbl: pd.DataFrame = None

    def _compute_pairwise_best_matches(self, rep_name_a: str, rep_name_b: str, df_a: pd.DataFrame, df_b: pd.DataFrame) -> pd.DataFrame:
        """Computes directional best matches specifically from Replicate A -> Replicate B."""
        overlaps = bf.overlap(
            df_a.loc[:, ['chrom', 'start', 'end', 'HDB_cluster']],
            df_b.loc[:, ['chrom', 'start', 'end', 'HDB_cluster']],
            return_overlap=True,
            how='inner',
            suffixes=['_a', '_b']
        )
        
        if overlaps.empty:
            return pd.DataFrame()
        # Compute Jaccard overlap
        overlaps = (
            overlaps
            .assign(
                inter_w=lambda df_: df_['overlap_end'] - df_['overlap_start'],
                max_end=lambda df_: np.maximum(df_['end_a'], df_['end_b']),
                min_start=lambda df_: np.minimum(df_['start_a'], df_['start_b'])
            )
            .assign(jaccard=lambda df_: df_['inter_w'] / (df_['max_end'] - df_['min_start']))
        )
        # Find best match from A -> B (for each cluster in A, find max Jaccard in B)
        best_a_to_b_idx = overlaps.groupby(['HDB_cluster_a'])['jaccard'].idxmax()
        match_a_b = overlaps.loc[best_a_to_b_idx].copy()
        match_a_b['direction'] = f"{rep_name_a}->{rep_name_b}"
        return match_a_b

    def build_graph(self,jaccard_thresh: float) -> None:
        """Iterates through all ordered replicate pairs, computes directed best matches, and populates the graph."""
        rep_names = list(self.replicates.keys())
        
        # 1. Add all cluster nodes to the graph
        for rep_name, obj in self.replicates.items():
            if obj.summary_tbl is None:
                raise RuntimeError(f"Replicate '{rep_name}' has no summary_tbl.")
            for _, row in obj.summary_tbl.iterrows():
                node_id = f"{rep_name}_cl_{row['HDB_cluster']}"
                self.graph.add_node(node_id, rep=rep_name, **row.to_dict())

        # 2. Evaluate all ordered pairs of replicates exhaustively (A -> B and B -> A separately)
        for name_a in rep_names:
            for name_b in rep_names:
                if name_a == name_b:
                    continue
                
                obj_a = self.replicates[name_a]
                obj_b = self.replicates[name_b]

                match_df = self._compute_pairwise_best_matches(name_a, name_b, obj_a.summary_tbl, obj_b.summary_tbl)
                
                if match_df.empty:
                    continue

                # 3. Add clean directed edges: Source is always in Replicate A, Destination is always in Replicate B
                for _, row in match_df.iterrows():
                    if row['jaccard'] < jaccard_thresh:
                        continue

                    src = f"{name_a}_cl_{row['HDB_cluster_a']}"
                    dst = f"{name_b}_cl_{row['HDB_cluster_b']}"

                    self.graph.add_edge(src, dst, jaccard=row['jaccard'], direction=row['direction'])
    def get_edge_table(self) -> pd.DataFrame:
        """
        Extracts all directed edges, Jaccard scores, and metadata from the graph 
        into a clean Pandas DataFrame.
        """
        if self.graph.number_of_edges() == 0:
            return pd.DataFrame(columns=['source', 'target', 'jaccard', 'direction', 'source_rep', 'target_rep'])
        edge_records = []
        for u, v, data in self.graph.edges(data=True):
            edge_records.append({
                "source": u,
                "target": v,
                "jaccard": data.get("jaccard"),
                'cvm_pvalue':data.get('cvm_pvalue'),
                "direction": data.get("direction"),
                "source_rep": self.graph.nodes[u].get("rep"),
                "target_rep": self.graph.nodes[v].get("rep")
            })
        return pd.DataFrame(edge_records)
    def extract_consensus_components(self, jaccard_thresh: float = 0.5, cvm_pval_thresh: float = 0.5) -> pd.DataFrame:
        """Extracts connected components from the directed best-match graph (treated as undirected for grouping)."""
        # Convert to undirected for component grouping (transitive connectivity across replicates)
        # 1. Create a working copy or subgraph with only high-quality edges that pass both thresholds
        valid_edges = [
            (u, v) for u, v, data in self.graph.edges(data=True)
            if data.get('jaccard', 0.0) >= jaccard_thresh 
            and data.get('cvm_pvalue', 1.0) >= cvm_pval_thresh
        ]
        
        # Build an undirected subgraph containing only nodes and edges that meet criteria
        filtered_g = self.graph.edge_subgraph(valid_edges).copy()
        # 3. Convert the directed subgraph to an undirected graph for transitive component grouping
        undirected_g = filtered_g.to_undirected()
        consensus_records = []
        for comp_id, component in enumerate(nx.connected_components(undirected_g)):        
            subgraph_data = [self.graph.nodes[node] for node in component]
            chrom = subgraph_data[0]['chrom']
            starts = [n['start'] for n in subgraph_data]
            ends = [n['end'] for n in subgraph_data]
            replicates_present = {n['rep'] for n in subgraph_data}
            consensus_records.append({
                "consensus_id": comp_id,
                "chrom": chrom,
                "start": int(min(starts)),
                "end": int(max(ends)),
                "replicate_support": len(replicates_present),
                "replicates": ",".join(sorted(replicates_present)),
                "cluster_nodes": list(component)
            })
        self.consensus_loci_tbl = pd.DataFrame(consensus_records).query('replicate_support > 1')
        return self.consensus_loci_tbl
    
    def compute_cvm_edge_weights(self) -> None:
        """
        Computes the two-sample Cramér-von Mises test for all existing edges using 
        the asymptotic limiting distribution, restricted to the shared overlap window.
        """
        for u, v, data in self.graph.edges(data=True):
            rep_u = self.graph.nodes[u]['rep']
            rep_v = self.graph.nodes[v]['rep']
            
            cluster_u_id = self.graph.nodes[u]['HDB_cluster']
            cluster_v_id = self.graph.nodes[v]['HDB_cluster']
            
            obj_u = self.replicates[rep_u]
            obj_v = self.replicates[rep_v]
            
            reads_u = obj_u.read_tbl.iloc[obj_u.cl_read_tbl.query("HDB_cluster == @cluster_u_id").read_id_set.to_numpy()[0],:].start.to_numpy()            
            reads_v = obj_v.read_tbl.iloc[obj_v.cl_read_tbl.query("HDB_cluster == @cluster_v_id").read_id_set.to_numpy()[0],:].start.to_numpy()
            
            if len(reads_u) < 2 or len(reads_v) < 2:
                self.graph[u][v]['cvm_statistic'] = np.nan
                self.graph[u][v]['cvm_pvalue'] = 1.0
                continue
                
            start_u, end_u = self.graph.nodes[u]['start'], self.graph.nodes[u]['end']
            start_v, end_v = self.graph.nodes[v]['start'], self.graph.nodes[v]['end']
            # 1. Define the union window (from global min start to global max end)
            union_start = min(start_u, start_v)
            union_end = max(end_u, end_v)
            union_span = union_end - union_start
            
            if union_span <= 0:
                self.graph[u][v]['cvm_statistic'] = np.nan
                self.graph[u][v]['cvm_pvalue'] = 1.0
                continue
                
            # 2. Retain reads falling within the union coordinate frame
            reads_u_filtered = reads_u[(reads_u >= union_start) & (reads_u <= union_end)]
            reads_v_filtered = reads_v[(reads_v >= union_start) & (reads_v <= union_end)]
            
            if len(reads_u_filtered) < 2 or len(reads_v_filtered) < 2:
                self.graph[u][v]['cvm_statistic'] = np.nan
                self.graph[u][v]['cvm_pvalue'] = 1.0
                continue
                
            # 3. Normalize read positions to [0, 1] relative to the union span
            norm_u = np.clip((reads_u_filtered - union_start) / union_span, 0.0, 1.0)
            norm_v = np.clip((reads_v_filtered - union_start) / union_span, 0.0, 1.0)
            
            # 4. Compute CvM test using the asymptotic limiting distribution
            res = stats.cramervonmises_2samp(norm_u, norm_v, method='asymptotic')            
            # Assign metrics
            self.graph[u][v]['cvm_statistic'] = res.statistic
            self.graph[u][v]['cvm_pvalue'] = res.pvalue   
    @staticmethod
    def _worker_compute_cvm(u, v, node_u_data, node_v_data, reads_u, reads_v):
        """
        Static worker function running in a separate process. 
        Takes only primitive arrays and dictionaries to minimize IPC overhead.
        """
        try:
            if len(reads_u) < 2 or len(reads_v) < 2:
                return u, v, np.nan, 1.0
                
            start_u, end_u = node_u_data['start'], node_u_data['end']
            start_v, end_v = node_v_data['start'], node_v_data['end']
            
            # Union window span
            union_start = min(start_u, start_v)
            union_end = max(end_u, end_v)
            union_span = union_end - union_start
            
            if union_span <= 0:
                return u, v, np.nan, 1.0
                
            # Filter reads to union span
            reads_u_filtered = reads_u[(reads_u >= union_start) & (reads_u <= union_end)]
            reads_v_filtered = reads_v[(reads_v >= union_start) & (reads_v <= union_end)]
            
            if len(reads_u_filtered) < 2 or len(reads_v_filtered) < 2:
                return u, v, np.nan, 1.0
                
            # Normalize to [0, 1] relative to union span
            norm_u = np.clip((reads_u_filtered - union_start) / union_span, 0.0, 1.0)
            norm_v = np.clip((reads_v_filtered - union_start) / union_span, 0.0, 1.0)
            
            # Asymptotic CvM test
            res = stats.cramervonmises_2samp(norm_u, norm_v, method='asymptotic')
            return u, v, res.statistic, res.pvalue
            
        except Exception:
            return u, v, np.nan, 1.0

    def parallel_compute_cvm_edge_weights(self, max_workers: int = None) -> None:
        """
        Computes CvM edge weights in parallel across CPU cores while minimizing memory footprint
        by pre-fetching only required cluster read arrays in the main process.
        """
        if max_workers is None:
            max_workers = max(1, os.cpu_count() - 1)  # Leave one core free for system stability

        edge_tasks = []
        
        # Pre-gather tasks and extract read arrays sequentially in main thread 
        # to avoid race conditions or locking issues on replicate objects.
        for u, v, data in self.graph.edges(data=True):
            rep_u = self.graph.nodes[u]['rep']
            rep_v = self.graph.nodes[v]['rep']
            
            cluster_u_id = self.graph.nodes[u]['HDB_cluster']
            cluster_v_id = self.graph.nodes[v]['HDB_cluster']
            
            obj_u = self.replicates[rep_u]
            obj_v = self.replicates[rep_v]
            
            # Extract target reads securely
            reads_u = obj_u.read_tbl.iloc[obj_u.cl_read_tbl.query("HDB_cluster == @cluster_u_id").read_id_set.to_numpy()[0],:].start.to_numpy()            
            reads_v = obj_v.read_tbl.iloc[obj_v.cl_read_tbl.query("HDB_cluster == @cluster_v_id").read_id_set.to_numpy()[0],:].start.to_numpy()
            
            node_u_data = self.graph.nodes[u]
            node_v_data = self.graph.nodes[v]
            
            edge_tasks.append((u, v, node_u_data, node_v_data, reads_u, reads_v))

        if not edge_tasks:
            return

        # Execute in parallel using ProcessPoolExecutor
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(self._worker_compute_cvm, u, v, n_u, n_v, r_u, r_v)
                for u, v, n_u, n_v, r_u, r_v in edge_tasks
            ]
            
            for future in concurrent.futures.as_completed(futures):
                u, v, stat, pval = future.result()
                self.graph[u][v]['cvm_statistic'] = stat
                self.graph[u][v]['cvm_pvalue'] = pval

