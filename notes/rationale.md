# TFHclust

- Adapt HDBScan to ChIP data at the single read level
    - merit:
        - hierarchical description of read aggregation along the genome
            - different scales of aggregation may highlight different processes determining the TF residency along the genome
        - more explicit and adaptive correspondence with control data using GAM
- Method overview:
    1. Per chromosome:
        - Perform HDBScan on single chromosome to describe read aggregation from whole chromosome to read triad
        - Estimate significance using GAM to predict observed read count in a given HDBScan cluster from the corresponding control read count
        - highlight culled subtree with significant aggregation content
        - highlight effective density by examining "local" enrichment: observed number reads within window spanning twice the width of the considered cluster (Poisson test) + chromosome-wide rate -> pick most conservative (MACS method)
- Payload:
    - How do different scale determine one another
        - What kind of regions get zoned in upon by increasingly smaller scale aggregations -> within scope of particular enriched subtree
            - TFBS content configuration as we consider different scales of read aggregate
                - how do successive scale of the same structure consolidate the enrichment in particular kinds of TFBS
        - "Fragmentation" pattern
    - Highlight aggregate event that are significant but missed by conventional peak detection
        - peak-less clusters -> low-affinity events
        - More broadly what is the broader and still significant affinity context beyond conventional peaks
    - Relation to known regulatory intervals
        - CRE
        - Gene annotation
## Viz/Fig
1. Ctrl vs Obs. + infered cluster tracks
    - incorporate hierarchy of signficant subtrees
2. Correspondence between Obs and Ctrl
3. relation between peaks and broader affinity landscape
    - Different embedding
4. Comparative analysis to illustrate change in landscape
    - same cell line but different factors
    - different cell lines but same factor
5. TFBS enrichment pattern
