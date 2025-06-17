---
marp: true
paginate: true
---

# ChIP-ARC: ChIP **A**daptive **R**ead **C**lustering

---

1. Chromatin binding/patterning is multi-scale
2. HDBScan to describe multi-scale aggregation
3. Preliminary results

---

![bg contain](./img/transcriptionalregulation.png)

---

## ChIPseq data

---

![bg contain](./img/CTCF_chip_rep1_vs_rep2.png)

---

## Read-level clustering

---
![bg contain](./img/HDBSCan.png)

---

## Leveraging controls to estimate significance

---

![bg contain](./img/CTCF_ctrl_vs_obs.png)

---
### A model to correct for confounding factors

- Model describing how confounding can predict observed read density
    - $obs.count_{i} \sim control.count_{i}$
- GAM regression
    - direction of "error" will indicate denser or sparser interaction than expected


---

![bg contain](./img/CTCF_ctrl_vs_obs_zscore.png)


---

## Cluster Replicability

---

### Quantifying cluster convergence
- Similar regions -> Jaccard
- Similar shape -> Cramer von Mises

---

![bg contain](./img/CTCF_zscore_replicability.png)

---

![bg contain](./img/CTCF_depth_replicability.png)

---

### Illustrative results

---
![bg contain](./img/result_track.png)

---

## aggregation along scale

---
### CTCF aggregates span different scales
![bg right contain](./img/CTCF_size_vs_depth.png)

---

### H3K27me3 aggregates remain at the same scale
![bg right contain](./img/H3K27me3_size_vs_depth.png)

---

## Perspective

- How do conventional TFBS populate this broader binding landscape ?
- Are there other sequence determinants shaping these novel aggregates ?
- How constitutive is the balance between narrow and broad binding across factors and samples?