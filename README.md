# oracle

This repositorry contain a method to compare metagenomics and metabarcoding.

## Workflow

![workflow](https://github.com/GhMathis/oracle/blob/main/image/workflow_metaB_metaG.png)


## To Do

Remaining tasks, clarifications, or questions:

- [x] make a flowchart representing our statistical analyses, 
- [ ] rarefaction curves: some curves do not reach the rarefaction threshold (vertical line). Why?
- [ ] there is a large difference (100x) in the number of assigned reads for MetaG R2 and MetaG R1. Knowing that R1 and R2 have identical number of reads at the start of the analysis, why is there much more R2 reads assigned?
- [ ] diversity indices are very different between MetaG and MetaB. Hypothesis?
- [ ] beta diversity: identify the outliers (samples `ST-10`, `YL-19`, `TS-3`, `SER-10`, more?). Are these samples the ones that were re-sequenced?
- [ ] simulate the effect of amplification on mean and variance of our species throughout many samples. The goal is to compare with our observations?
- [ ] metacoder: very large fold-change values for non-parametric tests. Why?
