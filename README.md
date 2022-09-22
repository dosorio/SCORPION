# Population-level comparisons of gene regulatory networks modeled on high-throughput single-cell transcriptomic data.
Single-cell technologies enable high-resolution studies of the molecular mechanisms that define phenotypes. 
However, modeling biological variability across multiple samples characterized using single-cell technologies for differential analyses is difficult due to data sparsity and cellular heterogeneity. 
As a result, during single-cell gene regulatory network analyses, cell transcriptomes from different samples are typically collapsed by experimental group and then used to construct a single network representing it. 
This network is then interrogated or compared to other networks to gain insight into the transcription factor-target gene interactions underlying the phenotype of interest. 
We present **SCORPION**, a tool that reconstructs gene regulatory networks from single cells/nuclei RNA-seq data using the same baseline priors in a network refinement approach based on message passing. 
This method allows for the reconstruction of comparable networks suitable for use in population-level studies. 
We tested **SCORPION** using synthetic data and found that it outperforms 12 other gene regulatory network reconstruction techniques. 
Additionally, using supervised experiments, we show that **SCORPION** can accurately identify biological differences in regulatory networks between wild-type cells and cells carrying transcription factors perturbations.
Furthermore, we demonstrate **SCORPION**'s scalability to population-level analyses by applying it to a single-cell RNA-seq atlas that includes 200,436 cells derived from different regions of colorectal cancer tumors and healthy adjacent tissues. 
We performed a comparative network analysis on these networks, which detected regulatory differences between healthy adjacent tissue and different tumor regions that are consistent with our understanding of disease progression, and networks between left-sided and right-sided tumors that provide insight into the regulators associated with the phenotypes and the differences in their survival rate.
