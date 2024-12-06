# Estimate the Treatment Effect of Simulated Data on an Outcome Variable Y

### Background
Cluster randomized trials (CRTs) are commonly used to evaluate interventions applied at the group level, such as in schools or clinics. Designing CRTs involves balancing the number of clusters and the number of observations per cluster to optimize statistical efficiency while staying within budget constraints. These decisions are further complicated by factors such as intra-cluster correlation and differing costs for initial and additional samples in each cluster. Simulation studies provide a powerful tool to systematically evaluate design options under varying budget constraints and data-generating conditions, allowing researchers to identify optimal designs for estimating treatment effects accurately and efficiently.

### Methods
In this study, we use a simulation-based approach to evaluate the design of cluster randomized trials (CRTs) under fixed budget constraints. The hierarchical model assumes that the outcome variable (Y) is influenced by fixed effects (treatment assignment) and random effects (cluster-level variability). We systematically vary key parameters, including the treatment effect size ($\beta$), cluster-level variance ($\gamma^2$), and within-cluster variance ($\sigma^2$), while considering realistic cost structures where the first sample in a cluster (c1) is more expensive than subsequent samples (c2). Performance metrics, such as bias, standard error, and coverage probability, are calculated to assess the accuracy and precision of the treatment effect estimates. Additionally, the study explores how varying the cost ratio (c1/c2) and data generation parameters impact the optimal allocation of resources.

### Results
Our simulation study evaluated the impact of underlying data generation parameters and cost ratios on the estimation of treatment effects. The results demonstrated that bias remained minimal across all parameter settings, while standard errors increased with larger $\beta$ and $\gamma^2$, reflecting heightened variability in these conditions. Coverage rates were consistent across most settings but slightly reduced under higher $\alpha^2$ values. Cost ratios significantly influenced the allocation of clusters and observations per cluster, with lower c1/c2 ratios favoring smaller clusters and more observations, resulting in reduced bias and variability.



## Files

### R 

`run.R`: Contains the code used for generating data. 

`generate.R`:  Contains all functions written for simulation study. 


### report

`simulation_report.Rmd`: The Rmarkdown version of the  report, which includes both written text interpretations and raw code used in the analysis. 

`simulation_report.pdf`: The PDF version of the report, which includes both written text interpretations and a Code Appendix with the raw code used in the analysis. 


## Dependencies

The following packages were used in this analysis: 

 - Analysis: `lmer`
 - Table Formatting: `gtsummary`, `knitr`, `kableExtra`
 - Data Visualization: `ggplot2`, `gridExtra`


