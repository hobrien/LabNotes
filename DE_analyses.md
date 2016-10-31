

Created by [gh-md-toc](https://github.com/ekalinin/github-markdown-toc.go)


# Running DEseq2 via SARTools (at FDR < 0.05)
- When M vs. F Wald test run on all samples (N=108), there are 68 DE genes (mostly on sex chromosomes)
    - When samples < 14 weeks are excluded (N=61), there are 633 DE genes
    - When samples < 14 and >= 20 weeks are excluded (N=58), there are 1747 DE genes
        - When samples labeled A are excluded (N=53), there are 2513 DE genes
- When LRT is run on all samples not labeled A (N=100), using Sex * Age in the full model and dropping Sex in the reduced model, there are 
    - When samples < 14 and >= 20 weeks and labeled A are excluded (N=53), there are 868 DE genes