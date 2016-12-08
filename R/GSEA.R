# GSEA 1.0 -- Gene Set Enrichment Analysis / Broad Institute 
#
# R script to run GSEA Analysis of the Leukemia ALL/AML vs C1 example (cut and paste into R console)

GSEA.program.location <- "/Users/heo3/BTSync/Code/GSEA-P-R/GSEA.1.0.R"   #  R source program (change pathname to the rigth location in local machine)
source(GSEA.program.location, verbose=T, max.deparse.length=9999)
doc_string = "MvsF_C3_tft"
output_directory <- paste0("/Users/heo3/BTSync/FetalRNAseq/GSEA/", doc_string, "/")
dir.create(output_directory)
GSEA(                                                                    # Input/Output Files :-------------------------------------------
 input.ds =  "/Users/heo3/BTSync/FetalRNAseq/Counts/MvsF_14_20_noA_excl_15641_18432_Cooks.75/tables/MvsF.gct",           # Input gene expression Affy dataset file in RES or GCT format
 input.cls = "/Users/heo3/BTSync/FetalRNAseq/Counts/MvsF_14_20_noA_excl_15641_18432_Cooks.75/tables/MvsF.cls",           # Input class vector (phenotype) file in CLS format
 gs.db =     "/Users/heo3/BTSync/FetalRNAseq/GSEA/MSigDB/c3.tft.v5.2.symbols.gmt",         # Gene set database in GMT format
 output.directory      = output_directory,        # Directory where to store output and results (default: "")
#  Program parameters :-------------------------------------------------------------------------------------------------------------------------
 doc.string            = doc_string,   # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
 non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
 reshuffling.type      = "sample.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
 nperm                 = 1000,            # Number of random permutations (default: 1000)
 weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
 nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
 fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
 fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
 topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
 adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
 gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
 gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
 reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
 preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
 random.seed           = 3338,            # Random number generator seed. (default: 123456)
 perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
 fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
 replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
 save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
 OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
 use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
)
#-----------------------------------------------------------------------------------------------------------------------------------------------

# Overlap and leading gene subset assignment analysis of the GSEA results

GSEA.Analyze.Sets(
   directory  = output_directory,        # Directory where to store output and results (default: "")
   topgs = 20,                                                           # number of top scoring gene sets used for analysis
   height = 16,
   width = 16
)
library(readr)
SUMMARY_REPORT_Female <- read_delim(paste0("~/BTSync/FetalRNAseq/GSEA/",
                                           doc_string,
                                           "/",
                                           doc_string,
                                           ".SUMMARY.RESULTS.REPORT.Female.txt"
                                           ),
                                    "\t", 
                                    escape_double = FALSE, 
                                    trim_ws = TRUE
                                  )
                                                    
SUMMARY_REPORT_Male <- read_delim(paste0("~/BTSync/FetalRNAseq/GSEA/",
                                           doc_string,
                                           "/",
                                           doc_string,
                                           ".SUMMARY.RESULTS.REPORT.Male.txt"
                                           ),
                                    "\t",
                                    escape_double = FALSE, 
                                    trim_ws = TRUE
                                    )
