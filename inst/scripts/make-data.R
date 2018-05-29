###############################################################################
# data/exSCE.rda
###############################################################################
expr <- lapply(seq_len(10), function(i){
  set.seed(1101 + i)
  simulate_exprs(n_features=25,
                 n_samples=10,
                 prefix_sample=paste0("C", i, "_"))})
expr <- do.call(cbind, expr)
meta <- rep(c("2.Mid", "2.Mid", "3.Late", "1.Early", NA,
              "3.Late", "3.Late", NA, NA, "2.Mid"), each = 10)
exSCE <- SingleCellExperiment(assays=list(logcounts=expr),
                              colData=data.frame(age=meta))
se <- embedSamples(exSCE)
d <- findSpectrum(se$eigenvalues, frac=30)
latentSpace(exSCE) <- se$components[, d]
states(exSCE) <- findStates(exSCE,
                            min_size=0.01,
                            min_feat=2,
                            max_pval=1e-4,
                            min_fc=2)
exSCE <- connectStates(exSCE, l=30)
exSCE <- fitTrajectory(exSCE)
fn <- system.file("exdata/exSCE.graphml", package="CellTrails")
trajLayout(exSCE, adjust=TRUE) <- read.ygraphml(fn)
exSCE <- addTrail(exSCE, from="H1", to="H2", name="Tr1")
exSCE <- addTrail(exSCE, from="H2", to="H4", name="Tr2")
save(exSCE, file='exSCE.rda', compress='xz')

###############################################################################
# inst/exdata/exSCE.graphml
###############################################################################
#Export graph to graphml
write.ygraphml(exSCE, file="exSCE.graphml", color_by="phenoName", name="state",
               node_label="state")

# 1. Open file in yEd (https://www.yworks.com/products/yed)
# 2. Layout -> Tree -> Balloon
# 3. Save file

###############################################################################
# inst/exdata/bundle.rds
###############################################################################
# SingleCellExperiment object containing normalized multiplex RT-qPCR
# transcript expression profiles of 183 genes measured in 1,008 cells of the
# chicken utricle sensory epithelium at embryonic day 15. Genes were expected
# to be expressed during sensory hair cell bundle maturation and function.
# Experimental metadata was generated during tissue preparation (cell origin)
# and cell sorting (uptake of FM1-43 dye indicating cell maturity).
#
# Reference:
# Ellwanger DC, Scheibinger M, Dumont RA, Barr-Gillespie PG, and Heller S.
# Transcriptional dynamics of hair-bundle morphogenesis revealed with
# CellTrails. Cell Reports, 2018 Jun 5;23

###############################################################################
# inst/exdata/bundle.graphml
###############################################################################
bundle <- readRDS(system.file("exdata", "bundle.rds", package="CellTrails"))
se <- embedSamples(bundle)
d <- findSpectrum(se$eigenvalues, frac=100) #Similar to Figure 1E
latentSpace(bundle) <- se$components[, d]
states(bundle) <- findStates(bundle, min_size=0.01,
                             min_feat=5, max_pval=1e-04, min_fc=2)
bundle <- connectStates(bundle, l=10)
bundle <- selectTrajectory(bundle, component=1)
bundle <- fitTrajectory(bundle)
write.ygraphml(bundle, file="/Users/ellwanger/tmp/bundle.graphml",
               color_by="phenoName", name="state",
               node_label="state")

# 1. Open file in yEd (https://www.yworks.com/products/yed)
# 2. Layout -> Tree -> Balloon
# Optional: 3.1 Tools -> Geometric Transformations -> Mirror on X-axis
# Optional: 3.2 Tools -> Geometric Transformations -> Mirror on Y-axis
# 4. Save file

###############################################################################
# inst/exdata/th2.rds
###############################################################################
# SingleCellExperiment object containing normalized and filtered RNA-Seq
# transcript counts of 7,063 genes expressed in at least two of 81 cells
# as provided in the supplement of Buettner et al. 2015 (Table S5 and S7).
# The dataset measured gene expression of murine T helper 2 cell (Th2)
# differentiation; cell cycle effects were identified as a major confounder
# in this dataset. Metainformation on 116 manually annotated Th2 marker
# genes is provided in this object.
#
# References:
# Buettner F, Natarajan KN, Casale FP, Proserpio V, Scialdone A, Theis FJ,
# Teichmann SA, Marioni JC, and Stegle O.
# Computational analysis of cell-to-cell heterogeneity in single-cell
# RNA-sequencing data reveals hidden subpopulations of cells. Nature
# Biotechnology, 2015;33(2) 155-60
#
# Mahata B, Zhang X, Kolodziejczyk AA, Proserpio V, Haim-Vilmovsky L,
# Taylor AE, Hebenstreit D, Dingler FA, Moignard V, Goettgens B, Arlt W,
# McKenzie AN, Teichmann SA,
# Single-cell RNA sequencing reveals T helper cells synthesizing steroids
# de novo to contribute to immune homeostasis. Cell Reports, 2014;7(4) 1130-42
