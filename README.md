## CellTrails: Inference of Temporal Gene Expression Dynamics of Branching Biological Processes from Single-cell Expression Data
_Daniel Christian Ellwanger_

Department of Otolaryngology - Head & Neck Surgery, and Institute for Stem Cell Biology and Regenerative Medicine, Stanford University School of Medicine, Stanford, CA 94305, USA

Package: CellTrails 0.99.8

  High-throughput single-cell technologies facilitate the generation of -omic readouts from thousands of cells captured at different cellular maturation stages during development, or other normal or pathological processes with unprecedented resolution. A single snapshot of an asynchronously developing specimen, for example, constitutes a time series in which individual cells represent distinct time points along a continuum. However, recoding of valuable cell-specific information, such as a cell's developmental age, its location in a tissue, or its functional phenotype is limited during sample preparation, and remains hidden in high dimensional cellular expression profiles. This formulates the computational challenge to infer the latent internal time axis of the biological process from the obtained expression matrix alone, while considering common parameters of single-cell measurements, such as noise, dropouts and redundancy. In other words, biological samples need to be placed by means of hidden information onto a non-linear trajectory, which might constitute of branching processes towards distinct functional cell types.

  This manual describes the practical use of the _CellTrails_ implementation, an unsupervised algorithm for the _de novo_ chronological ordering, visualization and analysis of single-cell expression data. _CellTrails_ makes use of a geometrically motivated concept of lower-dimensional manifold learning, which exhibits a multitude of virtues that counteract intrinsic noise of single cell data caused by drop-outs, technical variance, and redundancy of predictive variables. _CellTrails_ enables the reconstruction of branching trajectories and provides an intuitive graphical representation of expression patterns along all branches simultaneously. It allows the user to define and infer the expression dynamics of individual and multiple pathways towards distinct phenotypes.

  _CellTrails_ was developed with a 183-dimensional RT-qPCR gene expression panel of 1,008 cells collected from the developing chicken utricle, a balance organ. Key players in the utricle's function are cohorts of sensory hair cells that display mechanosensing organelles, called hair bundles, protruding from their apical surfaces. Bundle growth and maturation is dictated by an orchestration of distinct sequential and overlapping cellular processes. Our goal was to elucidate the temporal expression program of key hair bundle genes in subtypes of hair cells that occur with distinct spatial distribution. We showed that _CellTrails_ faithfully predicted expression patterns of hair cell maturation with unprecedented resolution.
  
  We confirmed that _CellTrails_ can be applied to analysis of single-cell RNA-Seq datasets. We are pleased that you consider using _CellTrails_ in your research. A detailed theoretical description of the algorithm and its application to biological uses has been published in:
  
  __Ellwanger DC, Scheibinger M, Dumont RA, Barr-Gillespie PG, and Heller S. "Transcriptional dynamics of hair-bundle morphogenesis revealed with CellTrails". _Journal_ Date;Issue.__

<!-- ---------------------------------- -->
## INSTALLATION
<!-- ---------------------------------- -->
*CellTrails* is an extension for _R_ (https://www.r-project.org), which is a free software environment for statistical computing and graphics. A simple yet efficient way to work with _R_ consists in writing source code with your favorite text editor and sending it to the _R_ console. It is suggested to use a development environment, such as _Rstudio_ (https://www.rstudio.com/), or a rich text editor with _R_ functionalities, such as _Emacs_ (https://www.gnu.org/software/emacs/), which greatly eases the work. 

The *CellTrails* package can be installed from this repository directly using the `devtools` package within an active _R_ session.

``` 
if(!require("devtools")) {
  install.packages("devtools")
} 
if(require("devtools")) {
  install_github("dcellwanger/CellTrails")
} else {
  stop("Could not load package 'devtools'.")
}
```
<!-- 
**If you are using macOS X**, please note that _CellTrails_ makes use of a library, which depends on _XQuartz_, a version of the X.Org X Window System (X11) that runs on OS X. If your system does not have _XQuartz_ installed yet, you need to download it from http://xquartz.org. Please note, that it also needs to be re-installed when upgrading your macOS to a new major version (https://cran.r-project.org/bin/macosx/). -->

**We also recommend** to download and install the graph visualization software _yEd_ (http://www.yworks.com/products/yed). It provides great capabilities to visualize and analyze a trajectory graph produced by *CellTrails*.

Before ready to use, the *CellTrails* library must be loaded into the _R_ environment:
```
library(CellTrails)
```

__Please, refer to the [vignette](https://dcellwanger.github.io/CellTrails/) for a detailed explanation and instruction on how to use CellTrails.__
