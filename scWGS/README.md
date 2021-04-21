# Single cell Whole Genome Sequencing data analysis

This folder includes the Rmd files used for the analysis of the scWGS data.

The data were firstly processed with our [NextFlow Ginkgo pipeline](https://github.com/UCL-BLIC/nf-ginkgo). In addition
to the standard NGS pre-processing and alignment, this pipeline runs Ginkgo which we use to bin the reads and assess the
results (QC plots). The binned data are then re-normalised (cnv_data_processing_normal.Rmd & cnv_data_processing.Rmd),
segmented with [copynumber R package](https://bioconductor.org/packages/release/bioc/html/copynumber.html) (cnv_copynumber_\*.Rmd))
and trees are using the Minimum Arborescence Tree (i.e. rooted Minimum Spanning Tree) algorithm from the [optrees R
package](https://cran.r-project.org/web/packages/optrees/index.html).

1. cnv_data_processing_normal.Rmd: Pre-process the chord blood normal samples, defines bad bins and re-normalisation factors
2. cnv_data_processing.Rmd: Pre-process the sample data
3. resampling_\*.Rmd: Generate the trees
4. _other files:_ common functions 
