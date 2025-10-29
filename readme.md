# `annoLinker`: Annotating Genomic Regions Through Chromatin Interaction Links

## Overview

`annoLinker` provides a framework for annotating genomic regions by integrating chromatin interaction data.
Traditional annotation approaches often link peaks to the nearest genes based solely on linear genomic distance, overlooking the three-dimensional genome architecture revealed by Hi-C, ChIA-PET, PLAC-seq, HiCAR, and related assays.

`annoLinker` bridges this gap by constructing chromatin interaction networks using the igraph package and propagating gene annotations through connected genomic elements.
This approach enables the functional annotation of distal regulatory regions such as enhancers and silencers that interact with target genes through chromatin looping.

## Key Features

* Annotates genomic regions using DNA interaction data rather than linear proximity.

* Compatible with GRanges, GInteractions, and Pairs objects.

* Supports annotation of promoter, gene body, or downstream regions.

* Includes visualization tools.

## Installation

### Install from GitHub
```{r}
if (!requireNamespace("remotes")) install.packages("remotes")

remotes::install_github("jianhong/annoLinker")
```

## Contributing

Contributions, suggestions, and bug reports are welcome!
Please open an issue or submit a pull request on [GitHub](https://github.com/jianhong/annoLinker/issues).
