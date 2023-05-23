Single-cell RNA sequencing by combinatorial indexing often involves the
use of oligodTs and random hexamers combined for higher efficiency of
the reverse transcription step. If two different barcodes are used for
oligodTs and random hexamers, each cell in the raw UMI counts matrix
will have two possible barcodes, and the counts need to be merged before
further analysis. This package contains a function called merge_barcodes
that takes a matrix (`matrix`), a list of barcodes (`bc1_list`) and an
output directory as input. It merges the barcodes based on the barcode
list and creates a new sparse matrix with collapsed barcodes and summed
counts, which is written to the output directory as a .mtx and .RDS
file.

Install the package with

``` r
install.packages("scBact_1.0.tar.gz", repos = NULL)
```

An example of `bc1_list` is given in the file `bc1_collapsing.xlsx`, and
an example of raw count matrix is given in the file
`example_matrix.RDS`.

``` r
library(scBact)
bc1_list<-xlsx:::read.xlsx("bc1_collapsing.xlsx", sheetIndex = 1)
ex_matrix <- readRDS("example_matrix.RDS")
merge_barcodes(matrix = ex_matrix,
               bc1_list = bc1_list,
               output_dir = "./merge_barcodes_test",
               bc1_processed = NULL)
```

