This package contains a function called merge_barcodes that takes a
matrix (`matrix`) and a list of barcodes (`bc1_list`). `matrix` is the
raw UMI/cell count matrix (with barcodes as column names in the format
BC3_BC2_BC1 and features as row names). An example of `bc1_list` is
given in the file `bc1_collapsing.xlsx`, and can be read by
`bc1_list<-xlsx:::read.xlsx("bc1_collapsing.xlsx", sheetIndex = 1)`. It
also takes two optional arguments: an output directory (`output`,
current directory is default) and a pre-processed list of barcodes (an
output of the function, can be provided when you don’t run the function
for the first time, `bc1_processed`).

The function first checks if `bc1_processed` is NULL. If it is, the
function proceeds to process the barcodes in the matrix’s column names
by extracting the last part of the column name (the barcode 1). Reading
in the barcode list from an Excel file using the `read.xlsx` function
from the `xlsx` package, the function checks if the barcode is present
in column `barcode1` or `barcode2` (representing the two possibilities
of barcode 1 in each well) of the barcode list and appends the
corresponding barcode to `to_return`. If the barcode is not found or
found in both columns, an error message is thrown.

If `bc1_processed` is not NULL, `to_return` is read from the RDS file
specified by `bc1_processed`.

The function then modifies the matrix’s column names by replacing the
last part of the name with the corresponding barcode from `to_return`.
Next, the function creates a new sparse matrix (`final_sparse_Matrix`)
and iterates over duplicated column names in `matrix`. For each
duplicated column name, the function extracts the corresponding column
from the matrix and creates a new sparse matrix (`test_tmp2`) that sums
the row values and sets the column name to the original column name’s
first instance. The resulting sparse matrix is then appended to
`final_sparse_Matrix`.

It will output 2 files: a file containing the list of new barcodes named
`colnames_BC1_processed.RDS` and a new count matrix with updated
colnames (half of the initial number of columns) and summed counts named
UniqueAndMult-Uniform_BC1_collapsed.mtx

The function finally removes the temporary objects and runs the garbage
collector.

Install the package with

``` r
install.packages("scBact.tar.gz")
```
