# This function takes a sparse matrix, a barcode list, and output directory as input.
# It collapses the barcodes based on the barcode list and creates a new sparse matrix
# with collapsed barcodes, which is written to the output directory as a .mtx file.

#' Barcode collapsing tool for single cell RNA seq by split-pool barcoding
#'
#' @param matrix a dgT type matrix with barcodes as column names, in the format "ATGCCTAA_AATCCGTC_AAACGATA" and genes as row names
#' @param bc1_list a table containing a 'barcode1' columm and a 'barcode2' column, filled with the two possible round 1 barcodes in each well
#' @param output output directory. Current directory is default
#' @param bc1_processed optional. If the function was run before, give the directory of the 'colnames_BC1_processed.RDS' output file
#'
#' @return A new dgT type matrix with updated column names and summed counts
#' @export 

merge_barcodes <- function(matrix, bc1_list, output = getwd(), bc1_processed = NULL){
  
  # Get the column names of the sparse matrix
  colnames<-colnames(matrix)
  
  # If bc1_processed is not provided, process the barcode list
  if (rlang::is_empty(bc1_processed) == TRUE){
    
    # Split the column names by underscore and extract the last part
    colnames2<-strsplit(colnames, "_")
    colnames2<-sapply(colnames2, utils::tail, 1)
    
    # Read the barcode list
    my_excel <- bc1_list
    
    # Create an empty vector to store the processed barcode names
    to_return<-c()
    
    # Loop over the barcode names and look for matches in the barcode list
    idx_tmp <- 1
    for (string in colnames2){
      
      # Check if the barcode is in barcode2 column
      find1<-any(my_excel$barcode1==string)
      
      # Check if the barcode is in barcode1 column
      find2<-any(my_excel$barcode2==string)
      
      # If the barcode is in barcode2 column, add the corresponding barcode1 to the to_return vector
      if (find1==FALSE & find2==TRUE){
        to_return<-append(to_return, my_excel$barcode2[my_excel$barcode2==string])
        
        # If the barcode is in barcode1 column, add the corresponding barcode2 to the to_return vector
      } else if (find1==TRUE & find2==FALSE){
        to_return<-append(to_return, my_excel$barcode2[my_excel$barcode1==string])
        
        # If the barcode is in both columns, raise an error
      } else if (find1==TRUE & find2==TRUE){
        stop(paste0('Please check your barcode list, barcode of cell number ', idx_tmp, ' is the same in both columns'))
        
        # If the barcode is not in either column, raise an error
      } else if (find1==FALSE & find2==FALSE){
        stop(paste0('Please check your barcode list, barcode of cell number ', idx_tmp, ' was not found'))
        
        # If there is an issue with the barcode list, raise an error
      } else {
        stop('Please check your barcode list')
      }
      
      # Increment the counter and print status every 10,000 columns
      idx_tmp <- idx_tmp+1
      if(idx_tmp %% 10000 == 0){
        print(paste0("Creating new list of barcodes, processing matrix column ", idx_tmp))
      }
    }
    
    # Save the processed barcode names as an RDS file
    saveRDS(to_return, paste0(output, '/colnames_BC1_processed.RDS'))
    
    # If bc1_processed is provided, read the processed barcode names from the RDS file
  } else {
    to_return <- readRDS(bc1_processed)
  }
  
  # Extract the first part of the column names
  new_colnames_tmp<-sub("(.*?_.*?)_.*", "\\1", colnames)
  
  # Combine the first part of the column names with the processed barcode names
  colnames_final<-paste(new_colnames_tmp, to_return, sep = "_")
  
  # Update the column names of the sparse matrix
	colnames(matrix)<- colnames_final
	
	# Create a new sparse matrix with the same number of rows as the input matrix, but only one column
	final_sparse_Matrix <- MatrixExtra::emptySparse(nrow = nrow(matrix), ncol = 1L, format = "T", dtype = "d")
	
	# Give the one column a name that will be used to remove it later
	colnames(final_sparse_Matrix) <- "colToRm"
	
	# Initialize a counter to keep track of progress
	idx_tmp <- 1
	
	# Loop over each column in the input matrix that has a duplicate column name
	for(x in colnames(matrix)[duplicated(colnames(matrix))]) {
	  
	  # Extract the column from the input matrix that matches the current duplicate name
	  test_tmp <- matrix[,which(colnames(matrix) == x)]
	  
	  # Calculate the row sums for the extracted column
	  test_tmp2 <- Matrix::Matrix(rowSums(test_tmp, na.rm = TRUE))
	  
	  # Convert the row sums to a sparse matrix of type "dgTMatrix"
	  test_tmp2 <- methods::as(test_tmp2, "dgTMatrix")
	  
	  # Give the sparse matrix the same column name as the input matrix column it represents
	  colnames(test_tmp2) <- colnames(test_tmp)[1]
	  
	  # Give the sparse matrix the same row names as the input matrix
	  rownames(test_tmp2) <- rownames(test_tmp)
	  
	  # Append the sparse matrix as a new column in the final sparse matrix
	  final_sparse_Matrix <- cbind(final_sparse_Matrix, test_tmp2)
	  
	  # Increment the progress counter
	  idx_tmp <- idx_tmp+1
	  
	  # Print a progress message every 5000 columns
	  if(idx_tmp %% 5000 == 0){
	    print(paste0("Creating new count matrix, processing matrix column ", idx_tmp))
	  }
	}
	
	# Remove the "colToRm" column from the final sparse matrix
	final_sparse_Matrix2 <- final_sparse_Matrix[,-c(which(colnames(final_sparse_Matrix) == "colToRm"))]
	
	# Write the final sparse matrix to disk in Matrix Market format
	Matrix::writeMM(final_sparse_Matrix2, file = paste0(output, '/UniqueAndMult-Uniform_BC1_collapsed.mtx'))
	saveRDS(final_sparse_Matrix2, file = paste0(output, '/UniqueAndMult-Uniform_BC1_collapsed.RDS'))
	
	
	# Remove temporary objects
	gc()
}




