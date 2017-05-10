# fTree
An R package that implements methods for growing functional regression trees

To install this package open R or Rstudio and type: devtools::install_github("ogru/fTree"). Make sure to have "devtools" package installed (install.packages("devtools")).

After installation completes view package vignette (NOT UPLOADED YET!) to get started: https://rawgit.com/ogru/fTree/master/vignette/fTree_get_started.html

If package installation fails make sure to update R and Rstudio to the latest version(s). 
This code relies on Rcpp and RcppArmadillo and it compiles cpp code upon installation. Focus your troubleshooting to that domain. 

If installation fails with an error: "‘to_string’ is not a member of ‘std’"

This means you need to tell g++ compiler to use C++11 with the following command (in R/Rstudio): Sys.setenv("PKG_CXXFLAGS"="-std=c++11")

Then repeat devtools::install_github("ogru/fTree")

Any other problems? Feel free to contact me. 


