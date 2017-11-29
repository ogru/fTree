# fTree

An R package that implements the methods for growing regression trees with functional output data.
   
## Package Installation
To install the package open `R` or `Rstudio` and type: `devtools::install_github("ogru/fTree")`. Make sure to have the `devtools` package installed (i.e. `install.packages("devtools")`).
  
## Tutorial/Vignette
After installing view the [vignette](https://rawgit.com/ogru/fTree/master/vignette/first_steps.html) to get started.
  
## Troubleshooting
 
If the package installation fails, make sure to update your `R` and `Rstudio` to the latest version(s). 
  
This code relies on the `Rcpp` and `RcppArmadillo` packages and it compiles `Cpp` code upon installation. In case of problems, focus your troubleshooting to that domain. 
  
If the installation fails with the following error: "‘to_string’ is not a member of ‘std’", this means that you need to set up your `g++` compiler to use `C++11`. You can accomplish this with the following command (in R/Rstudio): `Sys.setenv("PKG_CXXFLAGS"="-std=c++11")`. Then repeat `devtools::install_github("ogru/fTree")`.
  
## Contact
 
Any problems? Raise an issue or send me an email ognjengr@gmail.com. 

## License
The package is distributed under the terms of the GPL-2 license (see LICENSE.txt).

## Citation

If you find this package useful in your work please consider citing it in your publication(s) with the following Bibtex entry:
```
@Manual{fTree_2017,
  title  = {fTree, an \proglang{R} package},
  author = {Ognjen Grujic},
  year   = {2017},
  note   = {\proglang{R}~package version~1.01},
  url    = {https://github.com/ogru/fTree},
}
```
