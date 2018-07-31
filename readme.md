# MNuc package
R package containing tools for bioinformatics.  

## Installation
### Install MNuc package from gitHub
This is the fast way to make package ready for work, however, man pages will not be displayed.  

Print in RStudio command line:  
```{r}
library("devtools") 
install_github("suvarzz/MNuc")  
```

### Install MNuc package from downloaded archive.
Complete installation of MNuc package locally on your PC.  
1. Download MNuc, extract from archive, rename if necessury to "MNuc", and place it into a directory of choice (e.g. ~/Projects/MNuc).  
2. Install devtools and roxygen2 packages:  
```{r}
install.packages(c(devtools, roxygen2))
```

3. Print in RStudio command line:  
```{r}
library("devtools")  
setwd("~/Projects/MNuc/")  
document()  
setwd("~/Projects/")  
install("MNuc")  
```






