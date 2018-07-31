# How to install MNuc package from gitHub
1. Print in RStudion command line:  
```{r}
install_github("suvarzz/MNuc")  
```

# How to install MNuc package in RStudio
1. Download MNuc and place it in a directory (e.g ~/Projects/MNuc)  
2. Print in RStudio command line:  
```{r}
library("devtools")  
setwd("~/Projects/MNuc/")  
document()  
setwd("~/Projects/")  
install("MNuc")  
```






