# DiffCircaPipeline R package

## Install from github

* In R console

```{R}


if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
#install dependencies
devtools::install_github("diffCircadian/diffCircadian") 
devtools::install_github("ricardo-bion/ggradar") 
devtools::install_github("Caleb-Huo/differentialR2") 
devtools::install_github("Caleb-Huo/AWFisher")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

#install DiffCircaPipeline
devtools::install_github("DiffCircaPipeline/Rpackage")
# devtools::install_github("DiffCircaPipeline/Rpackage@main") #run this if the above line does not work. 
```

## Tutorial 

The tutorial can be found at https://diffcircapipeline.github.io. 
