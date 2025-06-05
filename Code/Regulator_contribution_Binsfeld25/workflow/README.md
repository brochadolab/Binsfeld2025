Compiled at 2025-05-29 19:58:09 UTC

``` r
here::i_am("README.Rmd", uuid = "ec7bb4fc-a51d-4a22-985b-20747e1ab09e")

# function to get path to previous data: path_source("99-publish", "sample.csv")
path_source <- projthis::proj_path_source("README")
```

In this workflow, we determine regulator contributions to observed
changes in gene expression through a hierarchical interaction model
(hierNet).

The workflow is divided into four steps:

1.  `00-import`: Imports the CPI scores from the parent directory.
2.  `01-data_exploration`: Filters experiments, and performs quantile
    normalization.
3.  `02-expression_modeling`: Water thresholding and hierarchical
    interaction modeling.  
4.  `03-cross_validation`: Performs 10-fold cross validation.

You can execute all steps by running the following command in the parent
directory.

``` r
projthis::proj_workflow_render("workflow")
```

The intermediate data files will be found in `data`.
