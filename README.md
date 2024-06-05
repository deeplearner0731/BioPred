# BioPred

The R package BioPred offers a suite of tools for subgroup and biomarker analysis in precision medicine. Leveraging Extreme Gradient Boosting (XGBoost) along with propensity score weighting and A-learning methods, BioPred facilitates the optimization of individualized treatment rules (ITR) to streamline sub-group identification. BioPred also enables the identification of predictive biomarkers and obtaining their importance rankings. Moreover, the package provides graphical plots tailored for biomarker analysis. This tool enables clinical researchers seeking to enhance their understanding of biomarkers and patient popula-tion in drug development. 
## Installation

You can install the `BioPred` package from GitHub using the `devtools` package. If you don't have `devtools` installed, you can install it using:

```r
install.packages("devtools")
devtools::install_github("deeplearner0731/BioPred")
```
## Main Functions

- `BioPred()`: Main function for fitting predictive models.
- `predict()`: Function to make predictions based on the fitted model.
- `plot()`: Function to visualize the model and results.
- `fixcut_con()`: Function to fix cutoff values and perform statistical tests.

## Dependencies

The `BioPred` package depends on the following R packages:

- `ggplot2`
- `dplyr`
- `stats`

Make sure to install these packages before using `BioPred`.

## Author

This package is maintained by Zihuan Liu. For any queries or issues, please contact me at <zihuan.liu@abbvie.com>.

