% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{gam_plot}
\alias{gam_plot}
\title{GAM Plot}
\usage{
gam_plot(
  yvar,
  censorvar = NULL,
  xvar,
  xvars.adj = NULL,
  sxvars.adj = NULL,
  type,
  data,
  k,
  pred.type = "iterms",
  link.scale = TRUE,
  title = "Trend Plot",
  ybreaks = NULL,
  xbreaks = NULL,
  rugcol.var = NULL,
  add.points = FALSE,
  prt.sum = TRUE,
  prt.chk = FALSE,
  outlier.rm = FALSE,
  newdat = NULL
)
}
\arguments{
\item{yvar}{Response variable name.}

\item{censorvar}{Censoring variable name for survival analysis (0-censored, 1-event).}

\item{xvar}{Biomarker name.}

\item{xvars.adj}{Potential confounding variables to adjust for using linear terms.}

\item{sxvars.adj}{Potential confounding variables to adjust for using curve terms.}

\item{type}{"c" for continuous, "s" for survival, and "b" for binary response.}

\item{data}{The dataset containing the variables.}

\item{k}{Upper limit on the degrees of freedom associated with an s smooth.}

\item{pred.type}{"iterms" for trend of xvar, "response" for Y at the original scale.}

\item{link.scale}{Whether to show the plot in the scale of the link function.}

\item{title}{Title of the plot.}

\item{ybreaks}{Breaks on the y-axis.}

\item{xbreaks}{Breaks on the x-axis.}

\item{rugcol.var}{Variable name defining the color of the rug and points.}

\item{add.points}{Whether to add data points to the plot.}

\item{prt.sum}{Whether to print summary or not.}

\item{prt.chk}{Whether to print model diagnosis.}

\item{outlier.rm}{Whether to remove outliers based on 1.5IQR.}

\item{newdat}{User-supplied customized data for prediction and plotting.}
}
\value{
A list containing p-table, s-table, GAM summary, GAM check, and the plot.
}
\description{
Generates a generalized additive model (GAM) plot for exploring the relationship between a response
variable and a biomarker.
}
