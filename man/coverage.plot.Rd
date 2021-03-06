% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coverage.plot.R
\name{coverage.plot}
\alias{coverage.plot}
\title{Coverage plot}
\usage{
coverage.plot(indir, outdir, names = c("HU 0 min", "", "", "", "", "",
  "", "", "", "HU 90 min", "NOC 0 min", "", "", "", "", "", "", "", "",
  "NOC 90 min"), num_files = c(1, 10, 11, 20), y_lim = 30,
  chromosome = "chrIII", start = 1, end = 30000)
}
\arguments{
\item{indir}{Directory containing bedgraph files with signals (e.g. after MACS2).}

\item{outdir}{Output directory.}
}
\description{
Coverage plot
}
