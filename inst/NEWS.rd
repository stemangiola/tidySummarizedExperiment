\name{NEWS}
\title{News for Package \pkg{tidySummarizedExperiment}}

\section{Changes in version 1.19.5}{
\itemize{
    \item Soft deprecated \code{bind_rows()} in favor of \code{append_samples()} from ttservice.
    \item Added \code{append_samples()} method for SummarizedExperiment objects.
    \item \code{bind_rows()} is not a generic method in dplyr and may cause conflicts.
    \item Users are encouraged to use \code{append_samples()} instead.
}}

\section{Changes in version 1.19.2, Bioconductor 3.22 Release}{
\itemize{
    \item Updated documentation to properly reflect S3 method structure.
    \item Simplified ggplot2 compatibility - S3 method continues to work with ggplot2 4.0.0.
    \item Users are now directed to library(tidyprint) for tidy visualization and https://github.com/tidyomics/tidyprint for more information.
}}

\section{Changes in version 1.4.0, Bioconductor 3.14 Release}{
\itemize{
    \item Improved join_*() functions.
    \item Changed special column names with a starting "." to avoid conflicts with pre-existing column names.
    \item Improved all method for large-scale datasets.
}}

\section{Changes in version 1.5.3, Bioconductor 3.15 Release}{
\itemize{
    \item Speed-up nest.
    \item Adaptation to Ranged-SummarizedExperiment.
}}

\section{Changes in version 1.7.3, Bioconductor 3.16 Release}{
\itemize{
    \item Fixed as_tibble edge case
    \item Fixed print for DelayedArray
    \item Improve performance for large-scale datasets
    \item Fixed filter is the result is a no-gene dataset, and improve performance of filtering
}}

