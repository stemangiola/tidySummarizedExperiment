\name{NEWS}
\title{News for Package \pkg{tidySummarizedExperiment}}

\section{Changes in version 1.16.0, Bioconductor 3.17 Release}{
\itemize{
    \item Prepared for ggplot2 4.0.0 compatibility with S7 methods.
    \item Added conditional S7 method for ggplot while maintaining S3 method for current ggplot2 versions.
    \item Added S7 as suggested dependency for future ggplot2 4.0.0 compatibility.
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

