#' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
    attached <- tidyverse_attach()
    
    # Print loading message about printing
    packageStartupMessage(
        "tidySummarizedExperiment says: Printing is now handled externally. ",
        "If you want to visualize the data in a tidy way, do library(tidyprint). ",
        "See https://github.com/tidyomics/tidyprint for more information."
    )
}