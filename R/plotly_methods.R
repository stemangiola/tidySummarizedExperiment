#' @name plot_ly
#' @rdname plot_ly
#' @inherit ttservice::plot_ly
#' @return `plotly`
#' 
#' @examples
#' data(se)
#' se |>
#'     plot_ly(x = ~counts)
#' 
#' @importFrom ttservice plot_ly
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' @export
plot_ly.tbl_df <- function(data=data.frame(), ..., type=NULL, name=NULL,
    color=NULL, colors=NULL, alpha=NULL,
    stroke=NULL, strokes=NULL, alpha_stroke=1,
    size=NULL, sizes=c(10, 100),
    span=NULL, spans=c(1, 20),
    symbol=NULL, symbols=NULL,
    linetype=NULL, linetypes=NULL,
    split=NULL, frame=NULL,
    width=NULL, height=NULL, source="A") {
    data |>

        # This is a trick to not loop the call
        drop_class("tbl_df") |>
        plotly::plot_ly(...,
            type=type, name=name,
            color=color, colors=colors, alpha=alpha,
            stroke=stroke, strokes=strokes, alpha_stroke=alpha_stroke,
            size=size, sizes=sizes,
            span=span, spans=spans,
            symbol=symbol, symbols=symbols,
            linetype=linetype, linetypes=linetypes,
            split=split, frame=frame,
            width=width, height=height, source=source
        )
}

#' @name plot_ly
#' @rdname plot_ly
#' @inherit ttservice::plot_ly
#' @return `plotly`
#' 
#' @examples
#' data(se)
#' se |>
#'     plot_ly(x = ~counts)
#' 
#' @importFrom ttservice plot_ly
#' @references
#' Hutchison, W.J., Keyes, T.J., The tidyomics Consortium. et al. The tidyomics ecosystem: enhancing omic data analyses. Nat Methods 21, 1166–1170 (2024). https://doi.org/10.1038/s41592-024-02299-2
#' @export
plot_ly.SummarizedExperiment <- function(data=data.frame(),
    ..., type=NULL, name=NULL,
    color=NULL, colors=NULL, alpha=NULL,
    stroke=NULL, strokes=NULL, alpha_stroke=1,
    size=NULL, sizes=c(10, 100),
    span=NULL, spans=c(1, 20),
    symbol=NULL, symbols=NULL,
    linetype=NULL, linetypes=NULL,
    split=NULL, frame=NULL,
    width=NULL, height=NULL, source="A") {
    data |>

        # This is a trick to not loop the call
        as_tibble() |>
        plotly::plot_ly(...,
            type=type, name=name,
            color=color, colors=colors, alpha=alpha,
            stroke=stroke, strokes=strokes, alpha_stroke=alpha_stroke,
            size=size, sizes=sizes,
            span=span, spans=spans,
            symbol=symbol, symbols=symbols,
            linetype=linetype, linetypes=linetypes,
            split=split, frame=frame,
            width=width, height=height, source=source
        )
}
