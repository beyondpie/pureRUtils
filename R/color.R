#' Get multiple distinct colors using viridis package.
#' @description
#' See details:
#' https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
#' @param types vector/array, which we want to assign each element with a distinct color
#' @param b numeric, [0,1] color begin value, 0.05 default
#' @param e numeric, [0,1] color end value, 0.2 default, no less than b
#' NOTE, both b and e are calculated from the right-side in the color scales.
#' @param pal literal, color pannel to use, default is turbo (no "" needed)
#' @param showColor bool, after generating the color
#' if need to show colors with scales::show_col, default FALSE
#' @return vector of characters with types as names. Each element represents one color.
#' @export
getMultipleColors <- function(types, b = 0.05, e = 0.2,
                              pal = turbo,
                              showColor = FALSE) {
  if (length(unique(types)) < length(types)) {
    warning("types has repeat elements. Only unique elements are considered")
  }
  types <- unique(types)
  if(b > e) {
    stop("b is larger than e.")
  }
  eval(substitute({
    colors <- viridis::pal(n = length(types),
    begin = b,
    end = e)
    if(length(unique(colors)) < length(types)) {
      warning("Get ", length(unique(colors)), " colors, which is less than ",
              length(types), " types.")
    }
    names(colors) <- types
    if(showColor) {
      scales::show_col(colors)
    }
    colors
  }))
}

#' @export
getContinuousColors <- function(b = 0, e = 1, n = 50, pal = "Zissou1", type = "continuous") {
  circlize::colorRamp2(breaks = seq(b, e, length = n),
                       wesanderson::wes_palette(name = pal,
                                                n = n,
                                                type = type))
}
