#' Get multiple distinct colors using viridis package.
#' @description
#' See details:
#' https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
#' @param types vector/array, which we want to assign each element with a distinct color
#' @param b numeric, [0,1] color begin value, 0.05 default
#' @param e numeric, [0,1] color end value, 0.2 default, no less than b
#' NOTE, both b and e are calculated from the right-side in the color scales.
#' @param pal literal, color palette to use, default is turbo (no "" needed)
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

#' Get continous color with linear interplotation
#' @description Use circlize::colorRamp2 to do the linear interplotation
#' Use wasanderson::was_palette to generate the explict colors needed,
#' which could be replaced in the future with other more common packages.
#' @param b numeric, region start, 0 as default
#' @param e numeric, region end, 1 as default
#' @param n integer, number of intervals for explicit color points withn [b, e]
#' @param pal characters, color palette, "Zissou1" as default
#' @param type characters, parameter for wesanderson pacakge, "continuous" as default
#' @return function, f(value within [b, e]) will return a color
#' @export
getContinuousColors <- function(b = 0,
                                e = 1,
                                n = 50,
                                pal = "Zissou1",
                                type = "continuous") {
  if (b > e)  {
    stop("b is larger than e.")
  }
  circlize::colorRamp2(breaks = seq(b, e, length = n),
                       wesanderson::wes_palette(name = pal,
                                                n = n,
                                                type = type))
}
