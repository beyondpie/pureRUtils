## ggplot2 wrapper from ArchR

#' A ggplot-based dot plot wrapper function
#'
#' This function is a wrapper around ggplot geom_point to allow for a more intuitive plotting of ArchR data.
#'
#' @param x A numeric vector containing the x-axis values for each point.
#' @param y A numeric vector containing the y-axis values for each point.
#' @param color A numeric/categorical vector used to determine the coloration for each point.
#' @param discrete A boolean value indicating whether the supplied data is discrete (`TRUE`) or continuous (`FALSE`).
#' @param discreteSet The name of a custom palette from `ArchRPalettes` to use for categorical/discrete color.
#' This argument is only used if `discrete` is set to `TRUE`.
#' @param continuousSet The name of a custom palette from `ArchRPalettes` to use for numeric color.
#' This argument is only used if `discrete` is set to `FALSE`.
#' @param labelMeans A boolean value indicating whether the mean of each categorical/discrete color should be labeled.
#' @param pal A custom palette used to override discreteSet/continuousSet for coloring vector.
#' @param defaultColor The default color for points that do not have another color applied (i.e. `NA` values).
#' @param highlightPoints A integer vector describing which points to hightlight. The remainder of points will be colored light gray.
#' @param colorDensity A boolean value indicating whether the density of points on the plot should be indicated by color.
#' If `TRUE`, continuousSet is used as the color palette.
#' @param size The numeric size of the points to be plotted.
#' @param xlim A numeric vector of two values indicating the lower and upper bounds of the x-axis on the plot.
#' @param ylim A numeric vector of two values indicating the lower and upper bounds of the y-axis on the plot.
#' @param extend A numeric value indicating the fraction to extend the x-axis and y-axis beyond the maximum and minimum
#' values if `xlim` and `ylim` are not provided. For example, 0.05 will extend the x-axis and y-axis by 5 percent on each end.
#' @param xlabel The label to plot for the x-axis.
#' @param ylabel The label to plot for the y-axis.
#' @param title The title of the plot.
#' @param randomize A boolean value indicating whether to randomize the order of the points when plotting.
#' @param seed A numeric seed number for use in randomization.
#' @param colorTitle A title to be added to the legend if `color` is supplied.
#' @param colorOrder A vector that allows you to control the order of palette colors associated with the values in `color`.
#' For example if you have `color` as `c("a","b","c")` and want to have the first color selected from the palette be used for
#' "c", the second color for "b", and the third color for "a", you would supply the `colorOrder` as `c("c", "b", "a")`.
#' @param colorLimits A numeric vector of two values indicating the lower and upper bounds of colors if numeric. Values
#' beyond these limits are thresholded.
#' @param alpha A number indicating the transparency to use for each point. See `ggplot2` for more details.
#' @param baseSize The base font size (in points) to use in the plot.
#' @param legendSize The size in inches to use for plotting the color legend.
#' @param ratioYX The aspect ratio of the x and y axes on the plot.
#' @param labelAsFactors A boolean indicating whether to label the `color` input as a numeric factor (`TRUE`) or with a character string (`FALSE`).
#' @param fgColor The foreground color of the plot.
#' @param bgColor The background color of the plot.
#' @param bgWidth The background relative width size of the halos in the labeling.
#' @param labelSize The numeric font size of labels.
#' @param addFit A string indicating the method to use for adding a fit/regression line to the plot (see `ggplot2::geom_smooth()` methods).
#' If set to `NULL`, no fit/regression line is added.
#' @param rastr A boolean value that indicates whether the plot should be rasterized using `ggrastr`. This does not rasterize
#' lines and labels, just the internal portions of the plot.
#' @param dpi The resolution in dots per inch to use for the plot.
#' @importFrom dplyr %>%
#' @import ggplot2
#' @export
ggPoint <- function(x = NULL,
                    y = NULL,
                    color = NULL,
                    discrete = TRUE,
                    discreteSet = "stallion",
                    continuousSet = "solarExtra",
                    labelMeans = TRUE,
                    pal = NULL,
                    defaultColor = "lightGrey",
                    highlightPoints = NULL,
                    colorDensity = FALSE,
                    size = 1,
                    xlim = NULL,
                    ylim = NULL,
                    extend = 0.05,
                    xlabel = "x",
                    ylabel = "y",
                    title = "",
                    randomize = FALSE,
                    seed = 1,
                    colorTitle = NULL,
                    colorOrder = NULL,
                    colorLimits = NULL,
                    alpha = 1,
                    baseSize = 10,
                    legendSize = 3,
                    ratioYX = 1,
                    labelAsFactors = TRUE,
                    fgColor = "black",
                    bgColor = "white",
                    bgWidth = 1,
                    labelSize = 3,
                    addFit = NULL,
                    rastr = FALSE,
                    dpi = 300,
                    legendDirection = "vertical",
                    legendPosition = "right",
                    legendFontSize = 5,
                    ...) {
  stopifnot(length(y) == length(x))
  if (length(x) < 5) {
    stop("x must be at least length 5 to plot!")
  }

  if (randomize) {
    set.seed(seed)
    idx <- sample(seq_along(x), length(x))
  } else {
    idx <- seq_along(x)
  }

  df <- data.frame(x = x, y = y)
  include <- which(is.finite(x) & is.finite(y))

  if (length(include) != length(x)) {
    message("Some values are not finite! Excluding these points!")
    df <- df[include, ]
    x <- x[include]
    y <- y[include]
    if (!is.null(color)) {
      color <- color[include]
    }
  }

  if (is.null(xlim)) {
    xlim <- range(df$x) %>% extendrange(f = extend)
  }

  if (is.null(ylim)) {
    ylim <- range(df$y) %>% extendrange(f = extend)
  }

  ratioXY <- ratioYX * diff(xlim) / diff(ylim)

  if (is.null(color) & !colorDensity) {
    p <- ggplot(df[idx, ], aes(x = x, y = y)) +
      coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = F) +
      xlab(xlabel) +
      ylab(ylabel) +
      ggtitle(title) +
      theme_ArchR(baseSize = baseSize)

    if (rastr) {
      p <- p + .geom_point_rast2(
        size = size, raster.dpi = dpi, alpha = alpha, color = defaultColor
      )
    } else {
      p <- p + geom_point(size = size, alpha = alpha, color = defaultColor)
    }
  } else {
    if (colorDensity) {
      discrete <- FALSE
      df <- .getDensity(x, y, n = 100, sample = NULL) # change
      df <- df[order(df$density), , drop = FALSE]
      df$color <- df$density

      if (is.null(colorTitle)) {
        colorTitle <- "density"
      }
    } else if (discrete) {
      if (!is.null(highlightPoints)) {
        if (length(highlightPoints) < length(color)) {
          color[-highlightPoints] <- "Non.Highlighted"
          idx <- c(idx[-highlightPoints], idx[highlightPoints])
        }
      }
      color <- paste0(color)

      if (!is.null(colorOrder)) {
        if (!all(color %in% colorOrder)) {
          stop("Not all colors are in colorOrder!")
        }
      } else {
        colorOrder <- gtools::mixedsort(unique(color))
      }

      if (is.null(colorTitle)) {
        colorTitle <- "color"
      }

      stopifnot(length(color) == nrow(df))
      df$color <- factor(color, levels = colorOrder)

      if (labelAsFactors) {
        df$color <- factor(
          x = paste0(paste0(match(paste0(df$color), paste0(levels(df$color)))), "-", paste0(df$color)),
          levels = paste0(seq_along(levels(df$color)), "-", levels(df$color))
        )
        if (!is.null(pal)) {
          # print(pal)
          # print(paste0(levels(df$color))[match(names(pal), colorOrder)])
          names(pal) <- paste0(levels(df$color))[match(names(pal), colorOrder)]
        }
        colorOrder <- paste0(levels(df$color))
      }
    } else {
      stopifnot(length(color) == nrow(df))
      if (!is.null(highlightPoints)) {
        if (length(highlightPoints) < length(color)) {
          color[-highlightPoints] <- NA
          idx <- c(idx[-highlightPoints], idx[highlightPoints])
        }
      }
      if (!is.null(colorLimits)) {
        color[color < min(colorLimits)] <- min(colorLimits)
        color[color > max(colorLimits)] <- max(colorLimits)
      }
      df$color <- color
    }

    p <- ggplot(df[idx, ], aes(x = x, y = y, color = color)) +
      coord_equal(ratio = ratioXY, xlim = xlim, ylim = ylim, expand = FALSE) +
      xlab(xlabel) +
      ylab(ylabel) +
      ggtitle(title) +
      theme_ArchR(baseSize = baseSize) +
      theme(legend.direction = legendDirection,
            legend.box.background = element_rect(color = NA)) +
      labs(color = colorTitle)

    if (rastr) {
      p <- p + .geom_point_rast2(
        size = size, raster.dpi = dpi, alpha = alpha,
        raster.width = min(par("fin")),
        raster.height = (ratioYX * min(par("fin")))
      )
    } else {
      p <- p + geom_point(size = size, alpha = alpha)
    }

    if (discrete) {
      if (!is.null(pal)) {
        p <- p + scale_color_manual(values = pal)
      } else {
        pal <- paletteDiscrete(set = discreteSet, values = colorOrder)
        if (!is.null(highlightPoints)) {
          pal[grep("Non.Highlighted", names(pal))] <- "lightgrey"
        }
        # print(pal)
        p <- p + scale_color_manual(values = pal) +
          guides(color = guide_legend(override.aes = list(size = legendSize, shape = 15)))
      }

      if (labelMeans) {
        dfMean <- split(df, df$color) %>%
          lapply(., function(x) {
            data.frame(x = median(x[, 1]), y = median(x[, 2]), color = x[1, 3])
          }) %>%
          Reduce("rbind", .)

        if (labelAsFactors) {
          dfMean$label <- stringr::str_split(paste0(seq_len(nrow(dfMean))), pattern = "\\-", simplify = TRUE)[, 1]
        } else {
          dfMean$label <- dfMean$color
        }
        dfMean$text <- stringr::str_split(dfMean$color, pattern = "-", simplify = TRUE)[, 1]

        # make halo layers, similar to https://github.com/GuangchuangYu/shadowtext/blob/master/R/shadowtext-grob.R#L43
        theta <- seq(pi / 8, 2 * pi, length.out = 16)
        xo <- bgWidth * diff(range(df$x)) / 300
        yo <- bgWidth * diff(range(df$y)) / 300
        for (i in theta) {
          p <- p +
            geom_text(
              data = dfMean,
              aes_q(
                x = bquote(x + .(cos(i) * xo)),
                y = bquote(y + .(sin(i) * yo)),
                label = ~text
              ),
              size = labelSize,
              color = bgColor
            )
        }

        if (is.null(fgColor)) {
          p <- p + geom_text(
            data = dfMean, aes(x = x, y = y, color = color, label = label),
            size = labelSize, show.legend = FALSE
          )
        } else {
          p <- p + geom_text(
            data = dfMean, aes(x = x, y = y, label = label), color = fgColor,
            size = labelSize, show.legend = FALSE
          )
        }
      }
    } else {
      if (!is.null(pal)) {
        if (!is.null(colorLimits)) {
          p <- p + scale_colour_gradientn(colors = pal, limits = colorLimits, na.value = "lightgrey")
        } else {
          p <- p + scale_colour_gradientn(colors = pal, na.value = "lightgrey")
        }
      } else {
        if (!is.null(colorLimits)) {
          p <- p + scale_colour_gradientn(
            colors = paletteContinuous(set = continuousSet),
            limits = colorLimits, na.value = "lightgrey"
          )
        } else {
          p <- p + scale_colour_gradientn(
            colors = paletteContinuous(set = continuousSet),
            na.value = "lightgrey"
          )
        }
      }
    }
  }
  if (!is.null(addFit)) {
    p <- p + geom_smooth(data = df, aes(color = NULL), method = addFit, color = "black") +
      ggtitle(paste0(
        title, "\nPearson = ", round(cor(df$x, df$y), 3),
        "\nSpearman = ", round(cor(df$x, df$y, method = "spearman"), 3)
      ))
  }
  p <- p + theme(legend.position = legendPosition,
                 ## legend.key = element_rect(size = 3),
                 legend.text = element_text(size = legendFontSize)
                 )
  if (!is.null(ratioYX)) {
    attr(p, "ratioYX") <- ratioYX
  }

  return(p)
}


.checkCairo <- function() {
  tryCatch({
    tmp <- dev.cur()
    Cairo::Cairo(type = 'raster')
    dev.off()
    dev.set(tmp)
    TRUE
  }, error = function(e) {
    FALSE
  })
}

## Adapted from
## https://github.com/tidyverse/ggplot2/blob/660aad2db2b3495ae0d8040915a40d247133ffc0/R/geom-point.r
## from https://github.com/VPetukhov/ggrastr/blob/master/R/geom-point-rast.R
## This funciton now handles issues with Cairo installation that can lead to plot errors
.geom_point_rast2 <- function(mapping = NULL,
                              data = NULL,
                              stat = "identity",
                              position = "identity",
                              ...,
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE,
                              raster.width = min(par("fin")),
                              raster.height = min(par("fin")),
                              raster.dpi = 300) {
  GeomPointRast <- tryCatch(
    {
      if (!.checkCairo()) {
        stop()
      }

      # Try to create a geom rast for points if not then just use normal geom_point
      ggplot2::ggproto(
        "GeomPointRast",
        ggplot2::GeomPoint,
        required_aes = c("x", "y"),
        non_missing_aes = c("size", "shape", "colour"),
        default_aes = aes(
          shape = 19, colour = "black", size = 1.5, fill = NA,
          alpha = NA, stroke = 0.5
        ),
        draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                              raster.width = min(par("fin")), raster.height = min(par("fin")), raster.dpi = 300) {

          # From ggrastr
          prevDevID <- dev.cur()

          p <- ggplot2::GeomPoint$draw_panel(data, panel_params, coord)

          devID <- Cairo::Cairo(
            type = "raster",
            width = raster.width * raster.dpi,
            height = raster.height * raster.dpi,
            dpi = raster.dpi,
            units = "px",
            bg = "transparent"
          )[1]

          grid::pushViewport(grid::viewport(width = 1, height = 1))

          grid::grid.points(
            x = p$x,
            y = p$y,
            pch = p$pch,
            size = p$size,
            name = p$name,
            gp = p$gp,
            vp = p$vp,
            draw = TRUE
          )

          grid::popViewport()
          gridCapture <- grid::grid.cap()

          dev.off(devID)

          dev.set(prevDevID)

          grid::rasterGrob(
            gridCapture,
            x = 0,
            y = 0,
            width = 1,
            height = 1,
            default.units = "native",
            just = c("left", "bottom")
          )
        }
      )
    },
    error = function(e) {
      if (.checkCairo()) {
        message("WARNING: Error found with trying to rasterize geom. Continuing without rasterization.")
      }
      ## else {
      ##   message("WARNING: Error found with Cairo installation. Continuing without rasterization.")
      ## }

      # Default geom_point
      ggplot2::ggproto(
        "GeomPoint",
        ggplot2::GeomPoint,
        required_aes = c("x", "y"),
        non_missing_aes = c("size", "shape", "colour"),
        default_aes = aes(
          shape = 19, colour = "black", size = 1.5, fill = NA,
          alpha = NA, stroke = 0.5
        ),
        draw_panel = function(data, panel_params, coord, na.rm = FALSE,
                              raster.width = min(par("fin")), raster.height = min(par("fin")), raster.dpi = 300) {
          if (is.character(data$shape)) {
            data$shape <- ggplot2:::translate_shape_string(data$shape) # Hidden ggplot2
          }

          coords <- coord$transform(data, panel_params)

          pGrob <- grid::pointsGrob(
            x = coords$x,
            y = coords$y,
            pch = coords$shape,
            gp = grid::gpar(
              col = scales::alpha(coords$colour, coords$alpha),
              fill = scales::alpha(coords$fill, coords$alpha),
              # Stroke is added around the outside of the point
              fontsize = coords$size * .pt + coords$stroke * .stroke / 2,
              lwd = coords$stroke * .stroke / 2
            )
          )

          pGrob
        },
        draw_key = ggplot2::draw_key_point
      )
    }
  )

  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomPointRast,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      raster.width = raster.width,
      raster.height = raster.height,
      raster.dpi = raster.dpi,
      ...
    )
  )
}


#' ggplot2 default theme for ArchR
#'
#' This function returns a ggplot2 theme that is black borded with black font.
#'
#' @param color The color to be used for text, lines, ticks, etc for the plot.
#' @param textFamily The font default family to be used for the plot.
#' @param baseSize The base font size (in points) to use in the plot.
#' @param baseLineSize The base line width (in points) to be used throughout the plot.
#' @param baseRectSize The base line width (in points) to use for rectangular boxes throughout the plot.
#' @param plotMarginCm The width in centimeters of the whitespace margin around the plot.
#' @param legendPosition The location to put the legend. Valid options are "bottom", "top", "left", and "right.
#' @param legendTextSize The base text size (in points) for the legend text.
#' @param axisTickCm The length in centimeters to be used for the axis ticks.
#' @param xText90 A boolean value indicating whether the x-axis text should be rotated 90 degrees counterclockwise.
#' @param yText90 A boolean value indicating whether the y-axis text should be rotated 90 degrees counterclockwise.
#' @import ggplot2
#' @export
theme_ArchR <- function(color = "black",
                        textFamily = "sans",
                        baseSize = 10,
                        baseLineSize = 0.5,
                        baseRectSize = 0.5,
                        plotMarginCm = 1,
                        legendPosition = "bottom",
                        legendTextSize = 5,
                        axisTickCm = 0.1,
                        xText90 = FALSE,
                        yText90 = FALSE) {
  theme <- theme_bw() + theme(
    text = element_text(family = textFamily),
    axis.text = element_text(color = color, size = baseSize),
    axis.title = element_text(color = color, size = baseSize),
    title = element_text(color = color, size = baseSize),
    plot.margin = unit(c(plotMarginCm, plotMarginCm, plotMarginCm, plotMarginCm), "cm"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = color, size = (4 / 3) * baseRectSize * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
    axis.ticks.length = unit(axisTickCm, "cm"),
    axis.ticks = element_line(color = color, size = baseLineSize * (4 / 3) * as.numeric(grid::convertX(grid::unit(1, "points"), "mm"))),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.text = element_text(color = color, size = legendTextSize),
    legend.box.background = element_rect(color = NA),
    # legend.box.background = element_rect(fill = "transparent"),
    legend.position = legendPosition,
    strip.text = element_text(size = baseSize, color = "black") # ,
    # plot.background = element_rect(fill = "transparent", color = NA)
  )

  if (xText90) {
    theme <- theme %+replace% theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }

  if (yText90) {
    theme <- theme %+replace% theme(axis.text.y = element_text(angle = 90, vjust = 1))
  }

  return(theme)
}


.getDensity <- function(x = NULL, y = NULL, n = 100, sample = NULL, densityMax = 0.95) {
  # modified from http://slowkow.com/notes/ggplot2-color-by-density/
  df <- data.frame(x = x, y = y)
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  df$density <- dens$z[ii]
  df$density[df$density > quantile(unique(df$density), densityMax)] <- quantile(unique(df$density), densityMax) # make sure the higher end doesnt bias colors
  if (!is.null(sample)) {
    df <- df[sample(nrow(df), min(sample, nrow(df))), ]
  }
  return(df)
}

#' @import ggplot2
#' @import grid
#' @import gridExtra
#' @export
.fixPlotSize <- function(p = NULL,
                         plotWidth = unit(6, "in"),
                         plotHeight = unit(6, "in"),
                         margin = 0.25,
                         height = 1,
                         it = 0.05,
                         newPage = FALSE) {
  if (!inherits(plotWidth, "unit")) {
    plotWidth <- unit(plotWidth, "in")
  }

  if (!inherits(plotHeight, "unit")) {
    plotHeight <- unit(plotHeight, "in")
  }

  # adapted from https://github.com/jwdink/egg/blob/master/R/set_panel_size.r
  g <- ggplotGrob(p)

  legend <- grep("guide-box", g$layout$name)
  if (length(legend) != 0) {
    gl <- g$grobs[[legend]]
    g <- ggplotGrob(p + theme(legend.position = "none"))
  } else {
    gl <- NULL
    g <- ggplotGrob(p)
  }

  panels <- grep("panel", g$layout$name)
  panel_index_w <- unique(g$layout$l[panels])
  panel_index_h <- unique(g$layout$t[panels])

  nw <- length(panel_index_w)
  nh <- length(panel_index_h)

  pw <- convertWidth(plotWidth, unitTo = "in", valueOnly = TRUE)
  ph <- convertWidth(plotHeight, unitTo = "in", valueOnly = TRUE)

  pw <- pw * 0.95
  ph <- ph * 0.95

  x <- 0
  width <- 1
  sm <- FALSE

  while (!sm) {
    x <- x + it

    w <- unit(x * width, "in")
    h <- unit(x * height / width, "in")
    m <- unit(x * margin / width, "in")

    g$widths[panel_index_w] <- rep(w, nw)
    g$heights[panel_index_h] <- rep(h, nh)

    sw <- convertWidth(
      x = sum(g$widths) + m,
      unitTo = "in",
      valueOnly = TRUE
    )

    sh <- convertHeight(
      x = sum(g$heights) + m,
      unitTo = "in",
      valueOnly = TRUE
    )

    sm <- sw > pw | sh > ph
  }

  if (length(legend) != 0) {
    sgh <- convertHeight(
      x = sum(g$heights),
      unitTo = "in",
      valueOnly = TRUE
    )
    sgw <- convertWidth(
      x = sum(g$widths),
      unitTo = "in",
      valueOnly = TRUE
    )
    slh <- convertHeight(
      x = sum(gl$heights),
      unitTo = "in",
      valueOnly = TRUE
    )
    slw <- convertWidth(
      x = sum(gl$widths),
      unitTo = "in",
      valueOnly = TRUE
    )
    size <- 6
    wh <- 0.1
    it <- 0

    while (slh > 0.2 * ph | slw > pw) {
      it <- it + 1

      if (it > 3) {
        break
      }
      size <- size * 0.8
      wh <- wh * 0.8
      gl <- ggplotGrob(
        p + theme(
          legend.key.width = unit(wh, "cm"),
          legend.key.height = unit(wh, "cm"),
          legend.spacing.x = unit(0, "cm"),
          legend.spacing.y = unit(0, "cm"),
          legend.text = element_text(size = max(size, 2))
        ) + guides(fill = guide_legend(ncol = 4), color = guide_legend(ncol = 4))
      )$grobs[[legend]]

      slh <- convertHeight(
        x = sum(gl$heights),
        unitTo = "in",
        valueOnly = TRUE
      )
      slw <- convertWidth(
        x = sum(gl$widths),
        unitTo = "in",
        valueOnly = TRUE
      )
    }
    p <- grid.arrange(g, gl,
      ncol = 1, nrow = 2,
      heights = unit.c(unit(sgh, "in"), unit(min(slh, 0.2 * pw), "in")),
      newpage = newPage
    )
  } else {
    p <- grid.arrange(g, newpage = newPage)
  }
  invisible(p)
}
