
#' @title Scores scatterplot matrix in the OODA style
#'
#' @description Scatterplot matrix of principal scores in the style of
#' Marron and Dryden (2021): off-diagonal panels show pairwise score
#' scatterplots, and diagonal panels show a kernel density estimate of each
#' score together with the jittered one-dimensional projections (a
#' "jitter plot"), so that clusters and outliers are visible in each
#' direction. Panel labels report the percentage of variance explained.
#'
#' @param scores matrix of scores (columns are components).
#' @param n_comp number of components to display.
#' @param col vector of colors for the observations.
#' @param var optional vector of variance percentages.
#' @param main overall title.
#' @param pch,cex point type and size.
#' @param bw bandwidth selector for the diagonal KDEs.
#' @param common_scale use the same axis range for all scores? This shows
#' the relative sizes of the modes of variation. Defaults to \code{FALSE}.
#' @return Nothing; produces a plot.
#' @references
#' Marron, J. S. and Dryden, I. L. (2021). \emph{Object Oriented Data
#' Analysis}. Chapman and Hall/CRC, Boca Raton.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 150, r = 3)
#' tp <- tangent_pca(X, d = rep(2, 3))
#' col <- hcl.colors(150, "Spectral")[rank(attr(X, "theta"))]
#' plot_scores_matrix(tp$scores, n_comp = 3, col = col, var = tp$var)
#' @export
plot_scores_matrix <- function(scores, n_comp = min(4, ncol(scores)),
                               col = 1, var = NULL, main = NULL, pch = 16,
                               cex = 0.5, bw = "SJ", common_scale = FALSE) {

  S <- scores[, seq_len(n_comp), drop = FALSE]
  col <- rep_len(col, nrow(S))
  op <- par(mfrow = c(n_comp, n_comp), mar = c(1.6, 1.6, 0.4, 0.4),
            oma = c(0.5, 0.5, if (is.null(main)) 0.5 else 2, 0.5),
            mgp = c(1.5, 0.3, 0), tcl = -0.2, cex.axis = 0.7)
  on.exit(par(op))
  lims <- apply(S, 2, function(s) range(s) + c(-1, 1) * 0.05 * diff(range(s)))
  if (common_scale) {

    lims[] <- rep(range(lims), ncol(S))

  }
  ord <- sample(nrow(S))  # random plotting order avoids overplot bias
  for (i in seq_len(n_comp)) {

    for (j in seq_len(n_comp)) {

      if (i == j) {

        kd <- tryCatch(density(S[, i], bw = bw),
                       error = function(e) density(S[, i]))
        yl <- c(-0.35, 1.05) * max(kd$y)
        plot(kd$x, kd$y, type = "l", xlim = lims[, i], ylim = yl, lwd = 1.5,
             xlab = "", ylab = "", yaxt = "n")
        jit <- -runif(nrow(S), 0.05, 0.3) * max(kd$y)
        points(S[ord, i], jit[ord], col = col[ord], pch = pch, cex = cex)
        lab <- colnames(S)[i]
        if (is.null(lab)) lab <- paste0("PC", i)
        if (!is.null(var)) lab <- sprintf("%s (%.1f%%)", lab, var[i])
        mtext(lab, side = 3, line = -1.1, adj = 0.04, cex = 0.6, font = 2)

      } else {

        plot(S[ord, j], S[ord, i], col = col[ord], pch = pch, cex = cex,
             xlim = lims[, j], ylim = lims[, i], xlab = "", ylab = "")

      }

    }

  }
  if (!is.null(main)) {

    mtext(main, outer = TRUE, line = 0.5, cex = 0.9, font = 2)

  }

}


#' @title Orthographic view of points on \eqn{\mathcal{S}^2}
#'
#' @description Draws a sphere outline with graticule and plots points
#' (front hemisphere solid, back hemisphere faded) after a rotation that
#' puts \code{view} at the center. Lightweight base-graphics alternative to
#' \pkg{rgl} for publication figures.
#'
#' @param x matrix of size \code{c(n, 3)} with points on \eqn{\mathcal{S}^2}.
#' @param col point colors.
#' @param view unit vector looking at the sphere.
#' @param lines optional list of \code{c(m, 3)} matrices drawn as curves.
#' @param lines_col,lwd colors and width of the curves.
#' @param add add to an existing plot?
#' @param cex point size.
#' @param main title.
#' @return Invisibly, the rotation matrix used.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 200, r = 1, height = 0.5)
#' plot_sphere(X, col = hcl.colors(200, "Spectral")[rank(attr(X, "theta"))],
#'             view = c(0, 0, 1))
#' @export
plot_sphere <- function(x, col = 1, view = c(1, 1, 1), lines = NULL,
                        lines_col = 1, lwd = 2, add = FALSE, cex = 0.5,
                        main = "") {

  view <- view / sqrt(sum(view^2))
  # Orthonormal frame (e1, e2 in the image plane; view is depth)
  e1 <- c(-view[2], view[1], 0)
  if (sum(e1^2) < 1e-12) e1 <- c(1, 0, 0)
  e1 <- e1 / sqrt(sum(e1^2))
  e2 <- c(view[2] * e1[3] - view[3] * e1[2], view[3] * e1[1] - view[1] * e1[3],
          view[1] * e1[2] - view[2] * e1[1])
  R <- cbind(e1, e2, view)
  if (!add) {

    plot(NA, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, axes = FALSE,
         xlab = "", ylab = "", main = main)
    tt <- seq(0, 2 * pi, l = 200)
    polygon(cos(tt), sin(tt), col = "grey97", border = "grey40")
    for (lat in seq(-60, 60, by = 30) * pi / 180) {

      g <- cbind(cos(lat) * cos(tt), cos(lat) * sin(tt), sin(lat)) %*% R
      g[g[, 3] < 0, 1:2] <- NA
      lines(g[, 1], g[, 2], col = "grey85", lwd = 0.5)

    }
    for (lon in seq(0, 150, by = 30) * pi / 180) {

      g <- cbind(cos(tt) * cos(lon), cos(tt) * sin(lon), sin(tt)) %*% R
      g[g[, 3] < 0, 1:2] <- NA
      lines(g[, 1], g[, 2], col = "grey85", lwd = 0.5)

    }

  }
  P <- x %*% R
  col <- rep_len(col, nrow(P))
  front <- P[, 3] >= 0
  points(P[!front, 1], P[!front, 2], col = adjustcolor(col[!front], 0.25),
         pch = 16, cex = cex)
  points(P[front, 1], P[front, 2], col = col[front], pch = 16, cex = cex)
  if (!is.null(lines)) {

    lines_col <- rep_len(lines_col, length(lines))
    for (l in seq_along(lines)) {

      L <- lines[[l]] %*% R
      Lf <- L
      Lf[L[, 3] < 0, 1:2] <- NA
      Lb <- L
      Lb[L[, 3] >= 0, 1:2] <- NA
      graphics::lines(Lb[, 1], Lb[, 2], col = adjustcolor(lines_col[l], 0.3),
                      lwd = lwd, lty = 2)
      graphics::lines(Lf[, 1], Lf[, 2], col = lines_col[l], lwd = lwd)

    }

  }
  invisible(R)

}


#' @title Plot method for SP-PCA fits
#'
#' @description Scores scatterplot matrix of an \code{"sppca"} fit in the
#' OODA style (see \code{\link{plot_scores_matrix}}), with the percentages
#' of variance explained on the scaled hypersphere in the diagonal labels.
#'
#' @param x an \code{"sppca"} object.
#' @param n_comp number of components to display.
#' @param ... further arguments passed to \code{\link{plot_scores_matrix}}
#' (e.g., \code{col}).
#' @return Nothing; produces a plot.
#' @examples
#' set.seed(1)
#' X <- r_circles_polysph(n = 150, r = 4, sigma = 0.05)
#' fit <- sppca(X, d = rep(2, 4), p = 4)
#' col <- hcl.colors(150, "Spectral")[rank(attr(X, "theta"))]
#' plot(fit, n_comp = 3, col = col)
#' @export
plot.sppca <- function(x, n_comp = min(4, ncol(x$scores)), ...) {

  plot_scores_matrix(x$scores, n_comp = n_comp, var = x$var_sphere, ...)

}
