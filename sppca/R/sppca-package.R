
#' @title \pkg{sppca}: Scaled Polyspherical Principal Component Analysis
#'
#' @description Dimension reduction for data on products of spheres
#' (polyspheres) \eqn{\mathcal{S}^{d_1} \times \cdots \times
#' \mathcal{S}^{d_r}}, such as the spoke directions of skeletal
#' representations (s-reps) of anatomical shapes.
#'
#' @details Scaled polyspherical PCA (\code{\link{sppca}}) proceeds in three
#' steps:
#' \enumerate{
#'   \item \emph{Scaling}: \code{\link{smds_polysph}} maps the sample to a
#'   hypersphere \eqn{\mathcal{S}^p_\rho} of estimated radius by a landmark
#'   scaled spherical multidimensional scaling with profiled radius, analytic
#'   gradients, spherical classical-MDS initialization, and spherical Gower
#'   interpolation refined by Levenberg--Marquardt.
#'   \item \emph{Nesting}: \code{\link{pns_fit}} computes principal nested
#'   spheres of the scaled sample.
#'   \item \emph{Inversion}: \code{\link{nw_sph_polysph}} and
#'   \code{\link{ll_scores_polysph}} map modes of variation and
#'   reconstructions back to the polysphere, with bandwidths selected by
#'   closed-form leave-one-out cross-validation.
#' }
#' Competing methods (\code{\link{tangent_pca}}, \code{\link{extrinsic_pca}},
#' \code{\link{cpns}}), a method-agnostic cross-validated explained variance
#' (\code{\link{residual_scores}}), samplers (\code{\link{r_circles_polysph}},
#' \code{\link{r_clusters_polysph}}), and plots in the style of object
#' oriented data analysis (\code{\link{plot_scores_matrix}},
#' \code{\link{plot_sphere}}) are also provided. The package interoperates
#' with \pkg{polykde} (e.g., the \code{\link[polykde]{hippocampus}} s-reps)
#' and \pkg{shapes}.
#' @author Michael Nisenzon, Eduardo García-Portugués, and J. S. Marron.
#' @references
#' Jung, S., Dryden, I. L. and Marron, J. S. (2012). Analysis of principal
#' nested spheres. \emph{Biometrika}, 99(3):551--568.
#'
#' Zoubouloglou, P., García-Portugués, E. and Marron, J. S. (2023). Scaled
#' torus principal component analysis. \emph{Journal of Computational and
#' Graphical Statistics}, 32(3):1024--1035.
#'
#' García-Portugués, E. and Meilán-Vila, A. (2025). Kernel density
#' estimation with polyspherical data and its applications. \emph{Journal of
#' the American Statistical Association}, to appear.
#' @importFrom stats as.dist cutree density hclust kmeans optim optimize
#'   prcomp quantile rnorm runif sd setNames
#' @importFrom graphics abline axis legend lines mtext par plot points
#'   polygon
#' @importFrom grDevices adjustcolor
#' @importFrom utils head
#' @keywords internal
"_PACKAGE"
