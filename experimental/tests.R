#tests

library(rotasym)
library(polykde)
library(waddR)
library(sphunif)
library(shapes)
library(smacof)
library(sdetorus)
library(rgl)
library(MASS)
library(sdetorus)
library(rgl)
library(plot3D)
library(testthat)
library(usethis)

source("C:/Users/minise/Downloads/STPCA/Script1_Setup_and_Functions.R")

set.seed(100)
n <- 90
d <- 3
Sigma <- 0.1 * diag(c(1, 1,1,1))
cluster1 <- mvrnorm(n = n / 3, mu = c(-1, -2,.3,1), Sigma = Sigma)
cluster2 <- mvrnorm(n = n / 3, mu = c(3, 0.4,1,-3), Sigma = Sigma)
cluster3 <- mvrnorm(n = n / 3, mu = c(-0.8, 2.5,-1,0), Sigma = Sigma)
cluster4 <- mvrnorm(n = n / 3, mu = c(1, 2.5,-.5,-2), Sigma = Sigma)
cluster5 <- mvrnorm(n = n / 3, mu = c(.3,-.5,.6,1), Sigma = Sigma)
cluster6 <- mvrnorm(n = n / 3, mu = c(2, 1.1,-1.1,2), Sigma = Sigma)
clusters <- rbind(cluster1, cluster2, cluster3,cluster4,cluster5,cluster6)
clusters <- (clusters + pi) %% (2 * pi) - pi
a = c(.1,.4,-.9)
b =c(-.3,.3,.6)
# test_that("angles to polysphere", {

expect_equal(anglesToSphere(c(pi,0)),c(-1, 0,0))
expect_equal(length(anglesToSphere(c(pi,0))),3)
expect_equal(dim(anglesToSphere(clusters)),c(180,5))
expect_equal(dim(polyAnglesToSphere(clusters,c(4))),c(180,5))
expect_equal(dim(polyAnglesToSphere(clusters,c(3,1))),c(180,6))
expect_equal(dim(polyAnglesToSphere(clusters,c(2,2))),c(180,6))
expect_equal(dim(polyAnglesToSphere(clusters,c(1,1,2))),c(180,7))
expect_equal(dim(polyAnglesToSphere(clusters,c(1,1,1,1))),c(180,8))
# })
# test_that("polysphere distance", {
expect_equal(polysphere_distance(a,b,dims=c(1,1,1))-torus_distance(a,b),0)
# })


