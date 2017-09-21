/*
 * Dirichlet Multinomial Clustering
 *
 * This approach to clustering counts models variation in depths while
 * essentially focusing on relative abundances, so might be a good compromise
 * between using relative abundances vs. raw data.
 *
 * reference: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0030126
 *
 * author: sankaran.kris@gmail.com
 * date: 09/20/2017
 */

data {
  int<lower=0> N; // Number of samples
  int<lower=0> P; // Number of species
  int<lower=0> K; // Number of clusters
  int<lower=0> x[N, P]; // species abundances
}
