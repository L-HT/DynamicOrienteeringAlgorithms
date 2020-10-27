#ifndef DISTANCE_MATRIX_H
#define DISTANCE_MATRIX_H

#include <Rcpp.h>
#include "ProblemData.h"

//To-Do: Distanzberechnung fehlt noch
Rcpp::NumericMatrix calculateDistanceMatrix(ProblemData pd);

void exportToDistanceMatrix(const Rcpp::DataFrame& nodeDf,
                       const Rcpp::DataFrame& arcDf,
                       const Rcpp::DataFrame& problemDf,
                       std::string problemName,
                       std::string fileSuffix = "",
                       std::string pathToChanges = "");

Rcpp::NumericMatrix readDistanceMatrix(std::string pathToMatrix);

bool distanceMatrixValid(const Rcpp::NumericMatrix& distanceMatrix);

#endif
