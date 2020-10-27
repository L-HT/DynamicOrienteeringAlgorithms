#ifndef LINKERN_WRAPPER_H
#define LINKERN_WRAPPER_H

#include <Rcpp.h>
#include "ProblemData.h"
#include "Solver.h"
#include "HelperFunctions.h"

std::vector<int> nodeVectorToIndexVector(const ProblemData& pd, const std::vector<MyGraph::Node>& solution);

std::vector<MyGraph::Node> callRepeatedLinKernighan(const ProblemData& pd, const std::vector<MyGraph::Node>& initialSolution,
                              Rcpp::NumericMatrix dist, Solver& solver);

#endif
