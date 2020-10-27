#include <vector>
#include <list_graph.h>

#include "Solver.h"

using namespace lemon;

std::vector<MyGraph::Node> apply2OptIteration(Solver& solver, const std::vector<MyGraph::Node>& solution, Rcpp::NumericMatrix* distanceMatrix = NULL);
