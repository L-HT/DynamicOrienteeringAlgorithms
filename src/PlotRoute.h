#include <Rcpp.h>
#include <vector>

#include "ProblemData.h"
#include "IO.h"

using namespace Rcpp;
using namespace lemon;


/*
 * Calls the plot function in R (from C++).
 */
void plotNetworkRoute_cpp(const ProblemData& pd, const std::vector<std::vector<MyGraph::Arc>>& totalPath = std::vector<std::vector<MyGraph::Arc>>());

/*
 * Calls the plot function from R if only a path of nodeIDs is given.
 */
//' @export
// [[Rcpp::export]]
void plotNetworkRoute_R(const Rcpp::DataFrame& nodeDf,
                        const Rcpp::DataFrame& arcDf,
                        const Rcpp::DataFrame& problemDf,
                        const std::vector<int>& nodeIDs);

