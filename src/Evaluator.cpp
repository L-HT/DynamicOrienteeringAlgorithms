// a simple evaluator that just evaluates a given solution

#include <Rcpp.h>
#include <lgf_writer.h>
#include <list_graph.h>
#include <algorithm>
#include <ctime>
#include <random>

#include "Solver.h"
#include "PlotRoute.h"
#include "IO.h"
#include "ProblemData.h"
#include "CalculateLength.h"
#include "HelperFunctions.h"

#include "linkern_wrapper.h"
#include "DistanceMatrix.h"

using namespace lemon;


struct Evaluator : public Solver{

    ResultData res_;

    Evaluator(ProblemData& problemData, Rcpp::Environment rEnvironment, std::string pathToDistanceMatrix)

        : Solver(problemData, "", 0, rEnvironment,
            "", "", "", pathToDistanceMatrix) {

    }

    void evaluate(const std::vector<int>& ids) {
        std::vector<MyGraph::Node> solution;
        Rcpp::Rcout << "searching nodes ... " << std::endl;
        for (const int& i : ids){
            for (MyGraph::Node& n : problemData_.destinations_){
                if (problemData_.nodeMap_[n].id_ == i) {
                    solution.push_back(n);
                    break;
                }
            }
        }
        Rcpp::Rcout << "evaluating: " ;
        printVector(ids);
        printNodeIdsOfVector(problemData_, solution);

        // Rcpp::Rcout << "print matrix with " << ev.distanceMatrix_.nrow() << " rows and " << ev.distanceMatrix_.ncol() << " columns:" << "\n";
        // for (int i = 0; i < ev.distanceMatrix_.nrow(); i++){
        //     for (int j = 0; j < ev.distanceMatrix_.ncol(); j++) {
        //         Rcpp::Rcout << ev.distanceMatrix_(i,j) << " ";
        //     }
        //     Rcpp::Rcout << "\n";
        // }
        res_ = getSolutionQuality(problemData_, solution, "value", false, &distanceMatrix_);

        // std::vector<MyGraph::Arc>& resultPath
    }
};



//' @export
// [[Rcpp::export]]
Rcpp::List callEvaluator(const Rcpp::DataFrame& nodeDf, const Rcpp::DataFrame& arcDf, const Rcpp::DataFrame& problemDf,
                      double budget, const std::vector<int>& indices, std::string pathToDistanceMatrix = ""){

    MyGraph graph;
    MyGraph::ArcMap<ArcData> arcMap(graph);
    MyGraph::NodeMap<NodeData> nodeMap(graph);

    MyGraph::Node startNode;
    std::vector<MyGraph::Node> destinations;

    Rcpp::NumericMatrix myMatrix = readDistanceMatrix(pathToDistanceMatrix);

    if (distanceMatrixValid(myMatrix)) {
        Rcpp::Rcout << "readDistanceMatrixIntoVariables...\n";
        readDistanceMatrixIntoVariables(nodeDf, arcDf, problemDf, graph, nodeMap, arcMap, myMatrix);
    } else {
        Rcpp::Rcout << "readMyGraphIntoVariables...\n";
        readMyGraphIntoVariables(nodeDf, arcDf, graph, nodeMap, arcMap);
        Rcpp::Rcout << "readProblemDataIntoMaps...\n";
        readProblemDataIntoMaps(problemDf, graph, nodeMap, arcMap);
    }

    fillDestinations(destinations, graph, nodeMap, problemDf);

    setStartNode(startNode, graph, nodeMap, problemDf);
    // destinations.push_back(startNode);

    ProblemData pd(graph, nodeMap, arcMap, startNode, destinations, budget, "");

    Rcpp::Environment myEnvironment(MYPACKAGE);
    Evaluator ev(pd, myEnvironment, pathToDistanceMatrix);

    // Rcpp::Rcout << "print matrix with " << ev.distanceMatrix_.nrow() << " rows and " << ev.distanceMatrix_.ncol() << " columns:" << "\n";
    // for (int i = 0; i < ev.distanceMatrix_.nrow(); i++){
    //     for (int j = 0; j < ev.distanceMatrix_.ncol(); j++) {
    //         Rcpp::Rcout << ev.distanceMatrix_(i,j) << " ";
    //     }
    //     Rcpp::Rcout << "\n";
    // }

    ev.evaluate(indices);

    return Rcpp::List::create(
        Rcpp::Named("length") = ev.res_.length_,
        Rcpp::Named("value") = ev.res_.value_
    );

}
