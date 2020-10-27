#include <Rcpp.h>
#include <string>
#include <vector>
#include <iostream>

#include "CalculateLength.h"
#include "PlotRoute.h"
#include "VectorConversion.h"

using namespace Rcpp;
using namespace lemon;


void plotNetworkRoute_cpp(const ProblemData& pd, const std::vector<std::vector<MyGraph::Arc>>& totalPath){
    Environment e(MYPACKAGE);
    Function plotNetworkRoute = e[PLOTFUNCTION];
    // waitForInput("list1", DEBUG_ENABLED);
    List graphAsList = routeGraphToList_R(pd, totalPath);
    // waitForInput("list2", DEBUG_ENABLED);
    DataFrame nodeDf = Rcpp::as<Rcpp::DataFrame>(graphAsList["nodeDf"]);
    DataFrame arcDf = Rcpp::as<Rcpp::DataFrame>(graphAsList["arcDf"]);
    plotNetworkRoute(R_NilValue, nodeDf, arcDf);
}



void plotNetworkRoute_R(const Rcpp::DataFrame& nodeDf,
                        const Rcpp::DataFrame& arcDf,
                        const Rcpp::DataFrame& problemDf,
                        const std::vector<int>& nodeIDs){

    /*
     * nodeIDs must contain the starting node at the beginning and at the end
     */

    MyGraph graph;
    MyGraph::ArcMap<ArcData> arcMap(graph);
    MyGraph::NodeMap<NodeData> nodeMap(graph);
    MyGraph::Node startNode;
    std::vector<MyGraph::Node> destinations;

    // read distance matrix (if available)
    // Rcpp::NumericMatrix myMatrix = readDistanceMatrix(pathToDistanceMatrix);
    // if (distanceMatrixValid(myMatrix)) {
    //     Rcpp::Rcout << "readDistanceMatrixIntoVariables...\n";
    //     readDistanceMatrixIntoVariables(nodeDf, arcDf, problemDf, graph, nodeMap, arcMap, myMatrix);
    // } else {
    //     Rcpp::Rcout << "readMyGraphIntoVariables...\n";
    //     readMyGraphIntoVariables(nodeDf, arcDf, graph, nodeMap, arcMap);
    //     Rcpp::Rcout << "readProblemDataIntoMaps...\n";
    //     readProblemDataIntoMaps(problemDf, graph, nodeMap, arcMap);
    // }

    readMyGraphIntoVariables(nodeDf, arcDf, graph, nodeMap, arcMap);
    readProblemDataIntoMaps(problemDf, graph, nodeMap, arcMap);
    fillDestinations(destinations, graph, nodeMap, problemDf);
    setStartNode(startNode, graph, nodeMap, problemDf);

    ProblemData pd(graph, nodeMap, arcMap, startNode, destinations, 0, "");

    Environment e(MYPACKAGE);
    Function plotNetworkRoute = e[PLOTFUNCTION];

    std::vector<MyGraph::Node> mySolution;



    for (int i : nodeIDs){
        for (MyGraph::Node& n : pd.destinations_){

            if (pd.nodeMap_[n].id_ == i){
                Rcpp::Rcout << "push-back: " << pd.nodeMap_[n].id_ << "\n";
                mySolution.push_back(n);//pd.graph_.nodeFromId(i));
                break;
            }
        }
        if (i == pd.nodeMap_[pd.startNode_].id_) {
            Rcpp::Rcout << "push-back: " << i << "\n";
            mySolution.push_back(pd.startNode_);
        }
    }


    std::vector<std::vector<MyGraph::Arc>> totalPath;
    double value = 0;
    double length = simulateTrip(pd, mySolution, totalPath, value);
    Rcpp::Rcout << "length/value: " << length << "/" << value << std::endl;
    plotNetworkRoute_cpp(pd, totalPath);
}

