#include <Rcpp.h>
#include <grid_graph.h>
#include <list_graph.h>
#include <stdlib.h>
#include <lgf_writer.h>
#include <lgf_reader.h>
#include <algorithm>
#include <random>

#include <chrono>
#include <thread>

#include "GraphToList.h"
#include "IO.h"
#include "HelperFunctions.h"
#include "CalculateLength.h"
#include "SimulatePath.h"

#include "ProblemData.h"

using namespace Rcpp;
using namespace lemon;

//To-Do: Distanzberechnung fehlt noch
Rcpp::NumericMatrix calculateDistanceMatrix(ProblemData pd){
    // pd.destinations_ does not contain the starting node, thus +1
    std::size_t n = pd.destinations_.size() + 1;
    Rcpp::NumericMatrix result_rcpp(n, n);

    //std::vector<std::vector<double>> result(pd.destinations_.size(), std::vector<double>(0.0, pd.destinations_.size()));
    // Rcpp::Rcout << "startNode: " << pd.nodeMap_[pd.startNode_].id_ << std::endl;
    // printNodeIdsOfVector(pd, pd.destinations_);

    for (std::size_t i = 0; i < n - 1; i++){
        // Rcpp::Rcout << "i+1: " << i+1 << ": " << pd.nodeMap_[pd.destinations_[i]].id_ << std::endl;

        for (std::size_t j = i + 1; j < n; j++){
            if (i == 0){
                // Rcpp::Rcout << "dist to: " << pd.nodeMap_[pd.startNode_].id_ << " -> "
                // << pd.nodeMap_[pd.destinations_[j-1]].id_ << std::endl;

                result_rcpp(i,j) = calculateShortestPathDistanceOnly(pd, pd.startNode_, pd.destinations_[j-1]);
            } else {
                // Rcpp::Rcout << "dist to: " << pd.nodeMap_[pd.destinations_[i-1]].id_ << " -> "
                //             << pd.nodeMap_[pd.destinations_[j-1]].id_ << std::endl;

                result_rcpp(i,j) = calculateShortestPathDistanceOnly(pd, pd.destinations_[i-1], pd.destinations_[j-1]);
            }

                //MinimumDetourDistanceOnly(pd.graph_, pd.destinations_[i], pd.destinations_[j], pd.forbiddenEdges_);
            result_rcpp(j,i) = result_rcpp(i,j);
        }
    }

    // for (std::size_t i = 0; i < n ; i++){
    //     if (i == 0){
    //         Rcpp::Rcout << "start: " << pd.nodeMap_[pd.startNode_].id_ << ": ";
    //     } else {
    //         Rcpp::Rcout << "i+1: " << i+1 << ": " << pd.nodeMap_[pd.destinations_[i-1]].id_ << ": ";
    //     }
    //     for (std::size_t j = 0;  j < n; j++){
    //         Rcpp::Rcout << result_rcpp(i,j) << " ";
    //     }
    //     Rcpp::Rcout << std::endl;
    // }
    return result_rcpp;
}

/*
* Read an instance with labeled nodes.
*/
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix distanceMatrixTest(std::string filePath){
    //std::string filePath = "./instances/1.lgf";

    MyGraph graph;
    MyGraph::NodeMap<NodeData> nodeMap(graph);
    MyGraph::NodeMap<bool> forbiddenNodes(graph, false);
    MyGraph::ArcMap<ArcData> arcMap(graph);
    MyGraph::ArcMap<int> forbiddenEdges(graph, 0);
    MyGraph::Node startNode;

    // readGridGraphFile(filePath, graph, nodeMap, arcMap, forbiddenNodes, forbiddenEdges, startNode);

    std::vector<MyGraph::Node> allDestinations;
    allDestinations.push_back(startNode);

    for (MyGraph::NodeIt it(graph); it != INVALID; ++it){
        if (nodeMap[it].type_ == 1){
            allDestinations.push_back(it);
        }
    }

    int budget = 25;
    ProblemData pd(graph, nodeMap, arcMap,
                   startNode, allDestinations, budget);

    return calculateDistanceMatrix(pd);
}

//' @export
// [[Rcpp::export]]
void exportToDistanceMatrix(const Rcpp::DataFrame& nodeDf,
                            const Rcpp::DataFrame& arcDf,
                            const Rcpp::DataFrame& problemDf,
                            std::string problemName,
                            std::string destination) {
    MyGraph graph;
    MyGraph::ArcMap<ArcData> arcMap(graph);
    MyGraph::NodeMap<NodeData> nodeMap(graph);
    MyGraph::Node startNode;
    std::vector<MyGraph::Node> destinations;

    Rcpp::Rcout << "readMyGraphIntoVariables...\n";
    readMyGraphIntoVariables(nodeDf, arcDf, graph, nodeMap, arcMap);

    Rcpp::Rcout << "readProblemDataIntoMaps...\n";
    readProblemDataIntoMaps(problemDf, graph, nodeMap, arcMap);

    Rcpp::Rcout << "fillDestinations...\n";
    fillDestinations(destinations, graph, nodeMap, problemDf);

    Rcpp::Rcout << "setStartNode...\n";
    setStartNode(startNode, graph, nodeMap, problemDf);

    Rcpp::Rcout << "construct...\n";
    ProblemData pd(graph, nodeMap, arcMap, startNode, destinations, 0, problemName);

    Rcpp::NumericMatrix result = calculateDistanceMatrix(pd);

    std::ofstream myFile;
    myFile.open(destination);
    for (std::size_t i = 0; i < result.nrow(); i++){
        for (std::size_t j = 0; j < result.ncol(); j++) {
            myFile << result(i,j);
            if (j < result.ncol()-1) {
                myFile << ",";
            }
        }
        if (i < result.nrow()-1) {
            myFile << "\n";
        }
    }
    myFile.close();

}

//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix readDistanceMatrix(std::string pathToMatrix) {
    /*
     * This code does not check if the file specified by pathToMatrix actually contains
     * valid information. If the file contains invalid data, it is well possible that
     * the program crashes.
     */
    std::ifstream myFile;
    std::string myLine;
    myFile.open(pathToMatrix);

    // assumption: matrix is quadratic
    Rcpp::NumericMatrix result;

    if (!myFile.fail()) {
        // parse first line in order to get matrix size
        std::getline(myFile, myLine);
        std::string delimiter = ",";
        std::vector<double> values;

        // Rcpp::Rcout << "start reading..." << "\n";
        int start = 0;
        int end = myLine.find(delimiter);
        while (end != std::string::npos) {
            values.push_back(std::stod(myLine.substr(start, end - start)));
            start = end + delimiter.length();
            end = myLine.find(delimiter, start);
        }
        values.push_back(std::stod(myLine.substr(start, end - start)));
        // printVector(values);

        std::size_t sizeOfMatrix = values.size();

        result = Rcpp::NumericMatrix(sizeOfMatrix, sizeOfMatrix);
        for (std::size_t i = 0; i < sizeOfMatrix; i++){
            result(0,i) = values[i];
        }

        int i = 1;
        int j = 0;

        Rcpp::Rcout << "Matrix with " << sizeOfMatrix << " elements per row constructed. Parse rest..." << "\n";
        // parse the rest of the matrix
        while (std::getline(myFile, myLine)) {
            std::string delimiter = ",";

            int start = 0;
            int end = myLine.find(delimiter);
            j = 0;
            while (end != std::string::npos){
                result(i, j++) = std::stod(myLine.substr(start, end - start));
                start = end + delimiter.length();
                end = myLine.find(delimiter, start);
            }
            result(i, j++) = std::stod(myLine.substr(start, end - start));
            i++;
        }

        // if (myFile.fail()) {
        //     Rcpp::Rcerr << "myFile.fail() returned true. Using empty matrix as distance matrix.\n";
        //     result = Rcpp::NumericMatrix(0,0);
        // }
        myFile.close();

        // print the matrix
        // Rcpp::Rcout << "print matrix with " << result.nrow() << " rows and " << result.ncol() << " columns:" << "\n";
        // for (i = 0; i < sizeOfMatrix; i++){
        //     for (j = 0; j < sizeOfMatrix; j++) {
        //         Rcpp::Rcout << result(i,j) << " ";
        //     }
        //     Rcpp::Rcout << "\n";
        // }
    } else {
        Rcpp::Rcerr << "Invalid distance matrix. Return empty matrix.\n";
        result = Rcpp::NumericMatrix(0,0);
    }
    return result;
}

bool distanceMatrixValid(const Rcpp::NumericMatrix& distanceMatrix) {
    return (distanceMatrix.ncol() != 0 && distanceMatrix.nrow() != 0);
}


