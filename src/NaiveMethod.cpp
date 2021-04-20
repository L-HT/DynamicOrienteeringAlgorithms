#include <Rcpp.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
#include <algorithm>
#include <ctime>

#include "Solver.h"
#include "PlotRoute.h"
#include "IO.h"
#include "ProblemData.h"
#include "CalculateLength.h"
#include "HelperFunctions.h"

#include "linkern_wrapper.h"
#include "DistanceMatrix.h"

using namespace lemon;


struct NaiveMethod : public Solver{

    // Rcpp::NumericMatrix distanceMatrix_;

    NaiveMethod(ProblemData& problemData,
              std::string logFileName, unsigned int runNumber, Rcpp::Environment rEnvironment,
              std::string targetCriterion, std::string fileSuffix, std::string pathToDistanceMatrix)

        : Solver(problemData, logFileName, runNumber, rEnvironment,
          targetCriterion, "greedy", fileSuffix, pathToDistanceMatrix) {

        // distanceMatrix_ = calculateDistanceMatrix(problemData);
    }

    /*
     * Naive method:
     * -go through destinations multiple times and always add the node that leads to the smallest
     * increase in total path length until the budget is exceeded. Then stop.
     */
    void run(){
        Rcpp::Rcout << "start " << algorithmName_ << std::endl;

        std::vector<MyGraph::Node> mySolution;
        bool possibleToAddNode = true;

        do{
            possibleToAddNode = false;
            std::vector<MyGraph::Node> tempSolution(mySolution.begin(), mySolution.end());
            MyGraph::Node closestNode = INVALID;

            double currentLength = std::numeric_limits<double>::max();
            double minLength = std::numeric_limits<double>::max();;

            for (MyGraph::Node n : problemData_.destinations_){
                if (std::find(tempSolution.begin(), tempSolution.end(), n) == tempSolution.end()){
                    additionalLogData_.numberOfTestedBitVectors_++;
                    tempSolution.push_back(n);

                    /*
                     * wenn zwischendurch bessere Lösung kommen, die aber das Verfahren dann nicht macht
                     */

                    ResultData res = evaluateSolutionMatrix(problemData_, tempSolution, targetCriterion_);
                    currentLength = res.length_;

                    if (currentLength <= problemData_.budget_ && currentLength < minLength){
                        closestNode = n;
                        minLength = currentLength;
                        possibleToAddNode = true;
                    }

                    tempSolution.pop_back();
                }
            }

            if (possibleToAddNode){
                mySolution.push_back(closestNode);
            }
            additionalLogData_.currentIteration_++;
        } while (possibleToAddNode);
        Rcpp::Rcout << bestSolutionQuality_ << std::endl;
        Rcpp::Rcout << "end " << algorithmName_ << std::endl;
    }

    ResultData evaluateSolutionMatrix(ProblemData& problemData,
                                      const std::vector<MyGraph::Node>& solution,
                                      std::string targetCriterion = "value",
                                      bool forceLogging = false,
                                      bool abortOnInvalidity = true){


        return evaluateSolution(problemData_, solution, targetCriterion, forceLogging, abortOnInvalidity, &distanceMatrix_);
    }
};



//' @export
// [[Rcpp::export]]
void naiveInitialSolution(const Rcpp::DataFrame& nodeDf, const Rcpp::DataFrame& arcDf, const Rcpp::DataFrame& problemDf, double budget,
                          std::string problemName = "", std::string fileSuffix = "", std::string pathToDistanceMatrix = ""){

    MyGraph graph;
    MyGraph::ArcMap<ArcData> arcMap(graph);
    MyGraph::NodeMap<NodeData> nodeMap(graph);

    MyGraph::Node startNode;
    std::vector<MyGraph::Node> destinations;

    // read distance matrix (if available)
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

    ProblemData pd(graph, nodeMap, arcMap, startNode, destinations, budget, problemName);

    Rcpp::Environment myEnvironment(MYPACKAGE);
    NaiveMethod nm(pd, problemName, 1, myEnvironment, "value", fileSuffix, pathToDistanceMatrix);
    nm.run();

    std::vector<MyGraph::Node> mySolution = nm.asCompleteSolution(nm.bestSolution_);

    // printNodeIdsOfVector(pd, mySolution);
    // std::vector<std::vector<MyGraph::Arc>> totalPath;
    // double value = 0;
    // double length = simulateTrip(pd, mySolution, totalPath, value);
    // Rcpp::Rcout << "length/value: " << length << "/" << value << std::endl;
    //
    // printNodeIdsOfVector(pd, mySolution);
    //
    // waitForInput("before plot", DEBUG_ENABLED);
    // plotNetworkRoute_cpp(pd, totalPath);

    nm.writeSolution(mySolution, true);

    // printNodeIdsOfVector(pd, mySolution);
    //
    // // waitForInput("Berechne Dist", true);
    // Rcpp::NumericMatrix dist = calculateDistanceMatrix(pd);
    // // waitForInput("Dist fertig. Jetzt LK...", true);
    // std::vector<MyGraph::Node> mySolution2 = callRepeatedLinKernighan(pd, mySolution, dist, nm);
    //
    // length = simulateTrip(pd, mySolution2, totalPath, value);
    // Rcpp::Rcout << "length2/value2: " << length << "/" << value << std::endl;
    //
    // //To-Do: Warum ändern sich die ausgewählten Knoten?
    // //forced logging
    // printNodeIdsOfVector(pd, mySolution2);
    //nm.evaluateSolution(nm.problemData_, mySolution2, nm.targetCriterion_, true, true);



    // waitForInput("end print...", true);
    // printNodeIdsOfVector(pd, mySolution2);
    // waitForInput("try to write2", true);
    // printNodeIdsOfVector(nm.problemData_, mySolution2);
    // waitForInput("try to write", true);

    // nm.writeSolution(mySolution2, false);
    //
    // waitForInput("after plot", DEBUG_ENABLED);
    // plotNetworkRoute_cpp(pd, totalPath);

}
