#include <Rcpp.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
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


struct RandomSolver : public Solver{

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> randomNumber; //Intervall [0,1)

    RandomSolver(ProblemData& problemData,
              std::string algorithmName,
              std::string logFileName, unsigned int runNumber, Rcpp::Environment rEnvironment,
              std::string targetCriterion, std::string fileSuffix, std::string pathToDistanceMatrix)

        : Solver(problemData, logFileName, runNumber, rEnvironment,
          targetCriterion, algorithmName, fileSuffix, pathToDistanceMatrix) {

        gen = std::mt19937(rd());
        randomNumber = std::uniform_real_distribution<double>(0.0, 1.0);

        // distanceMatrix_ = calculateDistanceMatrix(problemData);

        // Rcpp::stop("end");
    }

    /*
     * random solver: add all destinations to a solution; then remove destinations randomly until the path
     * length does not exceed the budget
     *
     * random improver: add random nodes until budget is exceeded, then remove random nodes until solution becomes
     * valid again
     */

    void run(){
        Rcpp::Rcout << "start " << algorithmName_ << std::endl;

        if (initialSolution_.empty()){
            Rcpp::Rcout << "start " << algorithmName_ << " (solver)" << std::endl;

            do {
                std::vector<MyGraph::Node> result(problemData_.destinations_.begin(), problemData_.destinations_.end());
                std::shuffle(result.begin(), result.end(), gen);
                ResultData currentResultQuality = evaluateSolutionMatrix(problemData_, result, "value");
                while(currentResultQuality.length_ > problemData_.budget_) {
                    int randomIndex = getRandomNumber(0, result.size()-1);
                    result.erase(result.begin() + randomIndex);
                    currentResultQuality = evaluateSolutionMatrix(problemData_, result, "value");
                }
            } while (bestSolutionQuality_.length_ <= 0);
            Rcpp::Rcout << bestSolutionQuality_ << std::endl;
        } else {
            Rcpp::Rcout << "start " << algorithmName_ << " (improver)" << std::endl;
            std::vector<MyGraph::Node> mySolution(initialSolution_.begin(), initialSolution_.end());
            ResultData currentResultQuality = evaluateSolutionMatrix(problemData_, mySolution, "value");

            std::vector<MyGraph::Node> nodesNotInSolution;
            for (MyGraph::Node n : problemData_.destinations_){
                if (std::find(mySolution.begin(), mySolution.end(), n) == mySolution.end()){
                    nodesNotInSolution.push_back(n);
                }
            }


            /*
             * 1) Remove nodes until solution becomes feasible;
             * 2) Add (randomly chosen) nodes at positions with minimum increase in length until budget is exceeded.
             * Repeat.
             */
            do {
                // remove nodes
                while(currentResultQuality.length_ > problemData_.budget_) {
                    int randomIndex = getRandomNumber(0, mySolution.size()-1);
                    nodesNotInSolution.push_back(*(mySolution.begin() + randomIndex));
                    mySolution.erase(mySolution.begin() + randomIndex);
                    currentResultQuality = evaluateSolutionMatrix(problemData_, mySolution, "value");
                }

                // add nodes
                while(currentResultQuality.length_ <= problemData_.budget_) {
                    int randomIndex = getRandomNumber(0, nodesNotInSolution.size()-1);

                    // just adding the nodes at random positions would not be very reasonable, I think...
                    // int randomPosition = getRandomNumber(0, mySolution.size());
                    // if (randomPosition == mySolution.size()){
                    //     mySolution.push_back(nodesNotInSolution[randomIndex]);
                    // } else {
                    //     mySolution.insert(mySolution.begin() + randomPosition);
                    // }

                    // so: greedy insertion

                    double shortestLength = std::numeric_limits<double>::max();
                    std::size_t minInsertIndex = -1;
                    for (std::size_t i = 0; i <= mySolution.size(); i++){
                        std::vector<MyGraph::Node> tempSolution(mySolution.begin(), mySolution.end());
                        if (i == mySolution.size()){
                            tempSolution.push_back(nodesNotInSolution[randomIndex]);
                        } else {
                            tempSolution.insert(tempSolution.begin() + i, nodesNotInSolution[randomIndex]);
                        }
                        ResultData tempResultQuality = evaluateSolutionMatrix(problemData_, tempSolution, "value");
                        if (tempResultQuality.length_ < shortestLength){
                            shortestLength = tempResultQuality.length_;
                            minInsertIndex = i;
                        }
                    }
                    if (minInsertIndex == mySolution.size()){
                        mySolution.push_back(nodesNotInSolution[randomIndex]);
                    } else {
                        mySolution.insert(mySolution.begin() + minInsertIndex, nodesNotInSolution[randomIndex]);
                    }

                    nodesNotInSolution.erase(nodesNotInSolution.begin() + randomIndex);
                    currentResultQuality = evaluateSolutionMatrix(problemData_, mySolution, "value");
                }
            } while (!terminationCriterionSatisfied());
        }
        evaluateSolutionMatrix(problemData_, bestSolution_, "value", true);
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
void callRandomSolver(const Rcpp::DataFrame& nodeDf, const Rcpp::DataFrame& arcDf, const Rcpp::DataFrame& problemDf,
                      double budget,
                      std::string algoName = "random",
                      std::string problemName = "", std::string fileSuffix = "", std::string pathToDistanceMatrix = ""){

    MyGraph graph;
    MyGraph::ArcMap<ArcData> arcMap(graph);
    MyGraph::NodeMap<NodeData> nodeMap(graph);

    MyGraph::Node startNode;
    std::vector<MyGraph::Node> destinations;

    /*
     * Verdammt, dieser Abschnitt muss etwas umstrukturiert werden:
     * -wenn Matrixpfad gegeben: Einlesen und checken, wie groß die Matrix wird und ob das eine gültige Matrix ist.
     *  -wenn gültige Matrix: Graphdata kann eingeschränkt werden
     *
     * -wenn Matrixpfad nicht gültig oder keine gültige Matrix: Dann wie normal einlesen
     *
     */

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
    RandomSolver rs(pd, algoName, problemName, 1, myEnvironment, "value", fileSuffix, pathToDistanceMatrix);
    rs.run();

    std::vector<MyGraph::Node> mySolution = rs.asCompleteSolution(rs.bestSolution_);

    rs.writeSolution(mySolution, true);
}

//' @export
// [[Rcpp::export]]
void callRandomImprover(const Rcpp::DataFrame& nodeDf, const Rcpp::DataFrame& arcDf, const Rcpp::DataFrame& problemDf,
                      double budget,
                      std::string algoName = "random",
                      std::string problemName = "",
                      std::string pathToInitialSolution ="",
                      std::string fileSuffix = "",
                      std::string pathToDistanceMatrix = ""){

    MyGraph graph;
    MyGraph::ArcMap<ArcData> arcMap(graph);
    MyGraph::NodeMap<NodeData> nodeMap(graph);

    MyGraph::Node startNode;
    std::vector<MyGraph::Node> destinations;

    /*
     * Verdammt, dieser Abschnitt muss etwas umstrukturiert werden:
     * -wenn Matrixpfad gegeben: Einlesen und checken, wie groß die Matrix wird und ob das eine gültige Matrix ist.
     *  -wenn gültige Matrix: Graphdata kann eingeschränkt werden
     *
     * -wenn Matrixpfad nicht gültig oder keine gültige Matrix: Dann wie normal einlesen
     *
     */
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

    Rcpp::Rcout << "fillDestinations...\n";
    fillDestinations(destinations, graph, nodeMap, problemDf);

    Rcpp::Rcout << "setStartNode...\n";
    setStartNode(startNode, graph, nodeMap, problemDf);

    ProblemData pd(graph, nodeMap, arcMap, startNode, destinations, budget, problemName);

    Rcpp::Environment myEnvironment(MYPACKAGE);

    Rcpp::Rcout << "construct...\n";
    RandomSolver rs(pd, algoName, problemName, 1, myEnvironment, "value", fileSuffix, pathToDistanceMatrix);

    rs.readInitialSolutionFromFile(pathToInitialSolution);
    rs.run();

    std::vector<MyGraph::Node> mySolution = rs.asCompleteSolution(rs.bestSolution_);

    rs.writeSolution(mySolution, false);
}
