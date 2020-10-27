#include <Rcpp.h>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <chrono>
#include <fstream>
#include <sstream>
#include <random>

#include "HelperFunctions.h"
#include "CalculateLength.h"
#include "ProblemData.h"
#include "VectorConversion.h"

std::vector<std::string> vectorToPermutation(std::vector<int> x){
  std::vector<std::string> result;
  result.reserve(x.size());
  for (int i : x){
    result.push_back("j" + std::to_string(i));
  }
  return result;
}



//Grenzen sind mit drin
int getRandomNumber(int min, int max){
  int result = min + (rand() % static_cast<int>(max - min + 1));
  //Rcpp::Rcout << result << " ";
  return result;
}



void waitForInput(std::string msg, bool enabled){
  if (enabled){
    Rcpp::Environment baseEnvironment = Rcpp::Environment("package:base");
    Rcpp::Function readLineR = baseEnvironment["readline"];
    Rcpp::Rcout << msg << " (press Enter)" << std::endl;
    readLineR("");
  }
}


// note: abortOnInvalidity is currently unused
ResultData getSolutionQuality(ProblemData& problemData,
                       const std::vector<MyGraph::Node>& solution,
                       std::string targetCriterion,
                       bool abortOnInvalidity,
                       Rcpp::NumericMatrix* distanceMatrix){

    double value = 0;
    double length = 0;

    std::vector<MyGraph::Node> completeSolution;
    completeSolution.push_back(problemData.startNode_);
    completeSolution.insert(completeSolution.end(), solution.begin(), solution.end());
    completeSolution.push_back(problemData.startNode_);

    //std::vector<std::vector<MyGraph::Arc>> totalPath;

    // here a distance matrix of size 1x1 indicates that no distance matrix is given
    // (hab auf die Schnelle keinen besseren Weg gefunden)
    if (distanceMatrix->ncol()==0 && distanceMatrix->nrow()==0) {

        length = simulateTripDistanceOnly(problemData,
                                  completeSolution, value);
        // Rcpp::Rcout << "use manual path: " << length << "\n";
    } else {



        for (const MyGraph::Node& n : solution){
            value += problemData.nodeMap_[n].value_;
            // Rcpp::Rcout << "value: " << value << std::endl;
        }

        std::vector<int> indexVector = nodeVectorToIndexVector(problemData, completeSolution);
        for (std::size_t i = 0; i < indexVector.size() - 1; i++){
            // Rcpp::Rcout << "length: " << length << std::endl;
            length += (*distanceMatrix)(indexVector[i], indexVector[i+1]);
        }

        // Rcpp::Rcout << "use matrix path: " << length << "\n";

    }


    // Rcpp::Rcout << std::endl;

    // if (length > problemData.budget_){
        // return ResultData(-1, -1);
    // } else {
    return ResultData(length, value);

    // }
}


ResultData getAndLogSolutionQuality(ProblemData& problemData,
                             const std::vector<MyGraph::Node>& solution,
                             std::string targetCriterion,
                             std::ofstream& logFile,
                             std::chrono::high_resolution_clock::time_point startTime,
                             AdditionalLogData& ald,
                             ResultData& currentBest,
                             std::vector<MyGraph::Node>& currentBestSolution,
                             bool forceLogging,
                             bool abortOnInvalidity,
                             Rcpp::NumericMatrix* distanceMatrix){

    std::chrono::high_resolution_clock::time_point currentTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> timeSpan = std::chrono::duration_cast<std::chrono::duration<double>>(currentTime - startTime);

    double timePassed = timeSpan.count();

    ResultData resultData = getSolutionQuality(problemData,
                                    solution,
                                    targetCriterion,
                                    abortOnInvalidity,
                                    distanceMatrix);


    // update bestSolution
    if (targetCriterion == "value"){ // maximize value
        if (resultData.length_ <= problemData.budget_){

            if (resultData.value_ > currentBest.value_){
                // Rcpp::Rcout << "update best solution value: " << currentBest << " -> " << resultData << std::endl;
                currentBest = resultData;
                currentBestSolution.assign(solution.begin(), solution.end());
            } else {
                if (resultData.value_ == currentBest.value_ && resultData.length_ < currentBest.length_){
                    // Rcpp::Rcout << "update best solution length: " << currentBest << " -> " << resultData << std::endl;
                    currentBest = resultData;
                    currentBestSolution.assign(solution.begin(), solution.end());
                }
            }

        }

    } else { // minimize path length (if the problem is not orienteering, but finding a shortest path or a tour)
        if (resultData.length_ < currentBest.length_&& resultData.length_ <= problemData.budget_){
            currentBest = resultData;
            currentBestSolution.assign(solution.begin(), solution.end());
        }
    }

    // if force logging, then the currentSolution is logged
    if (forceLogging){
        // Rcpp::Rcout << "f " << ald.numberOfCalculatedPaths_ << " ";
        logFile << ald << "," <<  timePassed << "," << currentBest.length_ << "," << currentBest.value_ << "\n";
    } else {
        if (ald.numberOfCalculatedPaths_ % LOGINTERVAL == 0){
            // Rcpp::Rcout << ald.numberOfCalculatedPaths_ << " ";
            logFile << ald << "," <<  timePassed << "," << currentBest.length_ << "," << currentBest.value_ << "\n";
        }
    }
    ald.numberOfCalculatedPaths_++;

    // Evaluation of a solution: A permutation of length n has
    // n-1 intermediate routes plus 1 (going back to the starting node)
    ald.numberOfShortestPathCalls_ += solution.size() + 1;

    // Rcpp::Rcout << "result: " << resultData << std::endl;

    return resultData;
}

bool findArcInTotalPath(const MyGraph::Arc& arc, const std::vector<std::vector<MyGraph::Arc>>& totalPath, int& innerIndex,
                        int& outerIndex){
    for (int outer = 0; outer < (int) totalPath.size(); outer++){
        for (int inner = 0; inner < (int) totalPath[outer].size(); inner++){
            if (arc == totalPath[outer][inner]){
                outerIndex = outer;
                innerIndex = inner;
                return true;
            }
        }
    }
    outerIndex = -1;
    innerIndex = -1;
    return false;
}

void writeSolutionIntoFile(const ProblemData& pd, const std::vector<MyGraph::Node>& solution, const std::string& fileName){
    std::ofstream outputFile(fileName);
    for (std::size_t i = 0; i < solution.size(); i++){
        // std::stringstream ss;
        // ss << "this is " << i << std::endl;
        // waitForInput(ss.str(), true);
        outputFile << pd.nodeMap_[solution[i]].id_;//pd.graph_.id(solution[i]);
        if (i < solution.size() - 1){
            outputFile << ",";
        }
    }
    outputFile.close();
}

void printNodeIdsOfVector(const ProblemData& problemData, const std::vector<MyGraph::Node>& solution){
    for (MyGraph::Node n : solution){
        Rcpp::Rcout << problemData.nodeMap_[n].id_ << ","; //problemData.graph_.id(n) << ",";
    }
    Rcpp::Rcout << std::endl;
}

int findIndexInCumulativeProbabilities(std::vector<double> cumulativeProbabilities, double value){

    bool valueFound = false;
    int counter = 0;
    do{
        if (value < cumulativeProbabilities[counter]){
            valueFound = true;
            return counter;
        } else {
            counter++;
        }

    } while (!valueFound && counter < (int) cumulativeProbabilities.size());
    std::stringstream ss;
    ss << "Could not find value " << value << " in cumulative probabilities." << std::endl;
    Rcpp::stop(ss.str());
    return -1;
}

