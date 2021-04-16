#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <Rcpp.h>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <chrono>
#include <random>
#include <list_graph.h>

#include "GraphToList.h"
#include "ProblemData.h"

#define MYPACKAGE "package:DynamicOrienteeringAlgorithms"
#define PLOTFUNCTION "plotNetworkRoute"

#define OP_INSERT 1
#define OP_EDGEINSERT 2
#define OP_SWAP 3

#define DEBUG_ENABLED true
#define LOGINTERVAL 500

// altes LOGINTERVAL: 500
using namespace lemon;

std::vector<std::string> vectorToPermutation(std::vector<int> x);

struct AdditionalLogData{
    // numberOfShortestPathCalls_ increases by length of permutations
    // it thus favors paths with fewer nodes and algorithms that
    // successively build up a solution instead of working with existing solutions

    long numberOfShortestPathCalls_; // how often a shortest path between 2 nodes has been calculated (here it is similar to how often distances are calculated)
    long numberOfCalculatedPaths_; // how often the length of a given path has been calculated
    long numberOfTestedBitVectors_; // how many subsets of all destinations have been used so far
    int currentPhase_;
    int currentIteration_;

    AdditionalLogData():
        numberOfShortestPathCalls_(0), numberOfCalculatedPaths_(0), numberOfTestedBitVectors_(0),
        currentPhase_(0), currentIteration_(0){

    }
    AdditionalLogData(long numberOfShortestPathCalls, long numberOfCalculatedPaths, long numberOfTestedBitVectors,
                      int currentPhase, int currentIteration)
        : numberOfShortestPathCalls_(numberOfShortestPathCalls),
          numberOfCalculatedPaths_(numberOfCalculatedPaths),
          numberOfTestedBitVectors_(numberOfTestedBitVectors),
          currentPhase_(currentPhase), currentIteration_(currentIteration){

    }
    friend std::ostream& operator<<(std::ostream& os, const AdditionalLogData& ald){
        os << ald.numberOfCalculatedPaths_ << "," << ald.numberOfShortestPathCalls_ << "," << ald.numberOfTestedBitVectors_ << ",";
        os << ald.currentIteration_ << "," << ald.currentPhase_;
        return os;
    }
};

struct ResultData{
    double length_;
    double value_;
    ResultData() : length_(-1), value_(-1) {}
    ResultData(double length, double value) : length_(length), value_(value) {}

    friend std::ostream& operator<<(std::ostream& os, const ResultData& resultData){
        os << "(length: " << resultData.length_;
        os << ", value: " << resultData.value_ << ")";
        return os;
    }
    ResultData(const ResultData& other)
    {
        length_ = other.length_;
        value_ = other.value_;
    }
    ResultData& operator=(const ResultData& other)
    {
        length_ = other.length_;
        value_ = other.value_;
        return *this;
    }

    // "a > b" means that a is better than b
    bool operator<(const ResultData& b)
    {
        return (value_ < b.value_) || ((value_ == b.value_) && (length_ > b.length_));
    }
    bool operator>(const ResultData& b)
    {
        return (value_ > b.value_) || ((value_ == b.value_) && (length_ < b.length_));
    }
};

ResultData getAndLogSolutionQuality(ProblemData& problemData,
                             const std::vector<MyGraph::Node>& solution,
                             std::string targetCriterion,
                             std::ofstream& logFile,
                             std::chrono::high_resolution_clock::time_point startTime,
                             AdditionalLogData& ald,
                             ResultData& currentBest,
                             std::vector<MyGraph::Node>& currentBestSolution,
                             bool forceLogging = false,
                             bool abortOnInvalidity = true,
                             Rcpp::NumericMatrix* distanceMatrix = NULL);


ResultData getSolutionQuality(ProblemData& problemData,
                       const std::vector<MyGraph::Node>& solution,
                       std::string targetCriterion,
                       bool abortOnInvalidity = true,
                       Rcpp::NumericMatrix* distanceMatrix = NULL);

bool findArcInTotalPath(const MyGraph::Arc& arc, const std::vector<std::vector<MyGraph::Arc>>& totalPath, int& innerIndex,
                        int& outerIndex);

void printNodeIdsOfVector(const ProblemData& problemData, const std::vector<MyGraph::Node>& solution);

template <typename T>
void printVector (const std::vector<T>& v) {
  Rcpp::Rcout << ">>";
  for (typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); it++){
    Rcpp::Rcout << *it << ",";
  }
  Rcpp::Rcout << "<<" << std::endl;
}

template <typename T>
void printList (const std::list<T>& v) {
  Rcpp::Rcout << "]]";
  for (typename std::list<T>::const_iterator it = v.begin(); it != v.end(); it++){
    Rcpp::Rcout << *it << ",";
  }
  Rcpp::Rcout << "[[" << std::endl;
}
template <typename T, typename C>
bool contains (const C& v, T value) {
  return (std::find(v.begin(), v.end(), value) != v.end());
}

int getRandomNumber(int min, int max);
void waitForInput(std::string msg, bool enabled);


template <typename T>
void insertOperator(std::vector<T>& permutation, int takePosition, int insertPosition){

    T valueToInsert = permutation[takePosition];
    if (takePosition < insertPosition){
        for (int i = takePosition; i <= (insertPosition - 1); i++){
            permutation[i] = permutation[i+1];
        }
        permutation[insertPosition] = valueToInsert;
    }
    if (insertPosition < takePosition){
        for (int i = takePosition; i >= (insertPosition + 1); i--){
            permutation[i] = permutation[i-1];
        }
        permutation[insertPosition] = valueToInsert;
    }
}

template <typename T>
void swapOperator(std::vector<T>& permutation, int position1, int position2){
    T temp = permutation[position1];
    permutation[position1] = permutation[position2];
    permutation[position2] = temp;
}

template <typename T>
void edgeInsertOperator(std::vector<T>& permutation, int takePosition, int insertPosition){
    T value1ToInsert = permutation[takePosition];
    T value2ToInsert = permutation[takePosition + 1];
    if (takePosition < insertPosition){
        for (int i = takePosition; i <= (insertPosition - 1); i++){
            permutation[i] = permutation[i+2];
        }

        permutation[insertPosition] = value1ToInsert;
        permutation[insertPosition+1] = value2ToInsert;
    }
    if (insertPosition < takePosition){
        for (int i = takePosition+1; i >= (insertPosition + 2); i--){
            permutation[i] = permutation[i-2];
        }
        permutation[insertPosition] = value1ToInsert;
        permutation[insertPosition+1] = value2ToInsert;
    }
}
// void insertOperator(std::vector<int>& permutation, int takePosition, int insertPosition);
// void swapOperator(std::vector<int>& permutation, int position1, int position2);
// void edgeInsertOperator(std::vector<int>& permutation, int takePosition, int insertPosition);
// int replaceOperator(std::vector<int>& permutation, int pos, std::size_t numberOfJobs, int forbiddenJob = -1);

void writeSolutionIntoFile(const ProblemData& pd, const std::vector<MyGraph::Node>& solution, const std::string& fileName);

int findIndexInCumulativeProbabilities(std::vector<double> cumulativeProbabilities, double value);

#endif
