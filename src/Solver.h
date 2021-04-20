#ifndef SOLVER_H
#define SOLVER_H

#include <Rcpp.h>
#include <lemon/grid_graph.h>
#include <fstream>
#include <chrono>

#include "GraphToList.h"
#include "ProblemData.h"

using namespace lemon;

#define BUFFERSIZE 131072

struct Solver{

  ProblemData& problemData_;
  const std::string logFileName_;
  const std::string algorithmName_;

  const unsigned int runNumber_;

  char myBuffer_[BUFFERSIZE];
  std::ofstream logFile_;
  Rcpp::Environment rEnvironment_;

  std::string fileSuffix_;
  std::string targetCriterion_;

  std::vector<MyGraph::Node> initialSolution_;
  std::vector<MyGraph::Node> bestSolution_;

  std::chrono::high_resolution_clock::time_point startTime_;
  // double timeLimit_; // in Sekunden
  bool calledAsImprover_;

  AdditionalLogData additionalLogData_;
  ResultData bestSolutionQuality_;

  std::string terminationCriterion_;

  Rcpp::NumericMatrix distanceMatrix_;
  bool distanceMatrixIsGiven_;

  Solver(ProblemData& problemData,
         std::string logFileName,
         unsigned int runNumber,
         Rcpp::Environment rEnvironment,
         std::string targetCriterion,
         std::string algorithmName,
         std::string fileSuffix,
         std::string pathToDistanceMatrix);

  virtual ResultData evaluateSolution(ProblemData& problemData,
                   const std::vector<MyGraph::Node>& solution,
                   std::string targetCriterion = "value",
                   bool forceLogging = false,
                   bool abortOnInvalidity = true,
                   Rcpp::NumericMatrix* distanceMatrix = NULL);

  std::vector<MyGraph::Node> appendStartNode(std::vector<MyGraph::Node>& solution);

  virtual bool terminationCriterionSatisfied();

  void readInitialSolutionFromFile(const std::string& pathToFile);
  void setInitialSolution(const std::vector<MyGraph::Node>& initialSolution);

  std::vector<MyGraph::Node> asCompleteSolution(std::vector<MyGraph::Node> solution, bool left = true, bool right = true);

  void writeSolution(const std::vector<MyGraph::Node>& solution, bool isInitialSolution);
  ~Solver();
};

#endif
