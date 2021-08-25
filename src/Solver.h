#ifndef SOLVER_H
#define SOLVER_H

#include <Rcpp.h>
#include <lemon/grid_graph.h>
#include <fstream>
#include <chrono>

#include "GraphToList.h"
#include "ProblemData.h"
#include "BudgetChangeHandling.h"

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
  std::vector<MyGraph::Node> initialSolutionBackup_;
  ResultData initialSolutionQuality_;
  
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
  
  ResultData evaluateSolutionMatrix(ProblemData& problemData,
                    const std::vector<MyGraph::Node>& solution,
                    std::string targetCriterion = "value",
                    bool forceLogging = false,
                    bool abortOnInvalidity = true);
      
  double getAndLogDistance(int i, int j);
  
  virtual bool terminationCriterionSatisfied();

  void readInitialSolutionFromFile(const std::string& pathToFile);
  void setInitialSolution(const std::vector<MyGraph::Node>& initialSolution);

  std::vector<MyGraph::Node> appendStartNode(std::vector<MyGraph::Node>& solution);
  std::vector<MyGraph::Node> asCompleteSolution(std::vector<MyGraph::Node> solution, bool left = true, bool right = true);
  void writeSolution(const std::vector<MyGraph::Node>& solution, bool isInitialSolution);
  

  ~Solver();
  
  /////////////////////////////////////////////////
  // Methods and members for handling changes in the budget
  // Implementation can be found in "BudgetChangeHandling.cpp"
  /////////////////////////////////////////////////
  
  int minBudget_;
  int maxBudget_;
  int tableSize_;
  int typeOfHandling_;
  
  static const int HANDLING_NONE = 0;
  static const int HANDLING_RESTART = 1;
  static const int HANDLING_TABLE = 2;
  static const int HANDLING_FIX_VIOLATION = 31;
  static const int HANDLING_FIX_SPARINGLY = 32;
  
  static const int REMOVAL_CRITERION_RATIO = 1;
  static const int REMOVAL_CRITERION_LENGTH = 2;
  static const int REMOVAL_CRITERION_VALUE = 3;
  
  std::map<int, EvaluatedSolution> bestSolutionsByBudget_;
  EvaluatedSolution getBestSolutionForBudget(int budget);
  
  void initializeSampledTable(int minBudget, int maxBudget, int tableSize);
  void updateTable(std::vector<MyGraph::Node> solution, ResultData resultData);
  // void fixInvalidSolutionViolation(std::vector<MyGraph::Node>& solution);
  // void fixInvalidSolutionSparingly(std::vector<MyGraph::Node>& solution);
  
  void repairSolutionByCriterion(int removeCriterion, std::vector<MyGraph::Node>& mySolution,
                                 ResultData& currentSolutionQuality);
  ResultData estimatePathLengthDecrease(const std::vector<int>& solutionIndices, 
                                        int i, 
                                        const std::vector<MyGraph::Node>& tempSolution,
                                        const ResultData& currentSolutionQuality);
  
  void printBudgetChangeHandlingInfo();
  virtual void resetSolver();
  
  /////////////////////////////////////////////////
  // Methods and members for restarting the algorithm
  /////////////////////////////////////////////////

};



#endif
