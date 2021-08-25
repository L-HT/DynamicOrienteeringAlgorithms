
#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <sstream>
//#include <ctime>
#include <limits>
#include <chrono>

#include "HelperFunctions.h"
#include "Solver.h"
#include "ProblemData.h"
#include "IDSearchFunctions.h"
#include "DistanceMatrix.h"

Solver::Solver(ProblemData& problemData,
         std::string logFileName,
         unsigned int runNumber,
         Rcpp::Environment rEnvironment,
         std::string targetCriterion,
         std::string algorithmName,
         std::string fileSuffix,
         std::string pathToDistanceMatrix)
    : problemData_(problemData), logFileName_(logFileName),
      algorithmName_(algorithmName),
      runNumber_(runNumber),
      rEnvironment_(rEnvironment),
      fileSuffix_(fileSuffix),
      targetCriterion_(targetCriterion),
      calledAsImprover_(false),
      terminationCriterion_(""){

    std::stringstream ss;
    ss << "./output/" << logFileName_ << "-" << algorithmName_ << "-" << runNumber_;
    if (fileSuffix != ""){
      ss << "-" << fileSuffix;
    }
    Rcpp::Rcout << "output to: " << ss.str() << std::endl;

    // if (targetCriterion_ == "value"){
    //     bestSolutionQuality_ = std::numeric_limits<int>::min();
    // }
    // timeLimit_ = 600;

    logFile_.rdbuf()->pubsetbuf(myBuffer_, BUFFERSIZE);
    logFile_.open (ss.str(), std::ios_base::trunc); //app für append

    logFile_ << "calculatedPaths,shortestPathCalls,testedBitVectors,iteration,phase,time,pathLength,value" << "\n";
    Rcpp::Rcout << "Time measurement starts -now-!" << std::endl;
    startTime_ = std::chrono::high_resolution_clock::now();

    additionalLogData_ = AdditionalLogData();

    // path to distance matrix: if empty, then it is calculated from scratch
    // if the variable points to a valid file containing a distance matrix, it is read from that file
    // otherwise (invalid path or "none") nothing happens
    if (pathToDistanceMatrix == "") {
        Rcpp::Rcout << "distance matrix: calculated\n";
        distanceMatrix_ = calculateDistanceMatrix(problemData_);
        additionalLogData_.numberOfCalculatedPaths_ += (problemData_.destinations_.size() + 1)*(problemData_.destinations_.size() + 1)/2;
    } else {
        Rcpp::Rcout << "distance matrix: read\n";
        distanceMatrix_ = readDistanceMatrix(pathToDistanceMatrix);
    }
    // Rcpp::Rcout << "print matrix2 with " << distanceMatrix_.nrow() << " rows and " << distanceMatrix_.ncol() << " columns:" << "\n";
    // for (int i = 0; i < distanceMatrix_.nrow(); i++){
    //     for (int j = 0; j < distanceMatrix_.ncol(); j++) {
    //         Rcpp::Rcout << distanceMatrix_(i,j) << " ";
    //     }
    //     Rcpp::Rcout << "\n";
    // }
    
    minBudget_ = 0;
    maxBudget_ = 0;
    tableSize_ = 0;
    typeOfHandling_ = 0;
}

std::vector<MyGraph::Node> Solver::appendStartNode(std::vector<MyGraph::Node>& solution){
    std::vector<MyGraph::Node> result;
    result.push_back(problemData_.destinations_[0]);
    result.insert(result.end(), solution.begin(), solution.end());
    result.push_back(problemData_.destinations_[0]);
    return result;
}

ResultData Solver::evaluateSolution(ProblemData& problemData,
                     const std::vector<MyGraph::Node>& solution,
                     std::string targetCriterion,
                     bool forceLogging,
                     bool abortOnInvalidity,
                     Rcpp::NumericMatrix* distanceMatrix){

    /*
     * check for changes in the problem instance:
     * Is there a change in the queue? Implement it and remove it from queue.
     * Update bestSolutionQuality_.
     * Check if newQuality is better. If yes, then it will be rewritten.
     */


    int oldBudget = problemData_.budget_;
    bool dataHasChanged = checkForAndImplementChanges(problemData, additionalLogData_, startTime_);
    bool budgetHasChanged = oldBudget != problemData_.budget_;
    
    if (dataHasChanged) {

        bestSolutionQuality_ = evaluateSolution(problemData_, bestSolution_, targetCriterion, false, true, distanceMatrix);
        
        if (budgetHasChanged) {
            // What if this solution is invalid after the change?
            // set bestSolution=c() and bestSolutionQuality_=0
            if (bestSolutionQuality_.length_ > problemData.budget_) {
                Rcpp::Rcout << "Best solution has become invalid!" << "\n";

                switch(typeOfHandling_){
                    case(HANDLING_FIX_VIOLATION): // HANDLING BY FIXING THE VIOLATION
                        repairSolutionByCriterion(REMOVAL_CRITERION_LENGTH, bestSolution_, bestSolutionQuality_);
                        Rcpp::Rcout << "Handling by HANDLING_FIX_VIOLATION: " << bestSolutionQuality_ << "\n";
                        break;
                    case(HANDLING_FIX_SPARINGLY): // HANDLING BY REMOVING SMALLEST VALUE ELEMENTS
                        repairSolutionByCriterion(REMOVAL_CRITERION_VALUE, bestSolution_, bestSolutionQuality_);
                        Rcpp::Rcout << "Handling by HANDLING_FIX_SPARINGLY: " << bestSolutionQuality_ << "\n";
                        break;
                    case(HANDLING_TABLE): // HANDLING BY TABLE
                        // invalidate the current solution; the actual fixing is done below
                        bestSolutionQuality_ = ResultData();
                        bestSolution_.clear();
                        break;
                    case(HANDLING_RESTART): // HANDLING BY RESTARTING THE ALGORITHM
                        resetSolver();
                        break;
                    default:
                        bestSolutionQuality_ = ResultData();
                        bestSolution_.clear();
                        break;
                }
            }
            
            // HANDLING BY TABLE (put here so that this is also called when the budget increases)
            // this code also deals with the case when the budget increases
            if (typeOfHandling_ == HANDLING_TABLE) {
                EvaluatedSolution bestEvaluatedSolution = getBestSolutionForBudget(problemData.budget_);
                if (bestSolutionQuality_ < bestEvaluatedSolution.resultData_) {
                    // if bestSolution_ is still valid and better, then do not change the solution
                    // if bestSolution_ is invalid, this code is always executed
                    bestSolution_ = bestEvaluatedSolution.solution_;
                    bestSolutionQuality_ = bestEvaluatedSolution.resultData_;
                    Rcpp::Rcout << "Handling by table: " << bestSolutionQuality_ << "\n";
                }
            }
            forceLogging = true;
            
        }
    }

    ResultData oldBestQuality = bestSolutionQuality_;
    ResultData newQuality = getAndLogSolutionQuality(problemData, solution,
                                                    targetCriterion,
                                                    logFile_, startTime_, additionalLogData_,
                                                    bestSolutionQuality_, bestSolution_,
                                                    forceLogging, abortOnInvalidity, distanceMatrix);

    // update table if budget changes are handled with a table
    if (typeOfHandling_ == HANDLING_TABLE) {
        updateTable(solution, newQuality);
    }
    
    // write best solution into a file every time it is updated
    if (newQuality > oldBestQuality && newQuality.length_ <= problemData.budget_){
        // writeSolution(solution, !calledAsImprover_);
    }
    return newQuality;
}
ResultData Solver::evaluateSolutionMatrix(ProblemData& problemData,
                                  const std::vector<MyGraph::Node>& solution,
                                  std::string targetCriterion,
                                  bool forceLogging,
                                  bool abortOnInvalidity){
    return evaluateSolution(problemData, solution, targetCriterion, forceLogging, abortOnInvalidity, &distanceMatrix_);
}

/*
 * default criterion: absolute time
 */
bool Solver::terminationCriterionSatisfied(){


    if (problemData_.terminationCriterion_ == "evaluation") {
        return additionalLogData_.numberOfCalculatedPaths_ >= 1000000;
    }
    if (problemData_.terminationCriterion_ == "bitVector") {
        return additionalLogData_.numberOfTestedBitVectors_ >= 1000000;
    }
    if (problemData_.terminationCriterion_ == "distanceEvaluation") {
        return additionalLogData_.numberOfShortestPathCalls_ >= 1000000;
    }

    // std::chrono::high_resolution_clock::time_point currentTime = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> timeSpan = std::chrono::duration_cast<std::chrono::duration<double>>(currentTime - startTime_);
    // double timePassed = timeSpan.count();
    //
    // return timePassed >= timeLimit_;
    // return numberOfEvaluation_ >= maxEvaluation_;
    return additionalLogData_.numberOfCalculatedPaths_ >= 1000000
        && additionalLogData_.numberOfTestedBitVectors_ >= 1000000
        && additionalLogData_.numberOfShortestPathCalls_ >= 1000000;
}
Solver::~Solver(){
    if (typeOfHandling_ == HANDLING_TABLE) {
        std::map<int, EvaluatedSolution>::iterator it;
        for (it = bestSolutionsByBudget_.begin(); it != bestSolutionsByBudget_.end(); it++){
            Rcpp::Rcout << "B=" << it->first << ": " << it->second.resultData_ << "\n";
        }
    }
    
    logFile_.close();
}

void Solver::setInitialSolution(const std::vector<MyGraph::Node>& initialSolution){
    initialSolution_.assign(initialSolution.begin(), initialSolution.end());
}
void Solver::readInitialSolutionFromFile(const std::string& pathToFile){
    std::vector<int> idVector;
    std::ifstream myFile(pathToFile);
    std::string line;
    std::getline(myFile, line);

    std::string delimiter = ",";
    std::size_t start = 0;
    std::size_t  end = line.find(delimiter, start);

    while (end != std::string::npos){
        std::string tempString = line.substr(start, end - start);
        idVector.push_back(std::stoi(tempString));

        start = end + delimiter.length();
        end = line.find(delimiter, start);
    }

    // the range start-end contains the last element
    std::string tempString = line.substr(start, end);
    idVector.push_back(std::stoi(tempString));

    initialSolution_.clear();

    // without first and last element (because they are the initial node)
    for (std::size_t i = 1; i < idVector.size() - 1; i++){
        // Rcpp::Rcout << idVector[i] << ", ";
        MyGraph::Node n = getNodeFromInternalID(problemData_.graph_, problemData_.nodeMap_, idVector[i]);
        // initialSolution_.push_back(problemData_.graph_.nodeFromId(idVector[i]));
        initialSolution_.push_back(n);
    }
}

std::vector<MyGraph::Node> Solver::asCompleteSolution(std::vector<MyGraph::Node> solution, bool left, bool right){
    if (left){
        solution.insert(solution.begin(), problemData_.startNode_);
    }
    if (right){
        solution.push_back(problemData_.startNode_);
    }
    return solution;
}

void Solver::writeSolution(const std::vector<MyGraph::Node>& solution, bool isInitialSolution){

    std::stringstream ss;
    if (isInitialSolution){
        ss << "./initialSolutions/";
    } else {
        ss << "./solutions/";
    }
    if (!fileSuffix_.empty()){
        ss << logFileName_ << "_" << fileSuffix_ << "_" << algorithmName_;
    } else {
        ss << logFileName_ << "_" << algorithmName_;
    }

    writeSolutionIntoFile(problemData_, solution, ss.str());
}

double Solver::getAndLogDistance(int i, int j){
    additionalLogData_.numberOfShortestPathCalls_++;
    return distanceMatrix_(i,j);
}

void Solver::resetSolver() {
    // don't do anything
}
// void Solver::writeBestSolution(){
//   std::ofstream myfile;
//   std::stringstream ss;
//   ss << "./solutions/" << logFileName_ << "-" << algorithmName_ << "-" << runNumber_;
//   if (fileSuffix_ != ""){
//     ss << "-" << fileSuffix_;
//   }
//   myfile.open (ss.str(), std::ios_base::trunc);
//
//   for (std::size_t i = 0; i < numberOfJobs_; i++){
//     myfile << "j" << bestSolution_[i];
//     if (i < numberOfJobs_ - 1){
//       myfile << ",";
//     }
//   }
//   myfile << std::endl;
//   for (std::size_t i = 0; i < numberOfJobs_; i++){
//     myfile << "j" << bestSolutionM2_[i];
//     if (i < numberOfJobs_ - 1){
//       myfile << ",";
//     }
//   }
//   myfile << std::endl;
//
//   // wenn nicht-Permutationsplan gefunden wurde
//   if (bestSolutionM2_ != bestSolution_){
//     myfile << "!";
//   }
//   myfile.close();
//   // printVector(bestSolution_);
//   // printVector(bestSolutionM2_);
// }

