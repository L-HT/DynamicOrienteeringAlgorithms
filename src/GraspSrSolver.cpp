#include <Rcpp.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
#include <algorithm>
#include <random>

#include "Solver.h"
#include "PlotRoute.h"
#include "IO.h"
#include "ProblemData.h"
#include "CalculateLength.h"
#include "HelperFunctions.h"
#include "SimpleMethods.h"
#include "DistanceMatrix.h"

using namespace lemon;

struct GraspSrSolver : public Solver{
    std::vector<MyGraph::Node> currentSolution_;
    ResultData currentSolutionRes_;
    
    // a second "best" solution which is needed in the algorithm but works differently from the best solution defined by Solver
    std::vector<MyGraph::Node> tempBestSolution_;
    ResultData tempBestSolutionRes_;

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> randomNumber; //Intervall [0,1)

    // Rcpp::NumericMatrix distanceMatrix_;

    GraspSrSolver(ProblemData& problemData,
                  std::string logFileName, unsigned int runNumber, Rcpp::Environment rEnvironment,
                  std::string targetCriterion, std::string fileSuffix, std::string pathToDistanceMatrix)

        : Solver(problemData, logFileName, runNumber, rEnvironment,
          targetCriterion, "gsr", fileSuffix, pathToDistanceMatrix) {


        gen = std::mt19937(rd());
        randomNumber = std::uniform_real_distribution<double>(0.0, 1.0);

        // distanceMatrix_ = calculateDistanceMatrix(problemData_);

        currentSolution_ = std::vector<MyGraph::Node>();
        currentSolutionRes_.length_ = 0;
        currentSolutionRes_.value_ = 0;
    }

    void run(){

        Rcpp::Rcout << "start " << algorithmName_  <<  std::endl;
        printBudgetChangeHandlingInfo();
        if (initialSolution_.empty() && !calledAsImprover_){
            Rcpp::Rcout << "start " << algorithmName_ << " (Solver)" << std::endl;
            // waitForInput("beforeAddVertex", DEBUG_ENABLED);

            // while (addVertex()) {
            //     // Rcpp::Rcout << "current: ";
            //     // printNodeIdsOfVector(problemData_, currentSolution_);
            //     additionalLogData_.currentIteration_++;
            // }
            callAddVertexRepeatedly();
            if (currentSolutionRes_.length_ > problemData_.budget_) {
                Rcpp::stop("Error1: addVertex() ended with invalid solution!");
            }
            tempBestSolution_.assign(currentSolution_.begin(), currentSolution_.end());
            tempBestSolutionRes_ = currentSolutionRes_;
            // check bestSolution vs. tempBestSolution
            if (tempBestSolutionRes_ > bestSolutionQuality_ && tempBestSolutionRes_.length_ < problemData_.budget_) {
                bestSolution_.assign(tempBestSolution_.begin(), tempBestSolution_.end());
                bestSolutionQuality_ = tempBestSolutionRes_;
            }
            printNodeIdsOfVector(problemData_, currentSolution_);

            // waitForInput("beforeLS", DEBUG_ENABLED);
            additionalLogData_.currentPhase_++;
            while (localSearch()) {
                Rcpp::Rcout << bestSolutionQuality_ << std::endl;
                additionalLogData_.currentIteration_++;
            }

            if (tempBestSolutionRes_.value_ != bestSolutionQuality_.value_ || tempBestSolutionRes_.length_ != bestSolutionQuality_.length_){
                double toleranceValue = bestSolutionQuality_.value_ * 0.001;
                double toleranceLength = bestSolutionQuality_.length_ * 0.001;
                if (std::abs(tempBestSolutionRes_.value_ - bestSolutionQuality_.value_) > toleranceValue ||
                    std::abs(tempBestSolutionRes_.length_ - bestSolutionQuality_.length_) > toleranceLength){

                    Rcpp::Rcout << "tempBest / bestSolution: " << tempBestSolutionRes_ << " / " << bestSolutionQuality_ << std::endl;
                    Rcpp::Rcout << "--WARNING: tempBest and bestSolution are rather different." << std::endl;
                } else {
                    Rcpp::Rcout << "tempBest / bestSolution: " << tempBestSolutionRes_ << " / " << bestSolutionQuality_ << std::endl;
                    Rcpp::Rcout << "--warning: tempBest and bestSolution are slightly different." << std::endl;
                }
            }
        } else {
            Rcpp::Rcout << "start " << algorithmName_ << " (improvement heuristic)" << std::endl;
            
            if (initialSolution_.empty()){
                Rcpp::Rcout << "GSR was called as improvement heuristic, but without initial solution. Call addVertex()...\n";
                // if currentSolution_ is still empty after the while-loop, it is
                // possible that no valid solution exists, but 
                // the budget might increase over time, so keep searching
                callAddVertexRepeatedly();
                if (currentSolutionRes_.length_ > problemData_.budget_) {
                    Rcpp::stop("Error2: callAddVertexRepeatedly() ended with invalid solution!");
                }
            } else {
                // Rcpp::Rcout << "assignments\n";
                currentSolution_.assign(initialSolution_.begin(), initialSolution_.end());
                currentSolutionRes_ = evaluateSolutionMatrix(problemData_, currentSolution_);
                bestSolution_.assign(initialSolution_.begin(), initialSolution_.end());
                bestSolutionQuality_ = currentSolutionRes_;
                // if (currentSolutionRes_.length_ > problemData_.budget_) {
                //     Rcpp::stop("Error3: Initial solution is invalid!");
                // }
            }
            
            tempBestSolution_.assign(currentSolution_.begin(), currentSolution_.end());
            tempBestSolutionRes_ = currentSolutionRes_;
            // printNodeIdsOfVector(problemData_, currentSolution_);
            /*
            *  anytime-variant for GRASP-SR:
            * -removes a random segment from the path
            */
            do {
                // waitForInput("start of big-do",true);
                // Rcpp::Rcout << "big do\n";
                // Rcpp::Rcout << "sizetempbest: " << tempBestSolution_.size() << "\n";
                while (localSearch()) {additionalLogData_.currentIteration_++;}

                if (currentSolution_.size() > 1){
                    int r1, r2;
                    do {
                        r1 = getRandomNumber(0, currentSolution_.size());
                        r2 = getRandomNumber(0, currentSolution_.size());

                        if (r1 > r2) {
                            std::swap(r1,r2);
                        }
                        // it is not allowed that all nodes from the solution are removed
                    } while (r1 == 0 && r2 == currentSolution_.size());

                    currentSolution_.erase(currentSolution_.begin() + r1, currentSolution_.begin() + r2);
                } else {
                    Rcpp::Rcerr << "Local search yielded empty solution or just a single node. Call addVertex()..." << std::endl;
                    // while (addVertex()) {
                    //     // Rcpp::Rcout << "current: ";
                    //     // printNodeIdsOfVector(problemData_, currentSolution_);
                    //     additionalLogData_.currentIteration_++;
                    // }
                    callAddVertexRepeatedly();
                    if (currentSolutionRes_.length_ > problemData_.budget_) {
                        Rcpp::stop("Error4: callAddVertexRepeatedly() (in subcase) ended with invalid solution!");
                    }
                    tempBestSolution_.assign(currentSolution_.begin(), currentSolution_.end());
                    tempBestSolutionRes_ = currentSolutionRes_;
                    // check bestSolution vs. tempBestSolution
                    if (tempBestSolutionRes_ > bestSolutionQuality_ && tempBestSolutionRes_.length_ < problemData_.budget_) {
                        bestSolution_.assign(tempBestSolution_.begin(), tempBestSolution_.end());
                        bestSolutionQuality_ = tempBestSolutionRes_;
                    }
                    printNodeIdsOfVector(problemData_, currentSolution_);
                }
            } while (!terminationCriterionSatisfied());
        }
        evaluateSolutionMatrix(problemData_, bestSolution_, "value", true);
        Rcpp::Rcout << "end " << algorithmName_ << std::endl;
    }

    bool addVertex(MyGraph::Node blockedNode = lemon::INVALID){
        // Paper: length of path (with start and end node) - 1; here: cS.size() is only path without start/end;
        // thus: (cS.size()+2)-1 equals l as found in paper; ergo +1 necessary
        // std::vector<MyGraph::Node> path_r(currentSolution.begin(), currentSolution.end());
        std::size_t l = currentSolution_.size()+1;

        std::vector<std::vector<MyGraph::Node>> candidateList;
        std::vector<double> candidateProfits;

        double bestImprovement = 0;
        // waitForInput("before W-loop", DEBUG_ENABLED);

        additionalLogData_.currentPhase_ = 1;
        for (MyGraph::Node n : problemData_.destinations_){
            // Rcpp::Rcout << "search " << problemData_.nodeMap_[n].id_ << " in ";
            // printNodeIdsOfVector(problemData_, currentSolution_);

            additionalLogData_.numberOfTestedBitVectors_++;
            additionalLogData_.currentPhase_ = 11;
            // waitForInput("search node", DEBUG_ENABLED);
            if (n != blockedNode && std::find(currentSolution_.begin(), currentSolution_.end(), n) == currentSolution_.end()){
                std::size_t bestInsertPosition;
                double minPathLength = std::numeric_limits<double>::max();
                ResultData tempBestRes;
                for (std::size_t insertPos = 0; insertPos < currentSolution_.size() + 1; insertPos++){
                    // Rcpp::Rcout << insertPos << " / " << currentSolution_.size() + 1 << std::endl;
                    // waitForInput("inner w-loop part 1", DEBUG_ENABLED);
                    std::vector<MyGraph::Node> tempSolution(currentSolution_.begin(), currentSolution_.begin() + insertPos);
                    tempSolution.push_back(n);

                    tempSolution.insert(tempSolution.end(), currentSolution_.begin() + insertPos, currentSolution_.end());
                    ResultData res = evaluateSolutionMatrix(problemData_, tempSolution);

                    if (res.length_ < minPathLength){
                        minPathLength = res.length_;
                        bestInsertPosition = insertPos;
                        tempBestRes = res;
                    }
                }

                additionalLogData_.currentPhase_ = 12;
                std::vector<MyGraph::Node> solution1(currentSolution_.begin(), currentSolution_.begin() + bestInsertPosition);
                solution1.push_back(n);
                solution1.insert(solution1.end(), currentSolution_.begin() + bestInsertPosition, currentSolution_.end());

                if (minPathLength <= problemData_.budget_){
                    // Rcpp::Rcout << "pushback easy: ";
                    // printNodeIdsOfVector(problemData_, solution1);
                    // Rcpp::Rcout << tempBestRes << std::endl;

                    // because of dynamics, it should be checked that this solution is actually better than the current solution

                    currentSolutionRes_ = evaluateSolutionMatrix(problemData_, currentSolution_);
                    if (tempBestRes.value_ > currentSolutionRes_.value_){

                        candidateList.push_back(solution1);
                        candidateProfits.push_back(tempBestRes.value_);

                    }
                } else {
                    std::size_t j = 0;
                    for (std::size_t i = 0; i < l; i++){
                        if (j < i) {
                            j = i;
                        }

                        // std::vector<MyGraph::Node> solution2(solution1.begin(), solution1.end());
                        // ResultData res;
                        // do {
                        //     solution2.erase(solution2.begin() + i);
                        //     res = evaluateSolutionMatrix(problemData_, solution2);
                        // } while (res.length_ > problemData_.budget_ && solution2.begin() + i != solution2.end());

                        std::vector<MyGraph::Node> solution2;
                        ResultData res;
                        // Rcpp::Rcout << "vector hier: ";
                        // printNodeIdsOfVector(problemData_, solution1);
                        do {
                            j++;
                            solution2.assign(solution1.begin(), solution1.begin()+i);
                            if (j < l){
                                solution2.insert(solution2.end(), solution1.begin()+j, solution1.end());
                            } else {
                                j = l;
                            }

                            // Rcpp::Rcout << "(i,j,l): " << i << ", " << j << ", " << l << std::endl;
                            additionalLogData_.numberOfTestedBitVectors_++;
                            res = evaluateSolutionMatrix(problemData_, solution2);
                        } while (res.length_ > problemData_.budget_ && j != l);

                        // printNodeIdsOfVector(problemData_, solution2);
                        // Rcpp::Rcout << "res: " << res << std::endl;

                        if (res.length_ <= problemData_.budget_) {
                            currentSolutionRes_ = evaluateSolutionMatrix(problemData_, currentSolution_);
                            if (res.value_ > currentSolutionRes_.value_ || (res.value_ == currentSolutionRes_.value_ && res.length_ < currentSolutionRes_.length_)){
                                // if (!solution2.empty()){
                                    // Rcpp::Rcout << "Push back." << std::endl;
                                    // Rcpp::Rcout << res << " ; " << currentSolutionRes_ << std::endl;
                                    candidateList.push_back(solution2);
                                    candidateProfits.push_back(res.value_);

                                // }
                            }
                        }
                        // waitForInput("check", true);
                    }
                }
            }
        }

        // Rcpp::Rcout << std::endl;

        additionalLogData_.currentPhase_ = 13;
        if (!candidateList.empty()){
            // Rcpp::Rcout << "candidates: " << std::endl;
            // for (std::vector<MyGraph::Node> s : candidateList){
            //     printNodeIdsOfVector(problemData_, s);
            // }
            // Rcpp::Rcout << "end cl: " << std::endl;
            if (candidateList.size() != candidateProfits.size()){
                Rcpp::stop("candidateList and candidateProfits have different number of elements. This should not happen...");
            }

            // waitForInput("candidates shown", DEBUG_ENABLED);

            double bestProfitDifference = 0;

            for (double d : candidateProfits){
                if (d - currentSolutionRes_.value_ > bestProfitDifference){
                    bestProfitDifference = d - currentSolutionRes_.value_;
                }
            }
            // Rcpp::Rcout << "bestDiff: " << bestProfitDifference << std::endl;
            // waitForInput("profitDiffs", DEBUG_ENABLED);

            std::vector<std::vector<MyGraph::Node>> restrictedCandidateList;
            for (std::size_t i = 0; i < candidateList.size(); i++){

                if (bestProfitDifference==0){
                    Rcpp::Rcout << "zero for: ";
                    printNodeIdsOfVector(problemData_, candidateList[i]);
                    Rcpp::Rcout << "cand_ " << i << ": " << candidateProfits[i] << std::endl;
                    Rcpp::Rcout << candidateProfits[i] - currentSolutionRes_.value_ << " >= " << 0.2 * bestProfitDifference << " ?" << std::endl;
                }
                if (candidateProfits[i] - currentSolutionRes_.value_ >= 0.2 * bestProfitDifference){
                    restrictedCandidateList.push_back(candidateList[i]);
                }
            }
            // Rcpp::Rcout << "r-candidates: " << std::endl;
            // for (std::vector<MyGraph::Node> s : restrictedCandidateList){
            //     printNodeIdsOfVector(problemData_, s);
            // }
            // Rcpp::Rcout << "end rcl " << std::endl;

            // waitForInput("rcl", DEBUG_ENABLED);

            if (restrictedCandidateList.size() == 0){
                 Rcpp::Rcout << "candidates: ";
                 for (std::size_t i = 0; i < candidateList.size(); i++){
                     printNodeIdsOfVector(problemData_, candidateList[i]);
                     Rcpp::Rcout << candidateProfits[i] << std::endl << std::endl;
                 }
                 // for (std::vector<MyGraph::Node> s : candidateList){
                 //
                 // }
                 Rcpp::Rcout << "bestProfitDifference: " << bestProfitDifference << std::endl;
                 Rcpp::Rcout << "r-candidates: " << std::endl;
                 for (std::vector<MyGraph::Node> s : restrictedCandidateList){
                     printNodeIdsOfVector(problemData_, s);
                 }
                 Rcpp::Rcout << "end rcl " << std::endl;
                 Rcpp::Rcout << "old quality: " << currentSolutionRes_ << std::endl;
                 Rcpp::Rcout << "new quality: " << evaluateSolutionMatrix(problemData_, currentSolution_);

                 Rcpp::Rcout << "INFO: There were so many changes, the candidate list is likely unusable (=worse then current solution). Restart local search..." << std::endl;
                 // Rcpp::stop("rcl empty");

                 // return true so that the local search is restarted (because there were too many changes)
                 return true;
             }

            additionalLogData_.currentPhase_ = 14;
            int r = getRandomNumber(0, restrictedCandidateList.size() - 1);
            // Rcpp::Rcout << "chosen: ";
            // printNodeIdsOfVector(problemData_, restrictedCandidateList[r]);
            // waitForInput("before sub-2opt", DEBUG_ENABLED);
            currentSolution_ = apply2OptIteration(*this, restrictedCandidateList[r], &distanceMatrix_);

            // Rcpp::Rcout << "after 2opt: ";
            // printNodeIdsOfVector(problemData_, currentSolution_);
            currentSolutionRes_ = evaluateSolutionMatrix(problemData_, currentSolution_);

            // waitForInput("after sub-2opt", DEBUG_ENABLED);
            return true;
        }
        return false;
    }

    bool localSearch(){
       
        if (tempBestSolution_.empty()){
            callAddVertexRepeatedly();
            if (currentSolutionRes_.length_ > problemData_.budget_) {
                Rcpp::Rcout << "length / budget: " << currentSolutionRes_.length_ << " / " << problemData_.budget_ << "\n";
                Rcpp::stop("Error5: callAddVertexRepeatedly() (in localSearch()) ended with invalid solution!");
            }
            tempBestSolution_.assign(currentSolution_.begin(), currentSolution_.end());
            tempBestSolutionRes_ = currentSolutionRes_;
        }
        
        std::size_t l = tempBestSolution_.size() - 1;
        bool improved = false;

        std::vector<MyGraph::Node> backUpPath(tempBestSolution_.begin(), tempBestSolution_.end());

        for (std::size_t i = 0; i < l; i++){
            additionalLogData_.currentPhase_ = 21; // before 2-OPT in local search
            Rcpp::Rcout << "i / l: " << i << " / " << l << std::endl;

            // waitForInput("LS loop", DEBUG_ENABLED);
            currentSolution_.assign(backUpPath.begin(), backUpPath.end());
            MyGraph::Node forbiddenNode = currentSolution_[i];

            // waitForInput("erase", DEBUG_ENABLED);
            currentSolution_.erase(currentSolution_.begin() + i);
            // printNodeIdsOfVector(problemData_, currentSolution_);

            // waitForInput("2opt", DEBUG_ENABLED);
            currentSolution_ = apply2OptIteration(*this, currentSolution_, &distanceMatrix_);


            currentSolutionRes_ = evaluateSolutionMatrix(problemData_, currentSolution_);
            // printNodeIdsOfVector(problemData_, currentSolution_);

            // waitForInput("addVertex", DEBUG_ENABLED);
            additionalLogData_.currentPhase_ = 22; // calling addVertex
            while (addVertex(forbiddenNode)){}

            // printNodeIdsOfVector(problemData_, currentSolution_);
            // waitForInput("iffy", DEBUG_ENABLED);
            if (currentSolutionRes_  > tempBestSolutionRes_ && currentSolutionRes_.length_ <= problemData_.budget_){
                // bestSolution_ is not overwritten here because this already happens during the calls to evaluateSolution
                tempBestSolution_.assign(currentSolution_.begin(), currentSolution_.end());
                tempBestSolutionRes_ = currentSolutionRes_;
                improved = true;
            }
            if (terminationCriterionSatisfied()){
                Rcpp::Rcout << "Termination criterion satisfied. Algorithm stopped." << std::endl;
                return false;
            }
            Rcpp::checkUserInterrupt();
        }
        return improved;
    }
    
    void callAddVertexRepeatedly() {
        // if currentSolution_ is still empty after the while-loop, it is
        // possible that no valid solution exists, but 
        // the budget might increase over time, so keep searching
        do{
            while (addVertex()) {
                // Rcpp::Rcout << "current: ";
                // printNodeIdsOfVector(problemData_, currentSolution_);
                additionalLogData_.currentIteration_++;
            }
        } while ((currentSolutionRes_.length_ > problemData_.budget_ || currentSolution_.empty()) && !terminationCriterionSatisfied());
        
        // terminationCriterion satisfied, but still empty solution
        // -> no valid solution exists 
        if (terminationCriterionSatisfied()) {
            Rcpp::Rcout << "currentSolutionQuality_:" << currentSolutionRes_ << "\n"; 
            Rcpp::Rcout << "Final budget: " << problemData_.budget_ << "\n";
            Rcpp::Rcout << additionalLogData_ << "\n";
            Rcpp::stop("Termination criterion satisfied, but no valid solution was found with addVertex(). It is possible that no valid solution exists.");
        }
    }
    void resetSolver() override{
        if (initialSolutionBackup_.empty()) {
            initialSolutionBackup_.assign(initialSolution_.begin(), initialSolution_.end());
            initialSolutionQuality_ =  evaluateSolutionMatrix(problemData_, initialSolutionBackup_);
        }
        if (initialSolutionQuality_.length_ > problemData_.budget_){
            initialSolution_.clear();
            bestSolution_.clear();
            bestSolutionQuality_ = ResultData();
            tempBestSolution_.clear();
            tempBestSolutionRes_ = ResultData();
            currentSolution_.clear();
            currentSolutionRes_ = ResultData();
            Rcpp::Rcout << "Restart algorithm with empty initial solution.\n";
        } else {
            initialSolution_.assign(initialSolutionBackup_.begin(), initialSolutionBackup_.end());
            bestSolution_.assign(initialSolution_.begin(), initialSolution_.end());
            bestSolutionQuality_ = initialSolutionQuality_;
            tempBestSolution_.assign(initialSolution_.begin(), initialSolution_.end());
            tempBestSolutionRes_ = bestSolutionQuality_;
            currentSolution_.assign(initialSolution_.begin(), initialSolution_.end());
            currentSolutionRes_ = bestSolutionQuality_;
            Rcpp::Rcout << "Restart algorithm with the given initial solution.\n";
        }
        Rcpp::Rcout << "Current log data: " << additionalLogData_ << "\n";
        run(); 
        
        std::vector<MyGraph::Node> mySolution = asCompleteSolution(bestSolution_);
        writeSolution(mySolution, !calledAsImprover_);
        Rcpp::Rcout << "Log data at the end:" << additionalLogData_ << "\n";
        Rcpp::stop("Algorithm used restarts, but the termination criterion is now satisfied. Algorithm is stopped...");
    }
    
    ResultData evaluateSolution(ProblemData& problemData,
                                const std::vector<MyGraph::Node>& solution,
                                std::string targetCriterion,
                                bool forceLogging,
                                bool abortOnInvalidity,
                                Rcpp::NumericMatrix* distanceMatrix) override{
        
        int oldBudget = problemData_.budget_;
        bool dataHasChanged = checkForAndImplementChanges(problemData, additionalLogData_, startTime_);
        bool budgetHasChanged = oldBudget != problemData_.budget_;
        
        if (dataHasChanged) {
            
            bestSolutionQuality_ = evaluateSolution(problemData_, bestSolution_, targetCriterion, false, true, distanceMatrix);
            
            if (budgetHasChanged) {
                // What if this solution is invalid after the change?
                // set bestSolution=c() and bestSolutionQuality_=0
                if (bestSolutionQuality_.length_ > problemData.budget_) {
                    Rcpp::Rcout << "Best solution has become invalid!!!" << "\n";
                    
                    /*
                     * TODO: Wird tempBestSolution manipuliert, obwohl die Lösung ungültig ist?
                     */
                    switch(typeOfHandling_){
                    case(HANDLING_FIX_VIOLATION): // HANDLING BY FIXING THE VIOLATION
                        repairSolutionByCriterion(REMOVAL_CRITERION_LENGTH, tempBestSolution_, tempBestSolutionRes_);
                        Rcpp::Rcout << "handle tempBestSolution_: " << tempBestSolutionRes_ << "\n";
                        repairSolutionByCriterion(REMOVAL_CRITERION_LENGTH, bestSolution_, bestSolutionQuality_);
                        Rcpp::Rcout << "Handling by HANDLING_FIX_VIOLATION: " << bestSolutionQuality_ << "\n";
                        break;
                    case(HANDLING_FIX_SPARINGLY): // HANDLING BY REMOVING SMALLEST VALUE ELEMENTS
                        repairSolutionByCriterion(REMOVAL_CRITERION_VALUE, tempBestSolution_, tempBestSolutionRes_);
                        repairSolutionByCriterion(REMOVAL_CRITERION_VALUE, bestSolution_, bestSolutionQuality_);
                        // repairSolutionByCriterion(REMOVAL_CRITERION_VALUE, currentSolution_, currentSolutionRes_);
                        // evaluateSolutionMatrix(problemData_, tempBestSolution_);
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
                        currentSolutionRes_ = ResultData();
                        currentSolution_.clear();
                        tempBestSolutionRes_ = ResultData();
                        tempBestSolution_.clear();
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
                        bestSolution_.assign(bestEvaluatedSolution.solution_.begin(), bestEvaluatedSolution.solution_.end());
                        bestSolutionQuality_ = bestEvaluatedSolution.resultData_;
                        tempBestSolution_.assign(bestEvaluatedSolution.solution_.begin(), bestEvaluatedSolution.solution_.end());
                        tempBestSolutionRes_ = bestEvaluatedSolution.resultData_;
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
        
        // if (tempBestSolutionRes_ > bestSolutionQuality_ && tempBestSolutionRes_.length_ < problemData_.budget_) {
        //     bestSolution_.assign(tempBestSolution_.begin(), tempBestSolution_.end());
        //     // Rcpp::Rcout << "temp/best/budget: " << tempBestSolutionRes_.length_ << "/" << bestSolutionQuality_.length_ << "/" <<
        //     //     problemData_.budget_ << "\n";
        //     bestSolutionQuality_ = tempBestSolutionRes_;
        // }
        // write best solution into a file every time it is updated
        if (newQuality > oldBestQuality && newQuality.length_ <= problemData.budget_){
            // writeSolution(solution, !calledAsImprover_);
        }
        return newQuality;
    }

    // ResultData evaluateSolutionMatrix(ProblemData& problemData,
    //                                   const std::vector<MyGraph::Node>& solution,
    //                                   std::string targetCriterion = "value",
    //                                   bool forceLogging = false,
    //                                   bool abortOnInvalidity = true){
    // 
    //     ResultData res = evaluateSolution(problemData_, solution, targetCriterion, forceLogging, abortOnInvalidity, &distanceMatrix_);
    //     return res;//evaluateSolution(problemData_, solution, targetCriterion, forceLogging, abortOnInvalidity, &distanceMatrix_);
    // }
};


//' @export
// [[Rcpp::export]]
void callGraspSrSolver(const Rcpp::DataFrame& nodeDf,
                       const Rcpp::DataFrame& arcDf,
                       const Rcpp::DataFrame& problemDf,
                       double budget,
                       std::string problemName,
                       int runNumber = 0,
                       std::string fileSuffix = "",
                       std::string pathToChanges = "",
                       std::string pathToDistanceMatrix = ""){

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
    if (pathToChanges != "") {
        readDynamicChangesFromFile(pd, pathToChanges);
    }
    /*
    * Calculate a naive starting solution and write it into initialSolutions.
    * Careful Method: Add the closest customer
    */

    Rcpp::Environment myEnvironment(MYPACKAGE);
    GraspSrSolver gsrs(pd, problemName, runNumber, myEnvironment, "value", fileSuffix, pathToDistanceMatrix);
    gsrs.run();


    /////////////////////////////

    std::vector<MyGraph::Node> mySolution = gsrs.asCompleteSolution(gsrs.bestSolution_, true, true);
    gsrs.writeSolution(mySolution, true);

}

//' @export
// [[Rcpp::export]]
void callGraspSrImprover(const Rcpp::DataFrame& nodeDf,
                         const Rcpp::DataFrame& arcDf,
                         const Rcpp::DataFrame& problemDf,
                         double budget,
                         std::string problemName,
                         unsigned int runNumber,
                         std::string pathToInitialSolution = "",
                         std::string fileSuffix = "",
                         std::string pathToChanges = "",
                         std::string pathToDistanceMatrix ="",
                         int budgetChangeHandlingMode = 0,
                         int minBudgetToHandle = 0,
                         int maxBudgetToHandle = 0,
                         int budgetChangeTableSize = 0){

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
    readDynamicChangesFromFile(pd, pathToChanges);
    /*
     * Calculate a naive starting solution and write it into initialSolutions.
     * Careful Method: Add the closest customer
     */

    Rcpp::Environment myEnvironment(MYPACKAGE);
    GraspSrSolver gsrs(pd, problemName, runNumber, myEnvironment, "value", fileSuffix, pathToDistanceMatrix);
    
    gsrs.typeOfHandling_ = budgetChangeHandlingMode;
    if (budgetChangeHandlingMode == Solver::HANDLING_TABLE) {
        gsrs.initializeSampledTable(minBudgetToHandle, maxBudgetToHandle, 
                                   budgetChangeTableSize);
    }
    gsrs.calledAsImprover_ = true;
    if (pathToInitialSolution != ""){
        gsrs.readInitialSolutionFromFile(pathToInitialSolution);
    }
    // waitForInput("beforeRun", DEBUG_ENABLED);
    gsrs.run();


    /////////////////////////////

    std::vector<MyGraph::Node> mySolution = gsrs.asCompleteSolution(gsrs.bestSolution_);
    gsrs.writeSolution(mySolution, false);

}
