#include <Rcpp.h>

#include <iostream>
#include <limits>
#include <cmath>

#include "Solver.h"
#include "VectorConversion.h"

EvaluatedSolution Solver::getBestSolutionForBudget(int budget) {
    /*
     * die Lösung, deren Länge noch am tiefsten unter dem Argument liegt
     * edge-cases: 
     *  it zeigt ans Ende: budget ist mehr als genug, 
     *      nimm das letzte gueltige Element
     *  it zeigt an den Anfang: es gibt keine gueltige Loesung
     *      nimm leere Loesung
     */
    std::map<int,EvaluatedSolution>::iterator itUpper = bestSolutionsByBudget_.upper_bound(budget);
    EvaluatedSolution result;
    
    /*
     * itUpper points to beginning: no valid solution -> return empty solution
     *  -also contains edge case where bestSolutionsByBudget_ is empty (end=begin)
     *  
     *  otherwise: there exists at least 1 valid solution, namely itUpper-1
     *      -also contains edge case where itUpper points to end (and end != begin)
     */
    if (itUpper != bestSolutionsByBudget_.begin()) {
        itUpper--;
        result = (itUpper->second);
        // result. assign(itUpper->second.solution_.begin(), itUpper->second.solution_.end());
    } else {
        Rcpp::Rcout << "Info: There are no stored solutions that are valid for budget=" << budget << ".\n";
    }
    return result;
}

// void Solver::fixInvalidSolutionViolation(std::vector<MyGraph::Node>& solution) {
//     repairSolutionByCriterion(REMOVAL_CRITERION_LENGTH, solution);
// }
// void Solver::fixInvalidSolutionSparingly(std::vector<MyGraph::Node>& solution) {
//     repairSolutionByCriterion(REMOVAL_CRITERION_VALUE, solution);
// }

void Solver::initializeSampledTable(int minBudget, int maxBudget, int tableSize) {
    /*
     * go from min_budget to max_budget, split the interval with tableSize=K points
     * and initialize with empty solutions
     * 
     * i=0: min
     * i=(K-1): max
     * 
     * for i=0,..,K-1:
     * min + i* (max-min)/(K-1)
     * 
     * It is also necessary to check whether K is larger than (max-min+1)
     */
    if (tableSize > (maxBudget - minBudget + 1)) {
        tableSize = maxBudget - minBudget + 1;
    }
    tableSize_ = tableSize;
    minBudget_ = minBudget;
    maxBudget_ = maxBudget;
    
    std::vector<int> samplePoints;
    samplePoints.push_back(minBudget_);

    for (int i = 0; i < tableSize_; i++) {
        double intervalPoint = minBudget + i*(maxBudget-minBudget) / (tableSize - 1);
        int currentBudget = std::floor(intervalPoint);
        EvaluatedSolution emptySolution;
        bestSolutionsByBudget_[currentBudget] = emptySolution;
    }
}
void Solver::updateTable(std::vector<MyGraph::Node> solution, ResultData resultData) {
    /*
     * groesse checken: sind schon genug Loesungen in der Map drin?
     * dann für gleichlange oder laengere Loesungen checken: ist solution besser?
     * 
     * look at all budget values that are equal or larger than resultData.length
     * and update their corresponding solutions, if necessary
     */
    std::map<int,EvaluatedSolution>::iterator itUpper = bestSolutionsByBudget_.lower_bound(resultData.length_);

    while (itUpper != bestSolutionsByBudget_.end()) {
        if (resultData > itUpper->second.resultData_) {
            itUpper->second.solution_.assign(solution.begin(), solution.end());
            itUpper->second.resultData_ = resultData;
        }
        itUpper++;
    }
}

void Solver::repairSolutionByCriterion(int removeCriterion, std::vector<MyGraph::Node>& mySolution,
                                       ResultData& currentSolutionQuality) {
    
    const int LARGE_NUMBER = 1000000;
    
    // ResultData currentSolutionQuality = evaluateSolution(problemData_, mySolution, "value", false, false, &distanceMatrix_);

    std::vector<int> solutionIndices = nodeVectorToIndexVector(problemData_, mySolution); // a +1 is already in there
    Rcpp::Rcout << "mySol0: " << mySolution.size() << "\n";
    Rcpp::Rcout << "solInd0: " << solutionIndices.size() << "\n";
    printVector(solutionIndices);
    Rcpp::Rcout << currentSolutionQuality << "\n";
    
    double minRemoveValue = 0;
    while (currentSolutionQuality.length_ > problemData_.budget_) {
        // remove node with lowest removeValue (which depends on removeCriterion)
        std::vector<MyGraph::Node> tempSolution;
        minRemoveValue = std::numeric_limits<double>::max();
        std::size_t minRemoveIndex = 0;
        
        Rcpp::Rcout << "mySol: " << mySolution.size() << " - " << minRemoveIndex << "\n";
        Rcpp::Rcout << "solInd: " << solutionIndices.size() << " - " << minRemoveIndex << "\n";
        for (std::size_t i = 0; i < mySolution.size(); i++){
            // Rcpp::Rcout << "fi/fsize: " << i << "/" << mySolution.size() << "-" << solutionIndices.size() << "\n";
            ResultData res;
            additionalLogData_.numberOfTestedBitVectors_++;
            
            double removeValue = 0;
            
            switch(removeCriterion) {
            case(REMOVAL_CRITERION_RATIO):
            {
                // res = evaluateSolutionMatrix(problemData_, tempSolution);
                res = estimatePathLengthDecrease(solutionIndices, i, mySolution, currentSolutionQuality);
                if (res.length_ - currentSolutionQuality.length_ != 0) {
                    removeValue =  problemData_.nodeMap_[mySolution[i]].value_/(currentSolutionQuality.length_ - res.length_);
                } else {
                    removeValue =  LARGE_NUMBER;
                }

            }
                break;
            case(REMOVAL_CRITERION_LENGTH):
                // res = evaluateSolutionMatrix(problemData_, tempSolution);
                res = estimatePathLengthDecrease(solutionIndices, i, mySolution, currentSolutionQuality);
                if (res.length_ - currentSolutionQuality.length_ != 0) {
                    removeValue =  1/(currentSolutionQuality.length_ - res.length_);
                } else {
                    removeValue =  LARGE_NUMBER;
                }
                break;
            case(REMOVAL_CRITERION_VALUE):
                // res = evaluateSolutionMatrix(problemData_, tempSolution);
                res = estimatePathLengthDecrease(solutionIndices, i, mySolution, currentSolutionQuality);
                removeValue =  problemData_.nodeMap_[mySolution[i]].value_;
                break;
            default:
                Rcpp::Rcerr << "invalid remove criterion number " << removeCriterion << ". Use CRITERION_RATIO.\n";
                // res = evaluateSolutionMatrix(problemData_, tempSolution);
                res = estimatePathLengthDecrease(solutionIndices, i, mySolution, currentSolutionQuality);
                if (res.length_ - currentSolutionQuality.length_ != 0) {
                    removeValue =  problemData_.nodeMap_[mySolution[i]].value_/(currentSolutionQuality.length_ - res.length_);
                } else {
                    removeValue =  LARGE_NUMBER;
                }
            break;
            }
            
            if (removeValue < minRemoveValue){
                minRemoveValue = removeValue;
                minRemoveIndex = i;
            }
        }
        
        mySolution.erase(mySolution.begin() + minRemoveIndex);
        solutionIndices.erase(solutionIndices.begin() + minRemoveIndex);
        Rcpp::Rcout << "mySol2: " << mySolution.size() << " - " << minRemoveIndex << "\n";
        Rcpp::Rcout << "solInd2: " << solutionIndices.size() << " - " << minRemoveIndex << "\n\n";
        currentSolutionQuality = evaluateSolutionMatrix(problemData_, mySolution, "value");
        
    } 
    // Rcpp::Rcout << "Handling done: " << currentSolutionQuality << "\n";
}

ResultData Solver::estimatePathLengthDecrease(const std::vector<int>& solutionIndices, int i, 
                                              const std::vector<MyGraph::Node>& tempSolution,
                                              const ResultData& currentSolutionQuality) {
    ResultData res;
    double pathLengthDecrease = 0.0;
    
    if (i == solutionIndices.size() - 1) {
        if (solutionIndices.size() == 1) {
            // edge case
            pathLengthDecrease += getAndLogDistance(solutionIndices[i], 0)
            +getAndLogDistance(0, solutionIndices[i]);
        } else {
            pathLengthDecrease += getAndLogDistance(solutionIndices[i], 0)
            +getAndLogDistance(solutionIndices[i-1], solutionIndices[i])
            -getAndLogDistance(solutionIndices[i-1], 0);
        }
    } else {
        if (i == 0) {
            pathLengthDecrease += getAndLogDistance(0, solutionIndices[0])
            +getAndLogDistance(solutionIndices[0], solutionIndices[1])
            -getAndLogDistance(0, solutionIndices[1]);
        } else {
            // Rcpp::Rcout << "i/size: " << i << "/" << solutionIndices.size() << "\n";
            pathLengthDecrease += getAndLogDistance(solutionIndices[i], solutionIndices[i+1])
            +getAndLogDistance(solutionIndices[i-1], solutionIndices[i])
            -getAndLogDistance(solutionIndices[i-1], solutionIndices[i+1]);
        }
    }
    
    int prev = (i==0) ? 0 : solutionIndices[i-1];
    int curr = solutionIndices[i];
    int next = (i == (int) solutionIndices.size()-1) ? 0 : solutionIndices[i+1];
    
    res.length_ = currentSolutionQuality.length_ - pathLengthDecrease;
    res.value_ = currentSolutionQuality.value_ - problemData_.nodeMap_[tempSolution[i]].value_;
    return res;
}



