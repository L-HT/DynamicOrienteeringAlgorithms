// #include <chrono>

#include "DynamicChange.h"
#include "ProblemData.h"
#include "HelperFunctions.h"
#include "IDSearchFunctions.h"


bool checkForAndImplementChanges(ProblemData& problemData, AdditionalLogData additionalLogData, std::chrono::high_resolution_clock::time_point startTime){

    std::chrono::high_resolution_clock::time_point currentTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> timeSpan = std::chrono::duration_cast<std::chrono::duration<double>>(currentTime - startTime);

    double timePassed = timeSpan.count();

    bool changeHappened = false;
    // check changes occuring by time
    auto it = problemData.dynamicChangesByTime_.begin();
    while (it != problemData.dynamicChangesByTime_.end() && (*it).atTime_ <= timePassed) {
        Change& change = *it;
        Rcpp::Rcout << "change by time: " << change << std::endl;
        implementChange(problemData, change);
        it = problemData.dynamicChangesByTime_.erase(it);
        changeHappened = true;
    }

    // check changes occuring by evaluation
    it = problemData.dynamicChangesByEvaluation_.begin();
    while (it != problemData.dynamicChangesByEvaluation_.end() && (*it).atEvaluation_ <= additionalLogData.numberOfCalculatedPaths_) {
        Change& change = *it;
        Rcpp::Rcout << "change by evaluation: " << change << std::endl;
        implementChange(problemData, change);
        it = problemData.dynamicChangesByEvaluation_.erase(it);
        changeHappened = true;
    }

    // check changes occuring by tested bit vectors
    it = problemData.dynamicChangesByTestedBitVectors_.begin();
    while (it != problemData.dynamicChangesByTestedBitVectors_.end() && (*it).atTestedBitVector_ <= additionalLogData.numberOfTestedBitVectors_) {
        Change& change = *it;
        Rcpp::Rcout << "change by tested bit vector: " << change << std::endl;
        implementChange(problemData, change);
        it = problemData.dynamicChangesByTestedBitVectors_.erase(it);
        changeHappened = true;
    }

    // check changes occuring by distance evaluations
    it = problemData.dynamicChangesByDistanceEvaluation_.begin();
    while (it != problemData.dynamicChangesByDistanceEvaluation_.end() && (*it).atDistanceEvaluation_ <= additionalLogData.numberOfShortestPathCalls_) {
        Change& change = *it;
        Rcpp::Rcout << "change by distance evaluation: " << change << std::endl;
        implementChange(problemData, change);
        it = problemData.dynamicChangesByDistanceEvaluation_.erase(it);
        changeHappened = true;
    }

    // // check mixed changes
    // it = problemData.dynamicChangesMixed_.begin();
    // while (it != problemData.dynamicChangesMixed_.end()) {
    //     Change& change = *it;
    //     Rcpp::Rcout << "mixed change: " << change << std::endl;
    //     if (change.atEvaluation_ <= additionalLogData.numberOfCalculatedPaths_ || change.atTime_ <= timePassed || change.atTestedBitVector_ <= additionalLogData.numberOfTestedBitVectors_) {
    //         implementChange(problemData, change);
    //         it = problemData.dynamicChangesMixed_.erase(it);
    //         changeHappened = true;
    //     } else {
    //         ++it;
    //     }
    // }
    // waitForInput("wait44!", true);
    return changeHappened;
}

void implementChange(ProblemData& problemData, Change& change) {

    if (change.isArcChange()) {
        /*
         * Not implemented yet
         * To-Do: If there is a change in arc length, does the distance matrix have to be recalculated?
         * That does not sound very good...
         */

    }
    if (change.isBudgetChange()) {
        double oldBudget = problemData.budget_;
        if (change.isRelativeChange_){
            problemData.budget_ *= change.change_;
        } else {
            problemData.budget_ += change.change_;
        }
        if (problemData.budget_ <= 0){
            Rcpp::Rcerr << "WARNING: Budget change not not applied: " << oldBudget << " -> " << problemData.budget_ << std::endl;
            problemData.budget_ = oldBudget;
        } else {
            Rcpp::Rcout << "change budget: " << oldBudget << " -> " << problemData.budget_ << std::endl;
        }
    }
    if (change.isNodeChange()) {


        MyGraph::Node n = getNodeFromInternalID(problemData.graph_, problemData.nodeMap_, change.nodeID1_);
        double oldValue = problemData.nodeMap_[n].value_;
        if (change.isRelativeChange_){
            problemData.nodeMap_[n].value_ *= change.change_;
        } else {
            problemData.nodeMap_[n].value_ += change.change_;
        }
        if (problemData.nodeMap_[n].value_ < 0){
            Rcpp::Rcerr << "WARNING -- value change of node id " << problemData.nodeMap_[n].id_ << " not applied: " << oldValue << " -> " << problemData.nodeMap_[n].value_ << std::endl;
            problemData.nodeMap_[n].value_ = oldValue;
        } else {
            Rcpp::Rcout << "change node value of id " << problemData.nodeMap_[n].id_ << ": " << oldValue << " -> " << problemData.nodeMap_[n].value_ << std::endl;
        }
    }

}
