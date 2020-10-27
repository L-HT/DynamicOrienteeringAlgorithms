#ifndef DYNAMIC_CHANGE_H
#define DYNAMIC_CHANGE_H

#include <chrono>
#include <iostream>
// #include "ProblemData.h"
// #include "HelperFunctions.h"

// -naher Ort ändert Wert (hoch/runter)
// -mittlerer Ort ändert Wert (hoch/runter)
// -ferner Ort ändert Wert (hoch/runter)
// -Budget sinkt/steigt um 10% (Standard: so, dass man ungefähr die Hälfte bekommen kann)
// -Kante ändert ihre Länge (wird erst mal verschoben, da viel Arbeit (Bogen suchen, Distanzmatrix aktualisieren, etc.))

/*
 * To-Do: implement dynamics changing length of arcs
 * To-Do: add dynamics for atDistanceEvaluation_ (operator><)
 */
struct ProblemData;

struct AdditionalLogData;

struct Change{
    int id_;
    long atEvaluation_;
    long atTestedBitVector_;
    long atDistanceEvaluation_;
    bool isRelativeChange_;
    double atTime_;
    int nodeID1_;
    int nodeID2_;
    double change_;

    bool isNodeChange(){
        return (nodeID1_ != -1 && nodeID2_ == -1);
    }
    bool isArcChange(){
        return (nodeID1_ != -1 && nodeID2_ != -1);
    }
    bool isBudgetChange(){
        return (nodeID1_ == -1 && nodeID2_ == -1);
    }

    bool operator<(const Change& b){
        if (atTime_ == -1 && b.atTime_ == -1 && atTestedBitVector_ == -1 && b.atTestedBitVector_ == -1 && atDistanceEvaluation_ == -1 && b.atDistanceEvaluation_ == -1){
            return atEvaluation_ < b.atEvaluation_;
        }
        if (atEvaluation_ == -1 && b.atEvaluation_ == -1 && atTestedBitVector_ == -1 && b.atTestedBitVector_ == -1 && atDistanceEvaluation_ == -1 && b.atDistanceEvaluation_ == -1){
            return atTime_ < b.atTime_;
        }
        if (atEvaluation_ == -1 && b.atEvaluation_ == -1 && atTime_ == -1 && b.atTime_ == -1 && atDistanceEvaluation_ == -1 && b.atDistanceEvaluation_ == -1){
            return atTestedBitVector_ < b.atTestedBitVector_;
        }
        if (atEvaluation_ == -1 && b.atEvaluation_ == -1 && atTime_ == -1 && b.atTime_ == -1 && atTestedBitVector_ == -1 && b.atTestedBitVector_ == -1){
            return atDistanceEvaluation_ < b.atDistanceEvaluation_;
        }
        return false;
    }
    bool operator>(const Change& b){
        if (atTime_ == -1 && b.atTime_ == -1 && atTestedBitVector_ == -1 && b.atTestedBitVector_ == -1 && atDistanceEvaluation_ == -1 && b.atDistanceEvaluation_ == -1){
            return atEvaluation_ > b.atEvaluation_;
        }
        if (atEvaluation_ == -1 && b.atEvaluation_ == -1 && atTestedBitVector_ == -1 && b.atTestedBitVector_ == -1 && atDistanceEvaluation_ == -1 && b.atDistanceEvaluation_ == -1){
            return atTime_ > b.atTime_;
        }

        if (atEvaluation_ == -1 && b.atEvaluation_ == -1 && atTime_ == -1 && b.atTime_ == -1 && atDistanceEvaluation_ == -1 && b.atDistanceEvaluation_ == -1){
            return atTestedBitVector_ > b.atTestedBitVector_;
        }

        if (atEvaluation_ == -1 && b.atEvaluation_ == -1 && atTime_ == -1 && b.atTime_ == -1 && atTestedBitVector_ == -1 && b.atTestedBitVector_ == -1){
            return atDistanceEvaluation_ > b.atDistanceEvaluation_;
        }

        return true; // by default?? Wait, why did I do this...?
    }

    friend std::ostream& operator<<(std::ostream& os, const Change& change){
        os << "id: " << change.id_ << ", time: " << change.atTime_  << ", eval: " << change.atEvaluation_ << ", bitVector: " << change.atTestedBitVector_ << ", ";
        os << "distEval: " << change.atDistanceEvaluation_ << ", ""relative: " << (change.isRelativeChange_ ? "relative" : "absolute") << "; id1: " << change.nodeID1_;
        os << "; id2: " << change.nodeID2_ << "; change: " << change.change_;
        return os;
    }
};

bool checkForAndImplementChanges(ProblemData& problemData, AdditionalLogData additionalLogData, std::chrono::high_resolution_clock::time_point elapsedTime);

void implementChange(ProblemData& problemData, Change& change);


#endif
