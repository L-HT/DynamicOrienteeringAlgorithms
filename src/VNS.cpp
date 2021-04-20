#include <Rcpp.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
#include <algorithm>
#include <random>
#include <cmath>

#include "Solver.h"
#include "PlotRoute.h"
#include "IO.h"
#include "ProblemData.h"
#include "CalculateLength.h"
#include "HelperFunctions.h"
#include "IDSearchFunctions.h"
#include "DistanceMatrix.h"
#include "linkern_wrapper.h"

#define CRITERION_RATIO 1
#define CRITERION_LENGTH 2
#define CRITERION_VALUE 3
#define CRITERION_RANDOM 4

#define LARGE_NUMBER 1000000

using namespace lemon;

struct VariableNeighborhoodSearch : public Solver{
    std::vector<MyGraph::Node> baseSolution_;
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> randomNumber; //Intervall [0,1)

    ResultData currentSolutionQuality_;

    int iterationsWithForcedCriterionRatio_;
    double perturbationProbability_;

    int initialRepetitions_;
    double probabilityForRandomInsertion_;

    bool withLocalSearch_;
    bool withLinKernighan_;

    VariableNeighborhoodSearch(ProblemData& problemData,
              std::string logFileName, unsigned int runNumber, Rcpp::Environment rEnvironment,
              std::string targetCriterion, std::string fileSuffix, std::string pathToDistanceMatrix)

        : Solver(problemData, logFileName, runNumber, rEnvironment,
          targetCriterion, "vns", fileSuffix, pathToDistanceMatrix) {

        gen = std::mt19937(rd());
        randomNumber = std::uniform_real_distribution<double>(0.0, 1.0);
        iterationsWithForcedCriterionRatio_ = (int) std::sqrt(problemData_.destinations_.size());
        perturbationProbability_ = 0; // oldValue: 0.15 // but now perturbation is done by insert/removeCriterion=4
        initialRepetitions_ = 1;//(int) std::sqrt(problemData_.destinations_.size());
    }

    void run(){
        Rcpp::Rcout << "start " << algorithmName_ << " (as improver: " << calledAsImprover_ << ")\n";
        std::vector<MyGraph::Node> mySolution;
        bool possibleToAddNode = true;
        additionalLogData_.currentPhase_ = 0;

        /*
        -tourImprov verbessert die Tourlänge
            -drop entfernt den Knoten, der die meiste Entfernung verbraucht
            -add fügt Knoten hinzu, die relativ zu ihrem (kleinstmöglichen) Weglängen-Zuwachs viel Wert bringen und weiterhin gültig sind
            -dieser Weglängenzuwachs wird aber nur geschätzt für die 3 nächsten Knoten
            -Idee: Ermittlung der 3 nächsten Knoten braucht schon O(n_sol) Distanzauswertungen.
            -gehe einfach nur greedy vor -> kann das zu lok. Optima führen? laut Paper scheint das langasam zu sein
            -es kann sein, dass vor drop ungültige Lösungen vorliegen, aber man trotzdem mit denen weiterarbeitet
            -> drop entfernt solange (nach Kriterium val/dist), bis die Lösung wieder gültig wird

         Jetzt hier neue Algorithmus:
            -nimm erst mal LK zur Pfad-Verbesserung, weil sie sehr wenige Pfad-Evals und meist nur Dist-Evals macht.
                -erst, wenn ein Durchlaufen mit LK keine Verbesserung bringt, kann man LS probieren
            -versuche am Anfang, möglichst gierig Knoten hinzuzufügen und sie zu testen
                -gierig bezüglich verschiedener Maße: nur Pfadverlängerung, nur Wert, nur Wert/Länge

                -Pfadverlängerung: füge den Knoten dort ein, wo der Pfad am wenigsten lang wird

            -bei BitVector ist aber AddVertex von GRASP-SR auch gut: Also auch mal dort reinschauen...
                -Idee dort ist: Für jeden unbesuchten Knoten: Füge ihn an der Stelle mit geringster Wegverlängerung ein
                    -Wenn über Budget, dann entferne von links ausgehend wachsende Segmente: 1, 12, 123,... 2, 23, 234,..., 3, 34,...
                        -solange, bis eine gültige Lösung entsteht
                    -alle gültigen Lösungen kommen in die Kandidatenliste rein

            -Idee hier: Speichere für jeden Knoten in der Lösung:
                -Distanz zu Nachbarknoten (Damit bei Entfernung dieses Knotens der Wegkürzung schnell ermittelt werden kann)

                -Speichere für jeden Knoten K außerhalb der Lösung:
                    -eine HashMap (KnotenID, Distanz), die sagt, wie lang der Weg von einem Knoten (in der Lösung) zu K ist
                    -irgendwie, um distanceEvaluations einzusparen... To-Do: vielleicht später mal...

            -also:
                -füge erst mal Knoten hinzu, die val/distIncrease maximieren
                -sobald Budget ausgeschöpft: verbessere mit LK
                -wenn das nicht hilft: verbessere mit LS
                -wenn das nicht hilft: entferne Knoten so lange, bis die Lösung wieder gültig wird
                -wiederhole, bis Abbruchkriterium erfüllt

            -immer wechseln

         */

        std::string operationSequence_ = "213";

        if (initialSolution_.empty() && calledAsImprover_){
            Rcpp::stop("Initial solution is missing");
        }

        bestSolution_.assign(initialSolution_.begin(), initialSolution_.end());
        bestSolutionQuality_ = evaluateSolutionMatrix(problemData_, bestSolution_, "value", true);
        Rcpp::Rcout << "bestSolutionQuality: " << bestSolutionQuality_ << std::endl;

        mySolution.assign(bestSolution_.begin(), bestSolution_.end());
        currentSolutionQuality_ = bestSolutionQuality_;

        ///////////////////////////////////////////////////////
        // Initialization (more like initial exploration)
        ///////////////////////////////////////////////////////
//
        additionalLogData_.currentPhase_ = 0;

        std::vector<MyGraph::Node> allNodes(problemData_.destinations_.begin(), problemData_.destinations_.end());
        callLinKernighanHeuristic(allNodes);
        ResultData allNodesResultData = evaluateSolutionMatrix(problemData_, allNodes);
        probabilityForRandomInsertion_ = std::sqrt(problemData_.budget_ / allNodesResultData.length_);
        Rcpp::Rcout << "probabilityForRandomInsertion: " << probabilityForRandomInsertion_ << "\n";

        // double p_estimate = 1-estimateInsertionProbability();
        // probabilityForRandomInsertion_ = 1-estimateInsertionProbability();
        // Rcpp::Rcout << "probabilityForRandomInsertion2: " << probabilityForRandomInsertion_ << "\n";
//         Rcpp::Rcout << "init probabilityForRandomInsertion: " << p_estimate << "\n";
//         for (int i = 0; i < initialRepetitions_; i++) {
//             /*
//              * Ein Intervall [p, 1-p] mit p \in (0,1) soll in n Intervalle aufgeteilt werden, sodass
//              * x_1 = p, x_n = 1-p und der Rest äquidistant verteilt ist.
//              * Antwort: x_i = p + (i-1) * (1-2p)/n
//              */
//             probabilityForRandomInsertion_ = p_estimate + i * (1.0-2*p_estimate) / initialRepetitions_;
//             std::vector<MyGraph::Node> tempSolution(mySolution.begin(), mySolution.end());
//             Rcpp::Rcout << "probabilityForRandomInsertion1: " << probabilityForRandomInsertion_ << "\n";
//             randomInsert(tempSolution, probabilityForRandomInsertion_);
//             additionalLogData_.numberOfTestedBitVectors_++;
//             callLinKernighanHeuristic(tempSolution);
//             // if better solution is found, it is automatically written to bestSolution
//             evaluateSolutionMatrix(problemData_, tempSolution, "value", true);
//
//         }
//         mySolution.assign(bestSolution_.begin(), bestSolution_.end());
//         evaluateSolutionMatrix(problemData_, bestSolution_, "value", true);
//         currentSolutionQuality_ = bestSolutionQuality_;

        Rcpp::Rcout << "Start main loop...\n";
        do {

            ///////////////////
            // adding nodes...
            ///////////////////
            additionalLogData_.currentPhase_ = 1;
            int insertCriterion = getRandomNumber(1,4);
            // if (additionalLogData_.currentIteration_ < 10) {
            //     insertCriterion = CRITERION_RATIO;
            //     greedyInsert(insertCriterion, mySolution);
            // } else {
            if (randomNumber(gen) < perturbationProbability_
                    || additionalLogData_.currentIteration_ < iterationsWithForcedCriterionRatio_
                    || mySolution.empty()) {
                insertCriterion = 4;
                randomInsert(mySolution, probabilityForRandomInsertion_);
            } else {
                greedyInsert(insertCriterion, mySolution);
            }
            // }
            Rcpp::checkUserInterrupt();

            ////////////////////////////////
            // call Lin-Kernighan-Heuristic
            ////////////////////////////////
            if (withLinKernighan_) {
                additionalLogData_.currentPhase_ = 2;
                callLinKernighanHeuristic(mySolution);
                currentSolutionQuality_ = evaluateSolutionMatrix(problemData_, mySolution, "value");
            }
            Rcpp::checkUserInterrupt();

            ////////////////
            // Local Search
            ////////////////
            if (currentSolutionQuality_.length_ > problemData_.budget_ &&  withLocalSearch_) {
                additionalLogData_.currentPhase_ = 3;
                baseSolution_.assign(mySolution.begin(), mySolution.end());
                std::shuffle(baseSolution_.begin(), baseSolution_.end(), gen);

                //////////// Lokale Suche //////////////
                bool lsSuccessful = false;
                // waitForInput("startLs", DEBUG_ENABLED);
                for (char c : operationSequence_){
                    //waitForInput("ls running...", DEBUG_ENABLED);
                    int op = c - '0';
                    switch(op){
                    case(OP_INSERT):
                        additionalLogData_.currentPhase_ = 31;
                        do{
                            // Rcpp::Rcout << "LS: Insert " << improvementFoundInCurrentIteration_ << std::endl;
                            lsSuccessful = localSearchByReferencePermutation(mySolution, "insert");
                        } while (lsSuccessful);
                        break;
                    case(OP_EDGEINSERT):
                        additionalLogData_.currentPhase_ = 32;
                        do{
                            // Rcpp::Rcout << "LS: EdgeInsert " << improvementFoundInCurrentIteration_ << std::endl;
                            lsSuccessful = localSearchByReferencePermutation(mySolution, "edgeInsert");
                        } while (lsSuccessful);
                        // Rcpp::Rcout << "LS: EdgeInsert " << improvementFoundInCurrentIteration_ << std::endl;
                        break;
                    case(OP_SWAP):
                        do{
                            additionalLogData_.currentPhase_ = 33;
                            // Rcpp::Rcout << "LS: Swap " << improvementFoundInCurrentIteration_ << std::endl;
                            lsSuccessful = localSearchByReferencePermutation(mySolution, "swap");
                        } while (lsSuccessful);
                        // Rcpp::Rcout << "LS: Swap " << improvementFoundInCurrentIteration_ << std::endl;
                        break;
                    default:
                        Rcpp::stop("op ist unbekannt");
                    break;

                    }
                    Rcpp::checkUserInterrupt();
                }
                currentSolutionQuality_ = evaluateSolutionMatrix(problemData_, mySolution, "value");
            }

            ////////////////////////////////
            // remove nodes
            ////////////////////////////////

            additionalLogData_.currentPhase_ = 4;
            int removeCriterion = getRandomNumber(1,4);

            // initial removal by CRITERION_RATIO + Logging
            if (additionalLogData_.currentIteration_ < iterationsWithForcedCriterionRatio_) {
                removeCriterion = 1;
                evaluateSolutionMatrix(problemData_, mySolution, "value", true);
            } //else {
                // if (randomNumber(gen) < perturbationProbability_|| additionalLogData_.currentIteration_ < 3) {
                //     removeCriterion = 4;
                // }
            // }
            if (currentSolutionQuality_.length_ > problemData_.budget_){
                greedyRemove(removeCriterion, mySolution);
            }
            // Rcpp::stop("myEnd");

            additionalLogData_.currentIteration_++;
            Rcpp::checkUserInterrupt();
            // Rcpp::Rcout << insertCriterion <<  removeCriterion << " ";
            //waitForInput("iteration ende", DEBUG_ENABLED);
        } while (!terminationCriterionSatisfied());

        evaluateSolutionMatrix(problemData_, bestSolution_, "value", true);
        Rcpp::Rcout << "end " << algorithmName_ << std::endl;
    }

    // randomly insert nodes (very exploration-focused procedure at the beginning)
    /*
     * n Distancepolls -> Schätzung einer zufälligen Route durch n Knoten = D~
     *
     * Wenn budget/D~ >= 1: Es kann gut sein, dass alle Knoten rein passen
     * Wenn budget/D~ < 1: Ein paar Knoten müssen raus.
     *
     * Laut EA-Paper: p=budget/D~ wäre aber eine Unterschätzung.
     * Wähle von daher sqrt(budget/D~)
     *
     */
    void randomInsert(std::vector<MyGraph::Node>& mySolution, double p) {
        for (MyGraph::Node n : problemData_.destinations_){
            if (std::find(mySolution.begin(), mySolution.end(), n) == mySolution.end()){
                if (randomNumber(gen) < p) {
                    std::size_t numberOfPossibleInsertionPositions = mySolution.size() + 1;
                    std::size_t k = getRandomNumber(0, numberOfPossibleInsertionPositions - 1);
                    if (k == numberOfPossibleInsertionPositions - 1){
                        mySolution.push_back(n);
                    } else {
                        mySolution.insert(mySolution.begin() + k, n);
                    }
                }
            }
        }
    }

    double estimateInsertionProbability() {
        std::size_t n = problemData_.destinations_.size();
        double distanceSum = 0.0;
        for (std::size_t i = 0; i < n; i++) {
            std::size_t randomRow = getRandomNumber(0, distanceMatrix_.nrow()-2);
            std::size_t randomCol = getRandomNumber(randomRow+1, distanceMatrix_.nrow()-1);
            distanceSum += distanceMatrix_(randomRow, randomCol);
        }
        double result = (double) problemData_.budget_ / distanceSum;
        return std::sqrt(result);
    }
    void greedyInsert(int insertCriterion, std::vector<MyGraph::Node>& mySolution) {
        insertCriterion = CRITERION_LENGTH;
        bool possibleToAddNode = false;
        do{
            // waitForInput("adding node begin", DEBUG_ENABLED);
            possibleToAddNode = false;
            std::vector<MyGraph::Node> tempSolution;//(mySolution.begin(), mySolution.end());
            MyGraph::Node bestNodeToInsert = lemon::INVALID;
            std::size_t bestInsertionPosition;
            double currentGoodness = 0;
            double maxGoodness = 0;
            std::size_t numberOfPossibleInsertionPositions = mySolution.size() + 1;

            std::vector<int> solutionIndices = nodeVectorToIndexVector(problemData_, mySolution); // a +1 is already in there

            for (std::size_t destinationIndex = 0; destinationIndex < problemData_.destinations_.size(); destinationIndex++){
                MyGraph::Node n = problemData_.destinations_[destinationIndex];
                if (std::find(mySolution.begin(), mySolution.end(), n) == mySolution.end()){

                    // test all possible insertion positions
                    additionalLogData_.numberOfTestedBitVectors_++;
                    double tempBestGoodness = 0;
                    std::size_t tempBestPosition = 0;

                    for (std::size_t i = 0; i < numberOfPossibleInsertionPositions; i++){
                        // tempSolution.assign(mySolution.begin(), mySolution.end());
                        //
                        // if (i == numberOfPossibleInsertionPositions - 1){
                        //     tempSolution.push_back(n);
                        // } else {
                        //     tempSolution.insert(tempSolution.begin() + i, n);
                        // }

                        // printNodeIdsOfVector(problemData_, tempSolution);


                        /*
                         * estimate increase in path length instead of fully evaluating a solution
                         */
                        ResultData res;
                        double pathLengthDifference = 0.0;
                        // if ((numberOfPossibleInsertionPositions - i < 3 || i < 3) && i > 0 && i < numberOfPossibleInsertionPositions) {
                        //     Rcpp::Rcout << "i/number: " << i << "/" << numberOfPossibleInsertionPositions << ", destIndex: " << destinationIndex << "\n";
                        //     Rcpp::Rcout << "solutionIndices[i-1]:  " << solutionIndices[i-1] << ", solutionIndices[i]: " << solutionIndices[i] << "\n";
                        // }
                        // if (i == numberOfPossibleInsertionPositions - 1) {
                        //     pathLengthDifference += getAndLogDistance(solutionIndices[i-1], destinationIndex+1)
                        //                             +getAndLogDistance(destinationIndex+1, 0)
                        //                             -getAndLogDistance(solutionIndices[i-1], 0);
                        // } else {
                        //     if (i == 0) {
                        //         pathLengthDifference += getAndLogDistance(destinationIndex+1, solutionIndices[i])
                        //         +getAndLogDistance(0, destinationIndex+1)
                        //         -getAndLogDistance(0, solutionIndices[i]);
                        //     } else {
                        //         pathLengthDifference += getAndLogDistance(solutionIndices[i-1], destinationIndex+1)
                        //                                 +getAndLogDistance(destinationIndex+1, solutionIndices[i])
                        //                                 -getAndLogDistance(solutionIndices[i-1], solutionIndices[i]);
                        //     }
                        // }
                        // res.length_ = currentSolutionQuality_.length_ + pathLengthDifference;
                        // res.value_ = currentSolutionQuality_.value_ + problemData_.nodeMap_[n].value_;

                        switch(insertCriterion) {
                        case(CRITERION_RATIO):

                            // res = evaluateSolutionMatrix(problemData_, tempSolution);
                            res = estimatePathLengthIncrease(solutionIndices, i, n, numberOfPossibleInsertionPositions, destinationIndex);
                            if (res.length_ - currentSolutionQuality_.length_ != 0) {
                                currentGoodness =  (res.value_ - currentSolutionQuality_.value_) / (res.length_ - currentSolutionQuality_.length_);
                            } else {
                                currentGoodness =  LARGE_NUMBER;
                            }
                            break;
                        case(CRITERION_LENGTH):
                            // res = evaluateSolutionMatrix(problemData_, tempSolution);
                            res = estimatePathLengthIncrease(solutionIndices, i, n, numberOfPossibleInsertionPositions, destinationIndex);
                            if (res.length_ - currentSolutionQuality_.length_ != 0) {
                                currentGoodness =  1/(res.length_ - currentSolutionQuality_.length_);
                            } else {
                                currentGoodness =  LARGE_NUMBER;
                            }
                            break;
                        case(CRITERION_VALUE):
                            // res = evaluateSolutionMatrix(problemData_, tempSolution);
                            res = estimatePathLengthIncrease(solutionIndices, i, n, numberOfPossibleInsertionPositions, destinationIndex);
                            currentGoodness =  (res.value_ - currentSolutionQuality_.value_);
                            break;
                        case(CRITERION_RANDOM):
                            currentGoodness = randomNumber(gen); // uniformly distributed between 0 and 1
                            break;
                        default:
                            Rcpp::Rcerr << "invalid insert criterion number " << insertCriterion << ". Use CRITERION_RATIO.\n";
                            // res = evaluateSolutionMatrix(problemData_, tempSolution);
                            res = estimatePathLengthIncrease(solutionIndices, i, n, numberOfPossibleInsertionPositions, destinationIndex);
                            if (res.length_ - currentSolutionQuality_.length_ != 0) {
                                currentGoodness =  (res.value_ - currentSolutionQuality_.value_)/(res.length_ - currentSolutionQuality_.length_);
                            } else {
                                currentGoodness =  LARGE_NUMBER;
                            }
                            break;
                        }
                        if (currentGoodness > tempBestGoodness){
                            tempBestPosition = i;
                            tempBestGoodness = currentGoodness;
                        }
                    }

                    if (tempBestGoodness > maxGoodness){
                        bestNodeToInsert = n;
                        maxGoodness = tempBestGoodness;
                        bestInsertionPosition = tempBestPosition;
                    }
                }
            }

            if (bestNodeToInsert != lemon::INVALID) {
                // solution might exceed budget at this point, but insert anyway
                mySolution.insert(mySolution.begin() + bestInsertionPosition, bestNodeToInsert);
                currentSolutionQuality_ = evaluateSolutionMatrix(problemData_, mySolution, "value");

                possibleToAddNode = currentSolutionQuality_.length_ <= problemData_.budget_;
            } else {
                // no node can be added since there are either no unvisited nodes or because
                // these nodes somehow make the solution worse
                Rcpp::Rcerr << "No node can be inserted. Abort greedyInsert...\n";
                possibleToAddNode = false;
            }
        } while (possibleToAddNode);
    }
    ResultData estimatePathLengthIncrease(const std::vector<int>& solutionIndices, int i,
                                          const MyGraph::Node n, std::size_t numberOfPossibleInsertionPositions,
                                          int destinationIndex) {
        ResultData res;
        double pathLengthDifference = 0.0;

        if (i == numberOfPossibleInsertionPositions - 1) {
            if (i == 0) { //special case: empty solution
                // Rcpp::Rcout << "Info: Empty solution has been given to greedyInsert. \n";
                pathLengthDifference += getAndLogDistance(0, destinationIndex+1)+getAndLogDistance(destinationIndex+1, 0);
            } else {
                pathLengthDifference += getAndLogDistance(solutionIndices[i-1], destinationIndex+1)
                +getAndLogDistance(destinationIndex+1, 0)
                -getAndLogDistance(solutionIndices[i-1], 0);
            }
        } else {
            if (i == 0) {
                pathLengthDifference += getAndLogDistance(destinationIndex+1, solutionIndices[i])
                +getAndLogDistance(0, destinationIndex+1)
                -getAndLogDistance(0, solutionIndices[i]);
            } else {
                pathLengthDifference += getAndLogDistance(solutionIndices[i-1], destinationIndex+1)
                +getAndLogDistance(destinationIndex+1, solutionIndices[i])
                -getAndLogDistance(solutionIndices[i-1], solutionIndices[i]);
            }
        }
        res.length_ = currentSolutionQuality_.length_ + pathLengthDifference;
        res.value_ = currentSolutionQuality_.value_ + problemData_.nodeMap_[n].value_;
        return res;
    }

    void greedyRemove(int removeCriterion, std::vector<MyGraph::Node>& mySolution) {

        if (mySolution.empty()){
            Rcpp::stop("mySolution is empty. Does a feasible solution even exist?");
        }
        std::vector<int> solutionIndices = nodeVectorToIndexVector(problemData_, mySolution); // a +1 is already in there

        double minRemoveValue = 0;
        do {
            // remove node with lowest removeValue (which depends on removeCriterion)
            std::vector<MyGraph::Node> tempSolution;
            minRemoveValue = std::numeric_limits<double>::max();
            std::size_t minRemoveIndex = 0;
            for (std::size_t i = 0; i < mySolution.size(); i++){
                // tempSolution.assign(mySolution.begin(), mySolution.end());
                // tempSolution.erase(tempSolution.begin() + i);

                /*
                 * estimate decrease in path length instead of fully evaluating a solution
                 */
                // if (mySolution.size()  - i < 3 || i < 3) {
                //     Rcpp::Rcout << "i/number: " << i << "/" << mySolution.size() << "\n";
                //     Rcpp::Rcout << "solutionIndices[i]:  " << solutionIndices[i] << ", solutionIndices[i+1]: " << solutionIndices[i+1] << "\n";
                // }
                ResultData res;

                // tempSolution.erase(tempSolution.begin() + i);
                additionalLogData_.numberOfTestedBitVectors_++;

                double removeValue = 0;

                switch(removeCriterion) {
                case(CRITERION_RATIO):
                    {
                        // res = evaluateSolutionMatrix(problemData_, tempSolution);
                        res = estimatePathLengthDecrease(solutionIndices, i, mySolution);
                        // { // for Debug purposes -- should be deleted
                        //     std::vector<MyGraph::Node> copySolution(mySolution.begin(), mySolution.end());
                        //     copySolution.erase(copySolution.begin() + i);
                        //     ResultData copyRes = evaluateSolutionMatrix(problemData_, copySolution);
                        //
                        //     double vonHand = currentSolutionQuality_.length_ - copyRes.length_;
                        //     double vonFormel =  currentSolutionQuality_.length_ - res.length_;
                        //     if (vonHand != vonFormel) {
                        //         Rcpp::Rcout << "i: " << i << "/" << mySolution.size() << " - " << vonHand << " / " << vonFormel << "\n";
                        //
                        //         // printNodeIdsOfVector(problemData_, mySolution);
                        //         // printNodeIdsOfVector(problemData_, copySolution);
                        //
                        //         Rcpp::Rcout << "old / new: " << currentSolutionQuality_ << " / " << copyRes << "\n";
                        //         Rcpp::Rcout << "Formel: " << getAndLogDistance(solutionIndices[i], solutionIndices[i+1]) << " + " <<
                        //             getAndLogDistance(solutionIndices[i-1], solutionIndices[i]) << " - " << getAndLogDistance(solutionIndices[i-1], solutionIndices[i+1]) << "\n";
                        //         Rcpp::stop("wtd");
                        //     }
                        // }
                        if (res.length_ - currentSolutionQuality_.length_ != 0) {
                            removeValue =  problemData_.nodeMap_[mySolution[i]].value_/(currentSolutionQuality_.length_ - res.length_);
                        } else {
                            removeValue =  LARGE_NUMBER;
                        }

                    }
                    break;
                case(CRITERION_LENGTH):
                    // res = evaluateSolutionMatrix(problemData_, tempSolution);
                    res = estimatePathLengthDecrease(solutionIndices, i, mySolution);
                    if (res.length_ - currentSolutionQuality_.length_ != 0) {
                        removeValue =  1/(currentSolutionQuality_.length_ - res.length_);
                    } else {
                        removeValue =  LARGE_NUMBER;
                    }
                    break;
                case(CRITERION_VALUE):
                    // res = evaluateSolutionMatrix(problemData_, tempSolution);
                    res = estimatePathLengthDecrease(solutionIndices, i, mySolution);
                    removeValue =  problemData_.nodeMap_[mySolution[i]].value_;
                    break;
                case(CRITERION_RANDOM):
                    removeValue =  randomNumber(gen); //uniformly distributed between 0 and 1
                    break;
                default:
                    Rcpp::Rcerr << "invalid remove criterion number " << removeCriterion << ". Use CRITERION_RATIO.\n";
                    // res = evaluateSolutionMatrix(problemData_, tempSolution);
                    res = estimatePathLengthDecrease(solutionIndices, i, mySolution);
                    if (res.length_ - currentSolutionQuality_.length_ != 0) {
                        removeValue =  problemData_.nodeMap_[mySolution[i]].value_/(currentSolutionQuality_.length_ - res.length_);
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
            // Rcpp::Rcout << minLength << " > " << problemData_.budget_ << std::endl;
            // Rcpp::Rcout << "remove1: " << problemData_.nodeMap_[mySolution[minRemoveIndex]].id_;
            mySolution.erase(mySolution.begin() + minRemoveIndex);
            solutionIndices.erase(solutionIndices.begin() + minRemoveIndex);
            currentSolutionQuality_ = evaluateSolutionMatrix(problemData_, mySolution, "value");
            // Rcpp::Rcout << " -> " << currentSolutionQuality_ << "\n";

        } while (currentSolutionQuality_.length_ > problemData_.budget_);
    }

    ResultData estimatePathLengthDecrease(const std::vector<int>& solutionIndices, int i, const std::vector<MyGraph::Node>& tempSolution) {
        ResultData res;
        double pathLengthDecrease = 0.0;

        if (i == solutionIndices.size() - 1) {
            pathLengthDecrease += getAndLogDistance(solutionIndices[i], 0)
            +getAndLogDistance(solutionIndices[i-1], solutionIndices[i])
            -getAndLogDistance(solutionIndices[i-1], 0);
        } else {
            if (i == 0) {
                pathLengthDecrease += getAndLogDistance(0, solutionIndices[0])
                +getAndLogDistance(solutionIndices[0], solutionIndices[1])
                -getAndLogDistance(0, solutionIndices[1]);
            } else {
                pathLengthDecrease += getAndLogDistance(solutionIndices[i], solutionIndices[i+1])
                +getAndLogDistance(solutionIndices[i-1], solutionIndices[i])
                -getAndLogDistance(solutionIndices[i-1], solutionIndices[i+1]);
            }
        }

        int prev = (i==0) ? 0 : solutionIndices[i-1];
        int curr = solutionIndices[i];
        int next = (i == (int) solutionIndices.size()-1) ? 0 : solutionIndices[i+1];

        res.length_ = currentSolutionQuality_.length_ - pathLengthDecrease;
        res.value_ = currentSolutionQuality_.value_ - problemData_.nodeMap_[tempSolution[i]].value_;
        return res;
    }

    bool localSearchByReferencePermutation(std::vector<MyGraph::Node>& permutation, std::string opString){
        std::size_t counter = 0;
        // std::size_t fallbackCounter = 0;
        double bestValue = std::numeric_limits<double>::max();
        double tempValue = std::numeric_limits<double>::max();
        std::size_t upperLimitForLoop = 0;
        //std::size_t bestTakePosition = 0;
        std::size_t bestOpPosition = 1;

        bool permutationHasImproved = false;

        int op = -1;

        if (opString == "insert"){
            op = OP_INSERT;
        }
        if (opString == "edgeInsert"){
            op = OP_EDGEINSERT;
        }
        if (opString == "swap"){
            op = OP_SWAP;
        }

        std::vector<MyGraph::Node> myCopy;
        ResultData resultData = evaluateSolutionMatrix(problemData_, permutation);

        double currentValue = resultData.length_;

        std::size_t numberOfNodes_ = permutation.size();

        for (counter = 0; counter < numberOfNodes_; counter++){
            std::size_t takePosition = 0;

            for (std::size_t i = 0; i < numberOfNodes_; i++){
                if (permutation[i] == baseSolution_[counter]){
                    takePosition = i;
                }
            }

            if (op == OP_EDGEINSERT && takePosition == numberOfNodes_ - 1){
                continue;
            }

            bestValue = std::numeric_limits<int>::max();

            if (op == OP_EDGEINSERT){
                upperLimitForLoop = numberOfNodes_ - 1;
            } else {
                upperLimitForLoop = numberOfNodes_;
            }
            for (std::size_t opPosition = 0; opPosition < upperLimitForLoop; opPosition++){

                if (takePosition != opPosition){
                    myCopy.assign(permutation.begin(), permutation.end());

                    switch(op){
                    case(OP_INSERT):

                        // Rcpp::Rcout << "Insert " << takePosition << ", " << opPosition << std::endl;
                        insertOperator(myCopy, takePosition, opPosition);
                        break;
                    case(OP_EDGEINSERT):
                        // Rcpp::Rcout << "EdgeInsert " << takePosition << ", " << opPosition << std::endl;
                        edgeInsertOperator(myCopy, takePosition, opPosition);
                        break;
                    case(OP_SWAP):
                        // Rcpp::Rcout << "Swap " << takePosition << ", " << opPosition << std::endl;
                        swapOperator(myCopy, takePosition, opPosition);
                        break;
                    default:
                        Rcpp::stop("op ist unbekannt");
                    break;
                    }
                    ResultData resultData = evaluateSolutionMatrix(problemData_, myCopy);
                    tempValue = resultData.length_;
                    if (tempValue < bestValue){
                        bestValue = tempValue;
                        bestOpPosition = opPosition;
                    }
                }
            }

            if (bestValue < currentValue){
                switch(op){
                case(OP_INSERT):
                    insertOperator(permutation, takePosition, bestOpPosition);
                    break;
                case(OP_EDGEINSERT):
                    edgeInsertOperator(permutation, takePosition, bestOpPosition);
                    break;
                case(OP_SWAP):
                    swapOperator(permutation, takePosition, bestOpPosition);
                    break;
                }
                if (bestValue < currentValue){
                    permutationHasImproved = true;
                }
                currentValue = bestValue;
            }
        }
        return permutationHasImproved;
    }

    ResultData evaluateSolutionMatrix(ProblemData& problemData,
                                      const std::vector<MyGraph::Node>& solution,
                                      std::string targetCriterion = "value",
                                      bool forceLogging = false,
                                      bool abortOnInvalidity = true){
        return evaluateSolution(problemData_, solution, targetCriterion, forceLogging, abortOnInvalidity, &distanceMatrix_);
    }
    double getAndLogDistance(int i, int j){
        additionalLogData_.numberOfShortestPathCalls_++;
        return distanceMatrix_(i,j);
    }

    void callLinKernighanHeuristic(std::vector<MyGraph::Node>& solution){
        // the paper uses Repeated Lin-Kernighan (Applegate et al., 2003)

        // Lin-Kernighan does not really work for paths with only two nodes, so a
        // simple exchange heuristic is used
        if (solution.size() < 2){
            return;
        }
        if (solution.size() == 2){
            std::vector<MyGraph::Node> tempSolution;
            tempSolution.push_back(solution[1]);
            tempSolution.push_back(solution[0]);
            ResultData tempResult = evaluateSolutionMatrix(problemData_, tempSolution);
            ResultData oldResult = evaluateSolutionMatrix(problemData_, solution);
            if (tempResult.length_ < oldResult.length_){
                solution.assign(tempSolution.begin(), tempSolution.end());
            }
        } else {
            // Lin-Kernighan requires the start nodes so they have to be added and later removed
            std::vector<MyGraph::Node> tempSolution;
            tempSolution.push_back(problemData_.startNode_);
            tempSolution.insert(tempSolution.end(), solution.begin(), solution.end());
            tempSolution.push_back(problemData_.startNode_);

            std::vector<MyGraph::Node> solution2 = callRepeatedLinKernighan(problemData_, tempSolution, distanceMatrix_, *this);
            solution.assign(solution2.begin()+1, solution2.end()-1);
        }
    }


};

//' @export
// [[Rcpp::export]]
void callVNSSolver(const Rcpp::DataFrame& nodeDf,
                                    const Rcpp::DataFrame& arcDf,
                                    const Rcpp::DataFrame& problemDf,
                                    double budget,
                                    std::string problemName,
                                    unsigned int runNumber,
                                    std::string fileSuffix = "",
                                    std::string pathToChanges = "",
                                    std::string pathToDistanceMatrix = "",
                                    bool withLocalSearch = false,
                                    bool withLinKernighan = true,
                                    double perturbationProbability = 0.15){

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

    Rcpp::Environment myEnvironment(MYPACKAGE);
    VariableNeighborhoodSearch vns(pd, problemName, runNumber, myEnvironment, "value", fileSuffix, pathToDistanceMatrix);

    vns.withLocalSearch_ = withLocalSearch;
    vns.withLinKernighan_ = withLinKernighan;
    // vns.iterationsWithForcedCriterionRatio_ = iterationsWithForcedCriterionRatio;
    vns.perturbationProbability_ = perturbationProbability;
    vns.calledAsImprover_ = false;
    vns.run();

    //////

    std::vector<MyGraph::Node> mySolution = vns.asCompleteSolution(vns.bestSolution_);
    vns.writeSolution(mySolution, true);
}


//' @export
// [[Rcpp::export]]
void callVNSImprover(const Rcpp::DataFrame& nodeDf,
                           const Rcpp::DataFrame& arcDf,
                           const Rcpp::DataFrame& problemDf,
                           double budget,
                           std::string problemName,
                           unsigned int runNumber,
                           std::string pathToInitialSolution,
                           std::string fileSuffix = "",
                           std::string pathToChanges = "",
                           std::string pathToDistanceMatrix = "",
                           bool withLocalSearch = false,
                           bool withLinKernighan = true,
                           double perturbationProbability = 0.15){

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

    Rcpp::Environment myEnvironment(MYPACKAGE);
    VariableNeighborhoodSearch vns(pd, problemName, runNumber, myEnvironment, "value", fileSuffix, pathToDistanceMatrix);

    vns.withLocalSearch_ = withLocalSearch;
    vns.withLinKernighan_ = withLinKernighan;
    // vns.iterationsWithForcedCriterionRatio_ = iterationsWithForcedCriterionRatio;
    vns.perturbationProbability_ = perturbationProbability;
    vns.calledAsImprover_ = true;
    vns.readInitialSolutionFromFile(pathToInitialSolution);
    vns.run();

    /////////////////////////////

    std::vector<MyGraph::Node> mySolution = vns.asCompleteSolution(vns.bestSolution_);
    vns.writeSolution(mySolution, false);
}
