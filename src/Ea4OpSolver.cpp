#include <Rcpp.h>
#include <lemon/lgf_writer.h>
#include <lemon/list_graph.h>
#include <algorithm>
#include <cmath>
#include <random>

#include "Solver.h"
#include "PlotRoute.h"
#include "IO.h"
#include "ProblemData.h"
#include "CalculateLength.h"
#include "HelperFunctions.h"
#include "SimpleMethods.h"
#include "DistanceMatrix.h"
#include "linkern_wrapper.h"

#define DEBUG1 false
#define DEBUG2 false
#define DEBUG_MAIN_ALGORITHM false
#define DEBUG_CROSSOVER false
#define DEBUG_ADD false

#define LARGE_NUMBER 100000

using namespace lemon;

struct NodeIndexWithDistance{
    int destinationIndex_;
    int indexInSolution_;
    int nodeId_;
    double distance_;
    NodeIndexWithDistance(): destinationIndex_(-1), indexInSolution_(-1), nodeId_(-1), distance_(-1){

    }
    NodeIndexWithDistance(int d1, int s1, int id, double d2): destinationIndex_(d1), indexInSolution_(s1), nodeId_(id), distance_(d2){

    }
    friend std::ostream& operator<<(std::ostream& os, const NodeIndexWithDistance& ald){
        os << "(destIndex, solIndex, id, distance): " << ald.destinationIndex_ << ", " << ald.indexInSolution_ << ", " << ald.nodeId_ << ", " << ald.distance_;
        return os;
    }
};
bool compareNodeWithDistance(NodeIndexWithDistance n1, NodeIndexWithDistance n2){
    return n1.distance_ < n2.distance_;
}

struct ThreeNearest{
    int thisNodeDestinationIndex_;
    int thisNodeIndexInUnvisited_;
    int thisNodeId_;
    NodeIndexWithDistance neigh1_;
    NodeIndexWithDistance neigh2_;
    NodeIndexWithDistance neigh3_;
    double addCost_;
    double addValue_;
    int minInsertPosition;
    int numberOfNeighborsFound_;

    ThreeNearest():
        thisNodeDestinationIndex_(-1), thisNodeIndexInUnvisited_(-1), thisNodeId_(-1),
        neigh1_(NodeIndexWithDistance()), neigh2_(NodeIndexWithDistance()), neigh3_(NodeIndexWithDistance()), addCost_(-1), addValue_(-1),
        minInsertPosition(-1), numberOfNeighborsFound_(-1) {

    }
    ThreeNearest(int nd, int nu, int id, NodeIndexWithDistance n1, NodeIndexWithDistance n2, NodeIndexWithDistance n3):
        thisNodeDestinationIndex_(nd), thisNodeIndexInUnvisited_(nu), thisNodeId_(id),
        neigh1_(n1), neigh2_(n2), neigh3_(n3), addCost_(0.0), addValue_(0.0),
        minInsertPosition(0), numberOfNeighborsFound_(3) {

    }

    friend std::ostream& operator<<(std::ostream& os, const ThreeNearest& tn){
        os << "thisNodeDestinationIndex: " << tn.thisNodeDestinationIndex_ << std::endl;
        os << "thisNodeIndexInUnvisited: " << tn.thisNodeIndexInUnvisited_ << std::endl;
        os << "thisNodeId: " << tn.thisNodeId_ << std::endl;
        os << "neigh1: " << tn.neigh1_ << std::endl;
        os << "neigh2: " << tn.neigh2_ << std::endl;
        os << "neigh3: " << tn.neigh3_ << std::endl;
        os << "numberOfNeighborsFound: " << tn.numberOfNeighborsFound_ << std::endl;
        os << "addCost: " << tn.addCost_ << std::endl;
        os << "addValue: " << tn.addValue_ << std::endl;
        os << "minInsertPosition: " << tn.minInsertPosition << " (actually -1) " << std::endl;
        return os;
    }
};

struct CommonNode{
    MyGraph::Node node_;
    std::vector<MyGraph::Node> connectedNodes_;
    int degree_;
    std::vector<std::vector<MyGraph::Node>> intermediatePaths_;

    CommonNode(MyGraph::Node node) : node_(node), degree_(0){

    }
};

struct Ea4OpSolver : public Solver{
    std::vector<MyGraph::Node> currentSolution_;

    // a second "best" solution which is needed in the algorithm but works differently from the best solution defined by Solver
    // std::vector<MyGraph::Node> tempBestSolution_;
    // ResultData tempBestSolutionRes_;

    std::vector<std::vector<MyGraph::Node>> population_;
    std::vector<ResultData> populationQualities_;

    std::size_t npop_;
    std::size_t ncand_;
    double p_;
    double mutationProbability_; // pmut in the paper
    int d2d_;
    double dominanceOfInitialSolution_;

    // Rcpp::NumericMatrix distanceMatrix_;

    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<double> randomNumber; //Intervall [0,1)

    int worstSolutionIndex_;
    ResultData worstSolutionRes_;
    // ResultData currentSolutionRes_;

    Ea4OpSolver(ProblemData& problemData,
                std::string logFileName, unsigned int runNumber, Rcpp::Environment rEnvironment,
                std::string targetCriterion, std::string fileSuffix, std::string pathToDistanceMatrix)

        : Solver(problemData, logFileName, runNumber, rEnvironment,
          targetCriterion, "ea", fileSuffix, pathToDistanceMatrix) {

        gen = std::mt19937(rd());
        p_ = 0.5; // Seite 13: aber konkrete Daten fehlen, also einfach 0.5 gesetzt
        npop_ = 100;
        ncand_ = 10;
        d2d_ = 50;
        randomNumber = std::uniform_real_distribution<double>(0.0, 1.0);
        // distanceMatrix_ = calculateDistanceMatrix(problemData_);

        worstSolutionIndex_ = -1;
        dominanceOfInitialSolution_ = 0;
        mutationProbability_ = 0.01;
    }

    void run(){

        Rcpp::Rcout << "start " << algorithmName_ << getRandomNumber(0,0) <<  std::endl;

        // calculate parameter p_

        std::vector<MyGraph::Node> allNodes(problemData_.destinations_.begin(), problemData_.destinations_.end());
        tourImprovementOperator(allNodes);
        ResultData allNodesResultData = evaluateSolutionMatrix(problemData_, allNodes);
        p_ = std::sqrt(problemData_.budget_ / allNodesResultData.length_);
        Rcpp::Rcout << "p_ calculated as: " << "sqrt(" << problemData_.budget_ << " / " << allNodesResultData.length_ << ") = " << p_ << "\n";

        if (initialSolution_.empty()){

            waitForInput("constructInitialPopulation", DEBUG_MAIN_ALGORITHM);
            additionalLogData_.currentPhase_ = 0;
            constructInitialPopulation();

            // int counter = 0;
            for (std::vector<MyGraph::Node> sol : population_){
                // Rcpp::Rcout << counter++ << " / " << population_.size() << std::endl;
                // Rcpp::Rcout << "current solution: ";
                // printNodeIdsOfVector(problemData_, sol);
                // Rcpp::Rcout << "ti-";
                waitForInput("tourImprovement", DEBUG_MAIN_ALGORITHM);
                additionalLogData_.currentPhase_ = 1;
                tourImprovementOperator(sol);


                // printNodeIdsOfVector(problemData_, sol);
                // printNodeIdsOfVector(problemData_, tempSol);

                // Rcpp::Rcout << "start: " << problemData_.nodeMap_[problemData_.startNode_].id_ << std::endl;
                // Rcpp::Rcout << "destinations: ";
                // printNodeIdsOfVector(problemData_, problemData_.destinations_);
                // Rcpp::Rcout << "dro-";
                waitForInput("drop", DEBUG_MAIN_ALGORITHM);
                additionalLogData_.currentPhase_ = 2;
                dropOperator(sol);
                // printNodeIdsOfVector(problemData_, sol);
                // Rcpp::Rcout << "ad-";
                waitForInput("add", DEBUG_MAIN_ALGORITHM);
                additionalLogData_.currentPhase_ = 3;
                addOperator(sol);
                // printNodeIdsOfVector(problemData_, sol);
            }

            waitForInput("before loop", DEBUG_MAIN_ALGORITHM);
            int i = 0;
            do{
                i++;
                // Rcpp::Rcout << "i, i mod: " << i << ", " << i % d2d_ << std::endl;
                if (i % d2d_ == 0){
                    waitForInput("mod 0 part", DEBUG_CROSSOVER);
                    std::vector<MyGraph::Node> parent1;
                    std::vector<MyGraph::Node> parent2;

                    // Rcpp::Rcout << "sp-";
                    waitForInput("select parents", DEBUG_CROSSOVER);
                    additionalLogData_.currentPhase_ = 4;
                    selectParents(parent1, parent2);

                    // Rcpp::Rcout << "cr-";
                    waitForInput("crossover", DEBUG_CROSSOVER);
                    additionalLogData_.currentPhase_ = 5;
                    std::vector<MyGraph::Node> childSolution = crossoverOperator(parent1, parent2);

                    // printNodeIdsOfVector(problemData_, childSolution);

                    double r = randomNumber(gen);
                    if (r < mutationProbability_){
                        Rcpp::Rcout << "mu-";
                        waitForInput("mutate", DEBUG_CROSSOVER);
                        additionalLogData_.currentPhase_ = 6;
                        mutationOperator(childSolution);
                    }

                    ResultData childResult = evaluateSolutionMatrix(problemData_, childSolution);
                    if (childResult > worstSolutionRes_){
                        waitForInput("replace worst", DEBUG_CROSSOVER);
                        worstSolutionRes_ = childResult;
                        population_[worstSolutionIndex_].assign(childSolution.begin(), childSolution.end());
                        populationQualities_[worstSolutionIndex_] = childResult;

                        // update worst solution in population
                        for (std::size_t i = 0; i < population_.size(); i++){
                            if (populationQualities_[i] < worstSolutionRes_){
                                worstSolutionRes_ = populationQualities_[i];
                                worstSolutionIndex_ = i;
                            }
                        }
                    }
                } else {
                    int counter = 0;
                    for (std::vector<MyGraph::Node> sol : population_){
                        // Rcpp::Rcout << "internal loop: " << counter++ << " / " << population_.size() << std::endl;
                        // Rcpp::Rcout << "ti-";
                        waitForInput("tourImprovement", DEBUG_MAIN_ALGORITHM);
                        additionalLogData_.currentPhase_ = 11;
                        tourImprovementOperator(sol);
                        
                        // Rcpp::Rcout << "dro-";
                        waitForInput("drop", DEBUG_MAIN_ALGORITHM);
                        additionalLogData_.currentPhase_ = 12;
                        dropOperator(sol);
                        // Rcpp::Rcout << "ad-";
                        waitForInput("add", DEBUG_MAIN_ALGORITHM);
                        additionalLogData_.currentPhase_ = 13;
                        addOperator(sol);
                        if (terminationCriterionSatisfied()) {
                            Rcpp::Rcout << "Info: Termination criterion satisfied. Break out of inner loop." << std::endl;
                            break;
                        }
                    }
                }

                additionalLogData_.currentIteration_++;
                Rcpp::Rcout << ".";
                Rcpp::checkUserInterrupt();
            } while(!terminationCriterionSatisfied());
            Rcpp::Rcout << bestSolutionQuality_ << std::endl;
        } else {
            // initial solution is given
            waitForInput("constructInitialPopulation", DEBUG_MAIN_ALGORITHM);
            additionalLogData_.currentPhase_ = 0;
            std::size_t old_npop = npop_;
            npop_ *= (1 - dominanceOfInitialSolution_);
            Rcpp::Rcout << "new npop: " << npop_ << std::endl;


            constructInitialPopulation();

            /*
            * initial population here: part of the population initialised as usual,
            * rest with the given initial solution and perturbed variants
            * (the ratio of these two partial populations is given by the parameter dominanceOfInitialSolution_)
            */

            npop_ = old_npop;
            Rcpp::Rcout << "old npop: " << npop_ << std::endl;
            population_.push_back(initialSolution_);
            ResultData tempResult = evaluateSolutionMatrix(problemData_, initialSolution_);

            if (tempResult < worstSolutionRes_){
                worstSolutionRes_.value_ = tempResult.value_;
                worstSolutionRes_.length_ = tempResult.length_;
                worstSolutionIndex_ = population_.size() - 1;
            }

            additionalLogData_.currentPhase_ = 10;
            while (population_.size() < npop_){
                std::vector<MyGraph::Node> tempSolution(initialSolution_.begin(), initialSolution_.end());
                mutationOperator(tempSolution);

                population_.push_back(tempSolution);
                ResultData tempResult = evaluateSolutionMatrix(problemData_, tempSolution);

                // Rcpp::Rcout << "tempResult: " << tempResult << std::endl;
                populationQualities_.push_back(tempResult);
                // remember worst solution so far
                if (tempResult < worstSolutionRes_){
                    worstSolutionRes_.value_ = tempResult.value_;
                    worstSolutionRes_.length_ = tempResult.length_;
                    worstSolutionIndex_ = population_.size() - 1;
                }
            }


            Rcpp::Rcout << "Population loop" << std::endl;
            for (std::vector<MyGraph::Node> sol : population_){
                // Rcpp::Rcout << "ti-";
                waitForInput("tourImprovement", DEBUG_MAIN_ALGORITHM);
                additionalLogData_.currentPhase_ = 1;
                tourImprovementOperator(sol);

                // Rcpp::Rcout << "DE1: " << additionalLogData_.numberOfShortestPathCalls_ << ", best1: " << bestSolutionQuality_.value_ << "\n";
                // evaluateSolutionMatrix(problemData_, sol, "value", true);
                // Rcpp::Rcout << "dro-";
                waitForInput("drop", DEBUG_MAIN_ALGORITHM);
                additionalLogData_.currentPhase_ = 2;
                dropOperator(sol);

                // evaluateSolutionMatrix(problemData_, sol, "value", true);
                // Rcpp::Rcout << "DE2: " << additionalLogData_.numberOfShortestPathCalls_ << ", best2: " << bestSolutionQuality_.value_ << "\n";
                // Rcpp::stop("myEnd");

                // Rcpp::Rcout << "ad-";
                waitForInput("add", DEBUG_MAIN_ALGORITHM);
                additionalLogData_.currentPhase_ = 3;
                addOperator(sol);

                Rcpp::Rcout << ".";
            }



            Rcpp::Rcout << std::endl << "Big loop" << std::endl;
            waitForInput("before loop", DEBUG_MAIN_ALGORITHM);
            int i = 0;
            do{
                i++;
                Rcpp::Rcout << "i, i mod: " << i << ", " << i % d2d_ << std::endl;
                if (i % d2d_ == 0){
                    waitForInput("mod 0 part", DEBUG_CROSSOVER);
                    std::vector<MyGraph::Node> parent1;
                    std::vector<MyGraph::Node> parent2;

                    // Rcpp::Rcout << "sp-";
                    waitForInput("select parents", DEBUG_CROSSOVER);
                    additionalLogData_.currentPhase_ = 4;
                    selectParents(parent1, parent2);

                    // Rcpp::Rcout << "cr-";
                    waitForInput("crossover", DEBUG_CROSSOVER);
                    additionalLogData_.currentPhase_ = 5;
                    std::vector<MyGraph::Node> childSolution = crossoverOperator(parent1, parent2);

                    // printNodeIdsOfVector(problemData_, childSolution);

                    double r = randomNumber(gen);
                    if (r < mutationProbability_){
                        Rcpp::Rcout << "mu-";
                        waitForInput("mutate", DEBUG_CROSSOVER);
                        additionalLogData_.currentPhase_ = 6;
                        mutationOperator(childSolution);
                    }

                    additionalLogData_.currentPhase_ = 7;
                    ResultData childResult = evaluateSolutionMatrix(problemData_, childSolution);
                    if (childResult > worstSolutionRes_){
                        waitForInput("replace worst", DEBUG_CROSSOVER);
                        worstSolutionRes_ = childResult;
                        population_[worstSolutionIndex_].assign(childSolution.begin(), childSolution.end());
                        populationQualities_[worstSolutionIndex_] = childResult;

                        // update worst solution in population
                        for (std::size_t i = 0; i < population_.size(); i++){
                            // Rcpp::Rcout << populationQualities_[i] << std::endl;

                            if (populationQualities_[i] < worstSolutionRes_){
                                worstSolutionRes_ = populationQualities_[i];
                                worstSolutionIndex_ = i;
                            }
                        }
                    }
                } else {
                    int counter = 0;
                    for (std::vector<MyGraph::Node> sol : population_){
                        //Rcpp::Rcout << "internal loop: " << counter++ << " / " << population_.size() << std::endl;

                        // Rcpp::Rcout << "ti-";
                        waitForInput("tourImprovement", DEBUG_MAIN_ALGORITHM);
                        additionalLogData_.currentPhase_ = 11;
                        tourImprovementOperator(sol);

                        // Rcpp::Rcout << "dro-";
                        waitForInput("drop", DEBUG_MAIN_ALGORITHM);
                        additionalLogData_.currentPhase_ = 12;
                        dropOperator(sol);

                        // Rcpp::Rcout << "ad-";
                        waitForInput("add", DEBUG_MAIN_ALGORITHM);
                        additionalLogData_.currentPhase_ = 13;
                        addOperator(sol);
                        if (terminationCriterionSatisfied()) {
                            Rcpp::Rcout << "Info: Termination criterion satisfied. Break out of inner loop." << std::endl;
                            break;
                        }
                    }
                }
                additionalLogData_.currentIteration_++;
                Rcpp::Rcout << ".";
                Rcpp::checkUserInterrupt();
            } while(!terminationCriterionSatisfied());
        }

        evaluateSolutionMatrix(problemData_, bestSolution_, "value", true);
        Rcpp::Rcout << "end " << algorithmName_ << std::endl;
    }

    void constructInitialPopulation(){
        worstSolutionRes_.length_ = 0;
        worstSolutionRes_.value_ = std::numeric_limits<double>::max();

        for (std::size_t i = 0; i < npop_; i++){
            std::vector<MyGraph::Node> tempSolution;
            for (MyGraph::Node n : problemData_.destinations_){
                double myRandomNumber = randomNumber(gen);
                if (myRandomNumber < p_){
                    tempSolution.push_back(n);
                }
            }
            std::shuffle(tempSolution.begin(), tempSolution.end(), gen);
            // Rcpp::Rcout << "push into initial population: ";
            // printNodeIdsOfVector(problemData_, tempSolution);
            population_.push_back(tempSolution);
            additionalLogData_.numberOfTestedBitVectors_++;
            ResultData tempResult = evaluateSolutionMatrix(problemData_, tempSolution, "value", true);

            // Rcpp::Rcout << "tempResult: " << tempResult << std::endl;
            populationQualities_.push_back(tempResult);
            // remember worst solution so far
            if (tempResult < worstSolutionRes_){
                worstSolutionRes_.value_ = tempResult.value_;
                worstSolutionRes_.length_ = tempResult.length_;
                worstSolutionIndex_ = i;
            }

        }


        ////// Test population

        // population_.clear();
        // std::vector<int> bla = {300,800};
        // std::vector<MyGraph::Node> tempSol;
        // for (int x : bla){
        //     for (MyGraph::Node d : problemData_.destinations_){
        //         if (problemData_.nodeMap_[d].id_ == x){
        //             tempSol.push_back(d);
        //         }
        //     }
        // }
        // population_.push_back(tempSol);
        // population_.push_back(tempSol);
        // population_.push_back(tempSol);
    }
    void mutationOperator(std::vector<MyGraph::Node>& solution){

        // Rcpp::Rcout << "before mutation: ";
        // printNodeIdsOfVector(problemData_, solution);
        std::size_t r = getRandomNumber(0, problemData_.destinations_.size() - 1);
        auto posOfElement = std::find(solution.begin(), solution.end(), problemData_.destinations_[r]);
        if (posOfElement == solution.end()){

            // Rcpp::Rcerr << "mutation: insert node " << problemData_.nodeMap_[problemData_.destinations_[r]].id_ << std::endl;

            // find 3 nearest nodes
            std::vector<int> solutionIndices = nodeVectorToIndexVector(problemData_, solution);
            std::vector<NodeIndexWithDistance> currentDistances;

            // add start node so that at least 3 closest nodes can be found
            solutionIndices.insert(solutionIndices.begin(), 0);

            // printVector(solutionIndices);
            // printNodeIdsOfVector(problemData_, solution);

            for (int j = 0; j < (int) solutionIndices.size(); j++){
                int nodeIdToWrite = (solutionIndices[j] == 0) ? problemData_.nodeMap_[problemData_.startNode_].id_ : problemData_.nodeMap_[problemData_.destinations_[solutionIndices[j]-1]].id_;

                currentDistances.push_back(
                    NodeIndexWithDistance(
                                          solutionIndices[j],//r,
                                          j,
                                          nodeIdToWrite,//problemData_.nodeMap_[problemData_.destinations_[r]].id_,
                                          getAndLogDistance(r, solutionIndices[j])
                                        )
                );
            }
            std::sort(currentDistances.begin(), currentDistances.end(), compareNodeWithDistance);

            // if the solution is (start, node, start), i.e. only 2 different nodes then add the start node
            if (currentDistances.size() < 3){
                currentDistances.push_back(NodeIndexWithDistance(-11, -11, -11, -11));
            }

            // check adjacency and insert
            ThreeNearest tn(r, -1,
                            problemData_.nodeMap_[problemData_.destinations_[r]].id_,
                            currentDistances[0],
                            currentDistances[1],
                            currentDistances[2]
            );

            if (tn.neigh3_.destinationIndex_ == -11){
                tn.numberOfNeighborsFound_ = 2;
            } else {
                tn.numberOfNeighborsFound_ = 3;
            }

            threeNearestCalculateMinCost(tn, solutionIndices);
            // Rcpp::Rcout << "tn.minInsertPosition: " << tn.minInsertPosition  << " / size: " << solutionIndices.size() << std::endl;

            solution.insert(solution.begin() + tn.minInsertPosition - 1, problemData_.destinations_[r]);
            // printNodeIdsOfVector(problemData_, solution);
        } else {
            // Rcpp::Rcout << "mutation: remove node " << problemData_.nodeMap_[*posOfElement].id_ << std::endl;
            solution.erase(posOfElement);
        }

        // Rcpp::Rcout << "after mutation: ";
        // printNodeIdsOfVector(problemData_, solution);

        additionalLogData_.numberOfTestedBitVectors_++;
    }

    void selectParents(std::vector<MyGraph::Node>& parent1, std::vector<MyGraph::Node>& parent2){
        double minFitness = std::numeric_limits<double>::max();
        std::vector<ResultData> resultDataContainer;

        std::vector<int> candidateIndices;
        while(candidateIndices.size() < ncand_){
            int randomIndex = getRandomNumber(0, population_.size()-1);
            if (std::find(candidateIndices.begin(), candidateIndices.end(), randomIndex) == candidateIndices.end()){
                candidateIndices.push_back(randomIndex);
            }
        }

        // find value of solution with minimum fitness (value)
        for (int i : candidateIndices){
            ResultData result = evaluateSolutionMatrix(problemData_, population_[i]);
            resultDataContainer.push_back(result);

            // update population qualities (in case there was a change in the instance) and worst individual
            /*
             * To-Do: this is only for the case where the solution after the change is still valid...
             */
            populationQualities_[i] = result;
            if (result < worstSolutionRes_) {
                worstSolutionRes_ = result;
                worstSolutionIndex_ = i;
            }

            if (result.value_ < minFitness){
                minFitness = result.value_;
            }
        }

        // calculate r_i
        std::vector<double> rValues;
        double totalR = 0;
        for (std::size_t i = 0; i < resultDataContainer.size(); i++){
            double r_i = resultDataContainer[i].value_ - minFitness + 1;
            rValues.push_back(r_i);
            totalR += r_i;
        }

        // calculate p_i
        std::vector<double> probabilities;
        std::vector<double> cumulativeProbabilities;

        for (std::size_t i = 0; i < resultDataContainer.size(); i++){
            double p_i = rValues[i] / totalR;
            if (cumulativeProbabilities.empty()){
                cumulativeProbabilities.push_back(p_i);
            } else {
                cumulativeProbabilities.push_back(cumulativeProbabilities.back() + p_i);
            }
            // Rcpp::Rcout << "cumulative: " << cumulativeProbabilities.back() << std::endl;
            probabilities.push_back(p_i);
        }
        cumulativeProbabilities.push_back(1.0);
        // Rcpp::Rcout << "cumulative: " << cumulativeProbabilities.back() << std::endl;

        // sample two parent solutions (which do not have to be different)
        double myRandomNumber = randomNumber(gen);
        int index = findIndexInCumulativeProbabilities(cumulativeProbabilities, myRandomNumber);
        // Rcpp::Rcout << "index1: " << index << std::endl;
        parent1.assign(population_[candidateIndices[index]].begin(), population_[candidateIndices[index]].end());

        myRandomNumber = randomNumber(gen);
        index = findIndexInCumulativeProbabilities(cumulativeProbabilities, myRandomNumber);
        // Rcpp::Rcout << "index2: " << index << std::endl;
        parent2.assign(population_[candidateIndices[index]].begin(), population_[candidateIndices[index]].end());

    }

    void tourImprovementOperator(std::vector<MyGraph::Node>& solution){
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
        evaluateSolutionMatrix(problemData_, solution);
    }

    std::vector<MyGraph::Node> crossoverOperator(std::vector<MyGraph::Node> parent1, std::vector<MyGraph::Node> parent2){

        waitForInput("startCrossover", DEBUG1);

        // Rcpp::Rcout << std::endl << "parents:" << std::endl;
        // printNodeIdsOfVector(problemData_, parent1);
        // printNodeIdsOfVector(problemData_, parent2);

        ////////////////////////////////////////
        // Example from the paper
        ////////////////////////////////////////
        // std::vector<int> parent1b = {4, 5, 6, 7, 13, 12, 11, 10, 9, 8};
        // std::vector<int> parent2b = {2, 4, 10, 14, 11, 6, 3, 7, 12, 16, 15, 9};
        //
        // parent1.clear();
        // parent2.clear();
        //
        // for (int n : parent1b){
        //     parent1.push_back(problemData_.destinations_[n]);
        // }
        //
        // for (int n : parent2b){
        //     parent2.push_back(problemData_.destinations_[n]);
        // }

        waitForInput("search common nodes", DEBUG1);

        // startNode has to be added
        parent1.insert(parent1.begin(), problemData_.startNode_);
        parent2.insert(parent2.begin(), problemData_.startNode_);

        std::vector<MyGraph::Node> nodesInBothGraphs;

        for (std::size_t i = 0; i < parent1.size(); i++){
            MyGraph::Node& n = parent1[i];
            auto positionInParent2 = std::find(parent2.begin(), parent2.end(), n);
            if (positionInParent2 != parent2.end()){
                nodesInBothGraphs.push_back(n);
            }
        }

        // special cases that are excluded
        if (nodesInBothGraphs.size() == 1 || parent1 == parent2){
            std::vector<MyGraph::Node> result;
            double r = randomNumber(gen);
            if (r < 0.5){
                result.assign(parent1.begin()+1, parent1.end()); // +1 since the start node is not part of the solution
            } else {
                result.assign(parent2.begin()+1, parent2.end());
            }
            // Rcpp::Rcout << "trivial case; crossover ends prematurely" << std::endl;
            return result;
        }


        // Rcpp::Rcout << "nodes in both graphs:";
        // for (MyGraph::Node n : nodesInBothGraphs) {
        //     Rcpp::Rcout << problemData_.nodeMap_[n].id_ << " ";
        // }
        // Rcpp::Rcout << std::endl;
        //
        // Rcpp::Rcout << "all destinations:";
        // for (MyGraph::Node n : problemData_.destinations_) {
        //     Rcpp::Rcout << problemData_.nodeMap_[n].id_ << " ";
        // }
        // Rcpp::Rcout << std::endl;
        // Rcpp::Rcout << std::endl;


        std::vector<CommonNode> edgeMap = constructEdgeMap(parent1, parent2, nodesInBothGraphs);



        /*
        * To-Do: Crossover operator
        * Beginne mit startNode, entferne aus edgeMap.
        * Wähle dann den nächsten adjazenten Knoten mit minimalen Grad.
        *  Nulldifferenzen werden zufällig aufgelöst.
        * Wenn keine adjazenten Knoten da: Wähle zufällig einen anderen Knoten in edgeMap
        *
        * Für jeden gewählten Knoten: nimm zufälligen intermediatePath zu diesem Knoten.
        *
        * Terminieren, wenn edgeMap leer ist. Dann wird letzter besuchter Knoten mit erstem Knoten verbunden
        */

        waitForInput("begin crossover", DEBUG2);
        std::vector<MyGraph::Node> result;
        CommonNode cn = edgeMap[0];
        int currentCommonNodePosition = 0;

        do {

            int minDegree = std::numeric_limits<int>::max();
            std::vector<CommonNode*> commonNodesWithMinDegree;
            std::vector<int> positionsOfCommonNodesWithMinDegree;
            int counter = 0;

            // Rcpp::Rcout << "edgeMap has: ";
            // for (CommonNode& nn : edgeMap){
            //     Rcpp::Rcout << problemData_.nodeMap_[nn.node_].id_ << ", ";
            // }
            // Rcpp::Rcout << std::endl;

            // Rcpp::Rcout << "Current commonNode:" << std::endl;
            // Rcpp::Rcout << problemData_.nodeMap_[cn.node_].id_ << std::endl;
            waitForInput("search connected common nodes", DEBUG2);
            for (CommonNode edgeMapNode : edgeMap){
                // suchen in commonNode
                auto posInConnectedNodes = std::find(cn.connectedNodes_.begin(), cn.connectedNodes_.end(), edgeMapNode.node_);
                if (posInConnectedNodes != cn.connectedNodes_.end()){

                    if (edgeMapNode.degree_ == minDegree){
                        commonNodesWithMinDegree.push_back(&edgeMapNode);
                        positionsOfCommonNodesWithMinDegree.push_back(counter);
                    }
                    if (edgeMapNode.degree_ < minDegree){
                        minDegree = edgeMapNode.degree_;
                        commonNodesWithMinDegree.clear();
                        commonNodesWithMinDegree.push_back(&edgeMapNode);
                        positionsOfCommonNodesWithMinDegree.clear();
                        positionsOfCommonNodesWithMinDegree.push_back(counter);
                    }
                }

                if (cn.node_ == edgeMapNode.node_){
                    currentCommonNodePosition = counter;
                }
                counter++;
            }


            int nextPosition = 0;
            int rNextNode = 0;
            if (commonNodesWithMinDegree.empty()) {
                waitForInput("no connected common nodes", DEBUG2);
                do {
                    nextPosition = getRandomNumber(0, edgeMap.size()-1);
                } while (nextPosition == currentCommonNodePosition);
                result.push_back(edgeMap[nextPosition].node_);
            } else {
                waitForInput("There are connected common nodes. Choose random intermediate path.", DEBUG2);
                rNextNode = getRandomNumber(0, commonNodesWithMinDegree.size() - 1);
                nextPosition = positionsOfCommonNodesWithMinDegree[rNextNode];

                MyGraph::Node& nextNode = edgeMap[nextPosition].node_;

                // Rcpp::Rcout << "nextNode: " << problemData_.nodeMap_[nextNode].id_ << std::endl;

                // choose an intermediate path
                std::vector<std::vector<MyGraph::Node>> pathsToConnectedNode;
                for (std::vector<MyGraph::Node> path : cn.intermediatePaths_){

                    if (*(path.end()-1) == nextNode){
                        pathsToConnectedNode.push_back(path);
                        // printNodeIdsOfVector(problemData_, path);
                    }
                }

                int rPath = getRandomNumber(0, pathsToConnectedNode.size() - 1);
                result.insert(result.end(), pathsToConnectedNode[rPath].begin()+1, pathsToConnectedNode[rPath].end());
                // Rcpp::Rcout << "chosenPath: >> ";
                // for (MyGraph::Node n : pathsToConnectedNode[rPath]){
                //     Rcpp::Rcout << problemData_.nodeMap_[n].id_ << " ";
                // }
                // Rcpp::Rcout << "<<" << std::endl;

                // printNodeIdsOfVector(problemData_, result);
            }

            waitForInput("recompute degrees of all connected nodes", DEBUG2);
            // recompute degrees of all connected nodes
            for (CommonNode edgeMapNode : edgeMap){

                // assumption: we have a graph where every arc can be traversed in both directions
                // so that we can conclude that cn is connected to edgeMapNode iff edgeMapNode is connected to cn
                auto findPos = std::find(edgeMapNode.connectedNodes_.begin(), edgeMapNode.connectedNodes_.end(), cn.node_);
                if (findPos != edgeMapNode.connectedNodes_.end()){
                    edgeMapNode.connectedNodes_.erase(findPos);
                    edgeMapNode.degree_--;
                }
            }
            waitForInput("reassign cn and delete", DEBUG2);
            cn = edgeMap[nextPosition];
            edgeMap.erase(edgeMap.begin() + currentCommonNodePosition);

        } while (edgeMap.size() > 1);
        waitForInput("crossover end", DEBUG2);
        // Rcpp::Rcout << "crossover result: ";
        // printNodeIdsOfVector(problemData_, result);

        // it is not necessary to go back to the starting node.
        // also, it is not possible to return to the starting node because it is immediately removed
        // as a common node

        // this is a new bit vector
        additionalLogData_.numberOfTestedBitVectors_++;
        return result;
    }

    std::vector<CommonNode> constructEdgeMap(std::vector<MyGraph::Node>& parent1, std::vector<MyGraph::Node>& parent2,
                                             std::vector<MyGraph::Node>& nodesInBothGraphs){

        // Rcpp::Rcout << "nodes in both graphs:";
        // for (MyGraph::Node n : nodesInBothGraphs) {
        //     Rcpp::Rcout << problemData_.nodeMap_[n].id_ << " ";
        // }
        // Rcpp::Rcout << std::endl;

        // Rcpp::Rcout << "all destinations:";
        // for (MyGraph::Node n : problemData_.destinations_) {
        //     Rcpp::Rcout << problemData_.nodeMap_[n].id_ << " ";
        // }
        // Rcpp::Rcout << std::endl;
        // Rcpp::Rcout << std::endl;

        std::vector<CommonNode> edgeMap;
        // for (MyGraph::Node n : nodesInBothGraphs){
        //     CommonNode commonNode(n);
        //
        //     auto positionInParent1 = std::find(parent1.begin(), parent1.end(), n);
        //     auto positionInParent2 = std::find(parent2.begin(), parent2.end(), n);
        //
        //     auto currentPos = positionInParent1;
        //     bool foundCommonNode = false;
        //     auto adjacentCommonNode = nodesInBothGraphs.begin();
        //
        //     /*
        //      * To-Do: Zyklen sind auch noch möglich?
        //      */
        //
        //
        //
        //     // parent1 -- go left
        //     do{
        //         currentPos--;
        //         adjacentCommonNode = std::find(nodesInBothGraphs.begin(), nodesInBothGraphs.end(), *currentPos);
        //         if (adjacentCommonNode != nodesInBothGraphs.end()) {
        //             foundCommonNode = true;
        //         }
        //     } while (!foundCommonNode && currentPos != parent1.begin());
        //
        //     if (foundCommonNode == true) {
        //         commonNode.connectedNodes_.push_back(*adjacentCommonNode);
        //         commonNode.degree_++;
        //     }
        //
        //     currentPos = positionInParent1;
        //     foundCommonNode = false;
        //
        //     // parent 1 -- go right
        //     do{
        //         currentPos++;
        //         adjacentCommonNode = std::find(nodesInBothGraphs.begin(), nodesInBothGraphs.end(), *currentPos);
        //         if (adjacentCommonNode != nodesInBothGraphs.end()) {
        //             foundCommonNode = true;
        //         }
        //     } while (!foundCommonNode && currentPos != parent1.end());
        //
        //     if (foundCommonNode == true) {
        //         commonNode.connectedNodes_.push_back(*adjacentCommonNode);
        //         commonNode.degree_++;
        //     }
        //
        //     currentPos = positionInParent2;
        //     foundCommonNode = false;
        //
        //     // parent 2 -- go left
        //     do{
        //         currentPos--;
        //         adjacentCommonNode = std::find(nodesInBothGraphs.begin(), nodesInBothGraphs.end(), *currentPos);
        //         if (adjacentCommonNode != nodesInBothGraphs.end()) {
        //             foundCommonNode = true;
        //         }
        //     } while (!foundCommonNode && currentPos != parent2.begin());
        //
        //     if (foundCommonNode == true) {
        //         commonNode.connectedNodes_.push_back(*adjacentCommonNode);
        //         commonNode.degree_++;
        //     }
        //
        //     currentPos = positionInParent2;
        //     foundCommonNode = false;
        //
        //     // parent 2 -- go right
        //     do{
        //         currentPos++;
        //         adjacentCommonNode = std::find(nodesInBothGraphs.begin(), nodesInBothGraphs.end(), *currentPos);
        //         if (adjacentCommonNode != nodesInBothGraphs.end()) {
        //             foundCommonNode = true;
        //         }
        //     } while (!foundCommonNode && currentPos != parent2.begin());
        //
        //     if (foundCommonNode == true) {
        //         commonNode.connectedNodes_.push_back(*adjacentCommonNode);
        //         commonNode.degree_++;
        //     }
        //
        //     currentPos = positionInParent2;
        //     foundCommonNode = false;
        //
        //     /*
        //      * intermediate paths: von jedem commonNode aus: suche ihn in parent1, parent2 und gehe in beiden links und
        //      * rechts. Schreibe alle besuchten Knoten auf, sofern sie noch nicht drin sind.
        //      * Das muss wohl irgendwie rekursiv gemacht werden...
        //      */
        //
        //
        //
        // }

        waitForInput("build edgemap", DEBUG1);
        for (MyGraph::Node n : nodesInBothGraphs){
            CommonNode commonNode(n);
            std::vector<MyGraph::Node> startPrefix;
            //startPrefix.push_back(n);

            // Rcpp::Rcout << "common node: " << problemData_.nodeMap_[n].id_ << " : " << problemData_.nodeMap_[commonNode.node_].id_ << std::endl;

            // calculate intermediate paths (as well as adjacent nodes and degrees)
            addIntermediatePathsRecursive(commonNode.intermediatePaths_,
                                          nodesInBothGraphs,
                                          startPrefix,
                                          parent1,
                                          parent2,
                                          n,
                                          commonNode);

            // Rcpp::Rcout << "common node (after): " << problemData_.nodeMap_[n].id_ << " : " << problemData_.nodeMap_[commonNode.node_].id_ << std::endl;

            edgeMap.push_back(commonNode);
        }

        /////////////////////////////////
        // Print edge map
        /////////////////////////////////

        // waitForInput("Start printing edge map", true);
        // Rcpp::Rcout << "----edge map--------" <<  std::endl;
        // for (CommonNode cn : edgeMap){
        //     Rcpp::Rcout << "this node id: " << problemData_.nodeMap_[cn.node_].id_ << std::endl;
        //     Rcpp::Rcout << "degree: " << cn.degree_ << std::endl;
        //
        //     Rcpp::Rcout << "connectedNodes: ";
        //     for (MyGraph::Node n : cn.connectedNodes_) {
        //         Rcpp::Rcout << problemData_.nodeMap_[n].id_ << ", ";
        //     }
        //     Rcpp::Rcout << std::endl;
        //
        //     Rcpp::Rcout << "intermediate paths: " << std::endl;
        //     for (std::vector<MyGraph::Node> path : cn.intermediatePaths_) {
        //         for (MyGraph::Node n : path) {
        //             Rcpp::Rcout << problemData_.nodeMap_[n].id_ << ", ";
        //         }
        //         Rcpp::Rcout << std::endl;
        //     }
        //     Rcpp::Rcout << std::endl;
        //     Rcpp::Rcout << "------------" <<  std::endl;
        // }
        // Rcpp::Rcout << std::endl;
        // Rcpp::Rcout << std::endl;
        // waitForInput("End printing edge map", true);
        // Rcpp::Rcout << "all destinations: ";
        // for (MyGraph::Node n : problemData_.destinations_) {
        //     Rcpp::Rcout << problemData_.nodeMap_[n].id_ << " ";
        // }
        // Rcpp::Rcout << std::endl;
        // Rcpp::Rcout << std::endl;

        return edgeMap;
    }

    void addIntermediatePathsRecursive(std::vector<std::vector<MyGraph::Node>>& setOfPaths,
                                       const std::vector<MyGraph::Node>& commonNodes,
                                       std::vector<MyGraph::Node> prefix,
                                       const std::vector<MyGraph::Node>& parent1,
                                       const std::vector<MyGraph::Node>& parent2,
                                       MyGraph::Node startNode,
                                       CommonNode& currentNode){

        // Rcpp::Rcout << "call recursive with prefix length " << prefix.size() << std::endl;
        /*
        * Wenn Startknoten schon markiert ist, dann Abbruch.
        * Jeder rekursive Aufruf gibt nicht nur einen Pfad, sondern eine Liste von Pfaden, die dann zusammengeführt werden müssen.
        * Input: Präfix - der Pfad zum aktuellen Knoten und Pfadmenge P
        * Output: Menge von Pfaden, die alle den gleichen Präfix haben
        *
        * für jeden zu startNode adjazenten Knoten:
        *  -wenn nicht in prefix drin, dann:
        *      -wenn startnode ein common node: P = P u {pfad+start}
        *      -wenn kein common node: P = P u rekursive
        */

        auto positionInParent1 = std::find(parent1.begin(), parent1.end(), startNode);
        auto positionInParent2 = std::find(parent2.begin(), parent2.end(), startNode);

        std::vector<MyGraph::Node> nodesToCheck;
        if (positionInParent1 != parent1.end()){
            if (positionInParent1 != parent1.begin()){
                nodesToCheck.push_back(*(positionInParent1 - 1));
            } else {
                nodesToCheck.push_back(*(parent1.end()-1));
            }
            if (positionInParent1 + 1 != parent1.end()){
                nodesToCheck.push_back(*(positionInParent1 + 1));
            } else {
                nodesToCheck.push_back(*parent1.begin());
            }
        }

        if (positionInParent2 != parent2.end()){
            if (positionInParent2 != parent2.begin()){
                nodesToCheck.push_back(*(positionInParent2 - 1));
            } else {
                nodesToCheck.push_back(*(parent2.end()-1));
            }
            if (positionInParent2 + 1 != parent2.end()){
                nodesToCheck.push_back(*(positionInParent2 + 1));
            } else {
                nodesToCheck.push_back(*parent2.begin());
            }
        }

        // Rcpp::Rcout << "startNode: " << problemData_.nodeMap_[startNode].id_ << std::endl;
        // Rcpp::Rcout << "nodes to check: ";
        // for (MyGraph::Node n : nodesToCheck) {
        //     Rcpp::Rcout << problemData_.nodeMap_[n].id_ << " ";
        // }
        // Rcpp::Rcout << std::endl;

        // Rcpp::Rcout << "all destinations:";
        // for (MyGraph::Node n : problemData_.destinations_) {
        //     Rcpp::Rcout << problemData_.nodeMap_[n].id_ << " ";
        // }
        // Rcpp::Rcout << std::endl;
        // Rcpp::Rcout << std::endl;

        waitForInput("loop begin", DEBUG1);

        prefix.push_back(startNode);

        for (MyGraph::Node tempNode : nodesToCheck){

            if (std::find(prefix.begin(), prefix.end(), tempNode) == prefix.end()){
                std::vector<MyGraph::Node> tempPrefix(prefix.begin(), prefix.end());

                // startNode is almost never a common node (exception: at the very beginning)
                // tempPrefix.push_back(startNode);
                // printNodeIdsOfVector(problemData_, tempPrefix);
                // Rcpp::Rcout << "tempNode: " << problemData_.nodeMap_[tempNode].id_ << std::endl;

                if (std::find(commonNodes.begin(), commonNodes.end(), tempNode) == commonNodes.end()){
                    waitForInput("not a common node", DEBUG1);
                    // not a common node
                    addIntermediatePathsRecursive(setOfPaths,
                                                  commonNodes,
                                                  tempPrefix,
                                                  parent1,
                                                  parent2,
                                                  tempNode,
                                                  currentNode);
                } else {
                    waitForInput("common node", DEBUG1);

                    tempPrefix.push_back(tempNode);
                    setOfPaths.push_back(tempPrefix);
                    // Rcpp::Rcout << "push back path: ";
                    // printNodeIdsOfVector(problemData_, tempPrefix);


                    if (std::find(currentNode.connectedNodes_.begin(), currentNode.connectedNodes_.end(), tempNode)
                            == currentNode.connectedNodes_.end() ){

                        // std::stringstream ss;
                        // ss << "add node " << problemData_.nodeMap_[startNode].id_ << " to "
                        //    << problemData_.nodeMap_[currentNode.node_].id_ << std::endl;
                        // Rcpp::Rcout << ss.str();
                        // waitForInput(ss.str(), DEBUG1);
                        currentNode.connectedNodes_.push_back(tempNode);
                        currentNode.degree_++;
                    }
                }
            }
        }

        // for (MyGraph::Node tempNode : nodesToCheck){
        //
        //     if (std::find(prefix.begin(), prefix.end(), tempNode) == prefix.end()){
        //         std::vector<MyGraph::Node> tempPrefix(prefix.begin(), prefix.end());
        //
        //         // startNode is almost never a common node (only at the very beginning)
        //         tempPrefix.push_back(startNode);
        //         printNodeIdsOfVector(problemData_, tempPrefix);
        //         Rcpp::Rcout << "tempNode: " << problemData_.nodeMap_[tempNode].id_ << std::endl;
        //
        //         if (std::find(commonNodes.begin(), commonNodes.end(), tempNode) == commonNodes.end()){
        //             waitForInput("not a common node", DEBUG1);
        //             // not a common node
        //             addIntermediatePathsRecursive(setOfPaths,
        //                                           commonNodes,
        //                                           tempPrefix,
        //                                           parent1,
        //                                           parent2,
        //                                           tempNode,
        //                                           currentNode);
        //         } else {
        //             waitForInput("common node", DEBUG1);
        //
        //             tempPrefix.push_back(tempNode);
        //             setOfPaths.push_back(tempPrefix);
        //
        //             // if-Bedingung prüfen
        //             Rcpp::Rcout << "erste Bed: " << (std::find(currentNode.connectedNodes_.begin(), currentNode.connectedNodes_.end(), startNode) == currentNode.connectedNodes_.end())
        //                         << " / " << (currentNode.node_ != tempNode) << std::endl;
        //             Rcpp::Rcout << "current / start: " << problemData_.nodeMap_[currentNode.node_].id_ << " " << problemData_.nodeMap_[startNode].id_ << std::endl;
        //
        //             if (std::find(currentNode.connectedNodes_.begin(), currentNode.connectedNodes_.end(), startNode) == currentNode.connectedNodes_.end()
        //                     && currentNode.node_ != tempNode){
        //
        //                 std::stringstream ss;
        //                 ss << "add node " << problemData_.nodeMap_[startNode].id_ << " to "
        //                    << problemData_.nodeMap_[currentNode.node_].id_ << std::endl;
        //                 Rcpp::Rcout << ss.str();
        //                 waitForInput(ss.str(), DEBUG1);
        //                 currentNode.connectedNodes_.push_back(tempNode);
        //                 currentNode.degree_++;
        //             }
        //         }
        //     }
        // }
    }

    // equivalent to the addOperator, but for the special case if solution is empty
    void greedyAddNode(std::vector<MyGraph::Node>& solution){
        std::vector<MyGraph::Node> tempBestSolution;
        double bestValue = 0;
        Rcpp::Rcout << "Info: greedyAddNode was called because addOperator started with empty vector!" << std::endl;

        /*
         * this do-while loop is necessary for the case where no valid solution exists
         * where the for-loop would end with tempBestSolution being empty (which is not compatible
         * with the add-operator that follows afterwards).
         * In this case the search has to be repeated until at least 1 valid node is found 
         * or until the termination criterion is satisfied. In the latter case, it is likely that
         * no valid solution exists.
         */
        do{
            for (MyGraph::Node n : problemData_.destinations_) {
                solution.clear();
                solution.push_back(n);
                additionalLogData_.numberOfTestedBitVectors_++;
                ResultData res = evaluateSolutionMatrix(problemData_, solution, "value");
                if (res.length_ <= problemData_.budget_){
                    if (res.value_ / (2*res.length_) > bestValue) {
                        tempBestSolution.assign(solution.begin(), solution.end());
                        bestValue = res.value_;
                    }
                }
            }
        } while (tempBestSolution.empty() && !terminationCriterionSatisfied());
        
        if (tempBestSolution.empty()){
            Rcpp::Rcout << "Final budget: " << problemData_.budget_ << "\n";
            Rcpp::Rcout << additionalLogData_ << "\n";
            Rcpp::stop("Termination criterion satisfied, but no valid solution was found with greedyAddNode() and addOperator(). It is possible that no valid solution exists.");
        } else {
            
            // Rcpp::Rcout << "tempBestSolution: " << evaluateSolutionMatrix(problemData_, tempBestSolution) << "\n";
        }
        solution.assign(tempBestSolution.begin(), tempBestSolution.end());
    }

    void addOperator(std::vector<MyGraph::Node>& solution){
        bool insertionSuccessful = false;
        bool isFirstIteration = true;

        if (solution.empty()){
            greedyAddNode(solution);
        }

        std::vector<MyGraph::Node> unvisitedNodes;
        for (MyGraph::Node n : problemData_.destinations_){
            if (std::find(solution.begin(), solution.end(), n) == solution.end()){
                unvisitedNodes.push_back(n);
            }
        }

        // calculate index vectors
        std::vector<int> unvisitedIndices = nodeVectorToIndexVector(problemData_, unvisitedNodes);
        std::vector<int> solutionIndices = nodeVectorToIndexVector(problemData_, solution);


        // Rcpp::Rcout << "(all nodes: " << problemData_.destinations_.size() << ", unvisited: " << unvisitedIndices.size() << std::endl;
        // printVector(unvisitedIndices);
        // printNodeIdsOfVector(problemData_, problemData_.destinations_);
        // Rcpp::Rcout << "solution (visited) (indices, then IDs): ";
        // printVector(solutionIndices);
        // printNodeIdsOfVector(problemData_, solution);


        // add start node so that at least 3 closest nodes can be found
        solutionIndices.insert(solutionIndices.begin(), 0);
        // solutionIndices.push_back(0);

        // drei nächste Nachbarn finden
        std::vector<ThreeNearest> threeNearestContainer;

        int lastInsertedNodeSolutionIndex = -1;
        int lastInsertedNodeDestinationIndex = -1;

        ThreeNearest *tnMax;

        waitForInput("add: before do-loop", DEBUG_ADD);
        do {
            double maxValue = 0;


            if (isFirstIteration){
                // Rcpp::Rcout << "solutionIndices: ";
                // printVector(solutionIndices);
                waitForInput("add: isFirstIteration", DEBUG_ADD);
                for (int i = 0; i < (int) unvisitedIndices.size(); i++){
                    int iu = unvisitedIndices[i];
                    std::vector<NodeIndexWithDistance> currentDistances;
                    for (int j = 0; j < (int) solutionIndices.size(); j++){
                        //Rcpp::Rcout << "solInd - 1: " << solutionIndices[j]-1 << std::endl;
                        int idToWrite = (solutionIndices[j] == 0) ? problemData_.nodeMap_[problemData_.startNode_].id_ : problemData_.nodeMap_[problemData_.destinations_[solutionIndices[j]-1]].id_;
                        currentDistances.push_back(NodeIndexWithDistance(solutionIndices[j],
                                                                         j,
                                                                         idToWrite,
                                                                         getAndLogDistance(iu, solutionIndices[j]))
                        );
                    }
                    std::sort(currentDistances.begin(), currentDistances.end(), compareNodeWithDistance);

                    // is it possible that there exist no three closest neighbors?
                    // this should not happen if solution contains at least 1 node, then the path would be (path, node, path)
                    // such that at least three possible neighbors can be found, I think...

                    int idToWrite = (iu == 0) ? problemData_.nodeMap_[problemData_.startNode_].id_ : problemData_.nodeMap_[problemData_.destinations_[iu-1]].id_;

                    if (iu == 0){
                        Rcpp::Rcerr << "in add operator (before creation of tempTn): iu has value 0" << std::endl;
                    }

                    // if the solution is (start, node, start), i.e. only 2 different nodes then add the start node
                    if (currentDistances.size() < 3){
                        currentDistances.push_back(NodeIndexWithDistance(-11, -11, -11, -11));
                    }
                    ThreeNearest tempTn(iu,
                                        i,
                                        idToWrite,
                                        currentDistances[0],
                                                        currentDistances[1],
                                                                        currentDistances[2]);

                    if (tempTn.neigh3_.destinationIndex_ == -11){
                        tempTn.numberOfNeighborsFound_ = 2;
                    } else {
                        tempTn.numberOfNeighborsFound_ = 3;
                    }

                    // Rcpp::Rcout << tempTn << std::endl;

                    threeNearestContainer.push_back(tempTn);//ThreeNearest(iu, i, currentDistances[0], currentDistances[1], currentDistances[2]));
                }
            } else {
                waitForInput("add: is not firstIteration", DEBUG_ADD);
                // check, if the newest inserted node is closer than the three previously saved nodes
                for (ThreeNearest &tn : threeNearestContainer){

                    // update positions since all nodes after the newly inserted node move to the right
                    if (tn.neigh1_.indexInSolution_ >= lastInsertedNodeSolutionIndex){
                        tn.neigh1_.indexInSolution_++;
                    }
                    if (tn.neigh2_.indexInSolution_ >= lastInsertedNodeSolutionIndex){
                        tn.neigh2_.indexInSolution_++;
                    }
                    if (tn.neigh3_.indexInSolution_ >= lastInsertedNodeSolutionIndex){
                        tn.neigh3_.indexInSolution_++;
                    }

                    if (tn.numberOfNeighborsFound_ <= 2) {
                        // if only two closest nodes have been found in the last iteration,
                        // the third (dummy) node is overwritten
                        double newDist = getAndLogDistance(lastInsertedNodeDestinationIndex, tn.thisNodeDestinationIndex_);

                        NodeIndexWithDistance newNeighbor(lastInsertedNodeDestinationIndex,
                                                          lastInsertedNodeSolutionIndex,
                                                          problemData_.nodeMap_[problemData_.destinations_[lastInsertedNodeDestinationIndex-1]].id_,
                                                          newDist);
                        tn.neigh3_ = newNeighbor;
                        tn.numberOfNeighborsFound_ = 3;
                    } else {

                        double maxDist = tn.neigh1_.distance_;

                        //int farthestNeighbor = 1;
                        NodeIndexWithDistance* farthestNeighbor = &tn.neigh1_;

                        // find the neighbor that has to be removed if the newly inserted node is closer to the current node

                        if (tn.neigh2_.distance_ > maxDist){
                            maxDist = tn.neigh2_.distance_;
                            //farthestNeighbor = 2;
                            farthestNeighbor = &tn.neigh2_;
                        }
                        if (tn.neigh3_.distance_ > maxDist){
                            maxDist = tn.neigh3_.distance_;
                            //farthestNeighbor = 3;
                            farthestNeighbor = &tn.neigh3_;
                        }
                        double newDist = getAndLogDistance(lastInsertedNodeDestinationIndex, tn.thisNodeDestinationIndex_);

                        NodeIndexWithDistance newNeighbor(lastInsertedNodeDestinationIndex,
                                                          lastInsertedNodeSolutionIndex,
                                                          problemData_.nodeMap_[problemData_.destinations_[lastInsertedNodeDestinationIndex-1]].id_,
                                                          newDist);

                        // Rcpp::Rcout << "maxDist: " << maxDist << std::endl;
                        // Rcpp::Rcout << "for: " << tn << std::endl;
                        // Rcpp::Rcout << "newly inserted node: " << newNeighbor << std::endl;
                        // Rcpp::Rcout << "farthestNeighbor: " << *farthestNeighbor << std::endl;
                        if (newDist < maxDist){
                            NodeIndexWithDistance newNeighbor(lastInsertedNodeDestinationIndex,
                                                              lastInsertedNodeSolutionIndex,
                                                              problemData_.nodeMap_[problemData_.destinations_[lastInsertedNodeDestinationIndex-1]].id_,
                                                              newDist);
                            // Rcpp::Rcout << "minNeigh replaced" << std::endl;
                            *farthestNeighbor = newNeighbor;

                        }
                    }
                }
            }



            // check the type of adjacency
            waitForInput("add: check type of adjacency", DEBUG_ADD);
            for (ThreeNearest &tn : threeNearestContainer) {
                threeNearestCalculateMinCost(tn, solutionIndices);

                // Rcpp::Rcout << tn << std::endl;
                waitForInput("add: check the length of the resulting solution", DEBUG_ADD);

                // check the length of the resulting solution
                std::vector<MyGraph::Node> tempSolution(solution.begin(), solution.end());
                // Rcpp::Rcout << "solution: ";
                // printNodeIdsOfVector(problemData_, tempSolution);
                ResultData before = evaluateSolutionMatrix(problemData_, tempSolution);
                // Rcpp::Rcout << "before: " << before << std::endl;

                MyGraph::Node nodeToInsert = problemData_.destinations_[tn.thisNodeDestinationIndex_-1];
                tempSolution.insert(tempSolution.begin() + tn.minInsertPosition - 1, nodeToInsert);

                additionalLogData_.numberOfTestedBitVectors_++;
                // printNodeIdsOfVector(problemData_, tempSolution);
                ResultData after = evaluateSolutionMatrix(problemData_, tempSolution);
                // Rcpp::Rcout << "after: " << after << std::endl;
                waitForInput("add: check the length of the resulting solution", DEBUG_ADD);

                double length = evaluateSolutionMatrix(problemData_, tempSolution).length_;
                if (length <= problemData_.budget_){
                    if (tn.addCost_ == 0){
                        // LARGE_NUMBER is used here so that multiple valid nodes with addCost=0 can be compared
                        tn.addValue_ = problemData_.nodeMap_[problemData_.destinations_[tn.thisNodeDestinationIndex_ - 1]].value_ * LARGE_NUMBER;
                    } else {
                        tn.addValue_ = problemData_.nodeMap_[problemData_.destinations_[tn.thisNodeDestinationIndex_ - 1]].value_ / tn.addCost_;
                    }
                    // Rcpp::Rcout << "-- Candidate for tnMax? length " << length << " / " << problemData_.budget_
                    // << " - value: " << tn.addValue_ << std::endl;
                    if (maxValue < tn.addValue_){
                        maxValue = tn.addValue_;
                        tnMax = &tn;
                        // Rcpp::Rcout << "tnMax updated. Now: " << *tnMax << std::endl;
                        insertionSuccessful = true;
                    }
                } else {
                    tn.addValue_ = 0;
                }
            }

            waitForInput("add: insert into solution", DEBUG_ADD);
            if (maxValue > 0){
                MyGraph::Node nodeToInsert = problemData_.destinations_[tnMax->thisNodeDestinationIndex_ - 1];
                solution.insert(solution.begin() + tnMax->minInsertPosition - 1, nodeToInsert);
                solutionIndices.insert(solutionIndices.begin() + tnMax->minInsertPosition, tnMax->thisNodeDestinationIndex_);
                lastInsertedNodeSolutionIndex = tnMax->minInsertPosition;
                lastInsertedNodeDestinationIndex = tnMax->thisNodeDestinationIndex_;


                /*
                * correct thisNodeIndexInUnvisited_ for all other TNs
                */
                for (ThreeNearest& tn : threeNearestContainer){
                    if (tn.thisNodeIndexInUnvisited_ > tnMax->thisNodeIndexInUnvisited_){
                        tn.thisNodeIndexInUnvisited_--;
                    }
                }

                unvisitedIndices.erase(unvisitedIndices.begin() + tnMax->thisNodeIndexInUnvisited_);
                unvisitedNodes.erase(unvisitedNodes.begin() + tnMax->thisNodeIndexInUnvisited_);
                threeNearestContainer.erase(threeNearestContainer.begin() + tnMax->thisNodeIndexInUnvisited_);
            } else {
                insertionSuccessful = false;
            }

            isFirstIteration = false;
            // Rcpp::Rcout << "new solution: ";
            // printNodeIdsOfVector(problemData_, solution);
            // Rcpp::Rcout << "new solution indices: ";
            // printVector(solutionIndices);
            ResultData after = evaluateSolutionMatrix(problemData_, solution);
            // Rcpp::Rcout << "after: " << after << std::endl;
        } while (insertionSuccessful);
    }

    // find, what type of adjacency there is and calculate minCost
    void threeNearestCalculateMinCost(ThreeNearest& tn, const std::vector<int>& solutionIndices){
        int numberOfDifferencesEqualTo1 = 0;

        // for the cases where the best insertion is at the very end, i.e. the start node is actually the end node
        bool diff12_reverse1 = false;
        bool diff12_reverse2 = false;
        bool diff13_reverse1 = false;
        bool diff13_reverse3 = false;
        bool diff23_reverse2 = false;
        bool diff23_reverse3 = false;

        int diff = tn.neigh1_.indexInSolution_ - tn.neigh2_.indexInSolution_;
        bool diff12Equal1, diff13Equal1, diff23Equal1;
        diff12Equal1 = (std::abs(diff) == 1);
        // Rcpp::Rcout << "diff1-2: " << diff << std::endl;

        // ThreeNearest backUpTN = tn;

        // edge case
        if ( (tn.neigh1_.indexInSolution_ == 0 && tn.neigh2_.indexInSolution_ == (int) solutionIndices.size()-1) ||
             (tn.neigh1_.indexInSolution_ == (int) solutionIndices.size()-1 && tn.neigh2_.indexInSolution_ == 0) ) {
            // Rcpp::Rcout << "diff1-2 reverse true" << std::endl;
            if (tn.neigh1_.indexInSolution_ == 0) {
                diff12_reverse1 = true;
                // tn.neigh1_.indexInSolution_ = solutionIndices.size();
            } else {
                diff12_reverse2 = true;
                // tn.neigh2_.indexInSolution_ = solutionIndices.size();
            }
            diff12Equal1 = true;
        }

        diff = tn.neigh1_.indexInSolution_ - tn.neigh3_.indexInSolution_;
        diff13Equal1 = (std::abs(diff) == 1);
        // Rcpp::Rcout << "diff1-3: " << diff << std::endl;

        // edge case
        if ( (tn.neigh1_.indexInSolution_ == 0 && tn.neigh3_.indexInSolution_ == (int) solutionIndices.size()-1) ||
             (tn.neigh1_.indexInSolution_ == (int) solutionIndices.size()-1 && tn.neigh3_.indexInSolution_ == 0) ) {
            // Rcpp::Rcout << "diff1-3 reverse true" << std::endl;
            if (tn.neigh1_.indexInSolution_ == 0) {
                diff13_reverse1 = true;
                // tn.neigh1_.indexInSolution_ = solutionIndices.size();
            } else {
                diff13_reverse3 = true;
                // tn.neigh3_.indexInSolution_ = solutionIndices.size();
            }
            diff13Equal1 = true;
        }

        diff = tn.neigh2_.indexInSolution_ - tn.neigh3_.indexInSolution_;
        diff23Equal1 = (std::abs(diff) == 1);
        // Rcpp::Rcout << "diff2-3: " << diff << std::endl;

        // edge case
        if ( (tn.neigh2_.indexInSolution_ == 0 && tn.neigh3_.indexInSolution_ == (int) solutionIndices.size()-1) ||
             (tn.neigh2_.indexInSolution_ == (int) solutionIndices.size()-1 && tn.neigh3_.indexInSolution_ == 0) ) {
            // Rcpp::Rcout << "diff2-3 reverse true" << std::endl;
            if (tn.neigh2_.indexInSolution_ == 0) {
                diff23_reverse2 = true;
                // tn.neigh2_.indexInSolution_ = solutionIndices.size();
            } else {
                diff23_reverse3 = true;
                // tn.neigh3_.indexInSolution_ = solutionIndices.size();
            }
            diff23Equal1 = true;
        }

        numberOfDifferencesEqualTo1 = diff12Equal1 + diff13Equal1 + diff23Equal1;

        if (DEBUG_ADD) {
            Rcpp::Rcout << "solutionIndices: ";
            printVector(solutionIndices);

            Rcpp::Rcout << "diffEqualTo1: " << numberOfDifferencesEqualTo1 << std::endl;
        }
        waitForInput("before numberOfDifferencesEqual1-Check", DEBUG_ADD);

        if (numberOfDifferencesEqualTo1 == 0){
            // no adjacency between the three nearest visited nodes
            double minCost = std::numeric_limits<double>::max();
            int minInsertionNeighbor = 0;
            bool minInsertionAfter = false;


            NodeIndexWithDistance neighbors[3] = {tn.neigh1_, tn.neigh2_, tn.neigh3_};

            for (int i = 0; i < tn.numberOfNeighborsFound_; i++) {
                NodeIndexWithDistance& currentNeighbor = neighbors[i];


                int destinationIndexOfPreviousNode = (currentNeighbor.indexInSolution_ > 0) ?
                solutionIndices[currentNeighbor.indexInSolution_-1] : (int) solutionIndices[solutionIndices.size()-1] ;

                int destinationIndexOfNextNode = (currentNeighbor.indexInSolution_ < (int) solutionIndices.size()-1) ?
                solutionIndices[currentNeighbor.indexInSolution_+1] : solutionIndices[0];

                /*
                * To-Do: Wenn links vom Startknoten eingesetzt wird, muss die Startposition angepasst werden
                * (minInsertionNeighbor)
                */

                double lengthChange = getAndLogDistance(destinationIndexOfPreviousNode, tn.thisNodeDestinationIndex_)
                    + getAndLogDistance(tn.thisNodeDestinationIndex_, currentNeighbor.destinationIndex_)
                    - getAndLogDistance(destinationIndexOfPreviousNode, currentNeighbor.destinationIndex_);

                    if (lengthChange < minCost) {
                        minCost = lengthChange;
                        minInsertionNeighbor = i;
                        minInsertionAfter = false;
                    }

                    lengthChange = getAndLogDistance(currentNeighbor.destinationIndex_, tn.thisNodeDestinationIndex_)
                        + getAndLogDistance(tn.thisNodeDestinationIndex_, destinationIndexOfNextNode)
                        - getAndLogDistance(currentNeighbor.destinationIndex_, destinationIndexOfNextNode);

                        if (lengthChange < minCost) {
                            minCost = lengthChange;
                            minInsertionNeighbor = i;
                            minInsertionAfter = true;
                        }
            }

            tn.addCost_ = minCost;
            // if the best insertion is after the neighbor node, then add minInsertionAfter=1 (true)
            tn.minInsertPosition = neighbors[minInsertionNeighbor].indexInSolution_ + minInsertionAfter;

            // if the best insertion is to the left of the start node, that is equivalent to inserting at the end
            if (tn.minInsertPosition <= 0){
                tn.minInsertPosition = solutionIndices.size();
            }
            // Rcpp::Rcout << "00 update minCost / minInsertPos: " << tn.addCost_ << " / " << tn.minInsertPosition << std::endl;
        }

        // if (numberOfDifferencesEqualTo1 >= 1){ // Is that case not equivalent to the case >= 2??
        //     // only 2 nodes are adjacent; then addCost is equal to the length increase after inserting
        //     // the new node between these 2 nodes
        //     if (diff12Equal1){
        //         tn.addCost_ = getAndLogDistance(tn.thisNodeDestinationIndex_, tn.neigh1_.destinationIndex_)
        //                     + getAndLogDistance(tn.neigh2_.destinationIndex_, tn.thisNodeDestinationIndex_)
        //                     - getAndLogDistance(tn.neigh1_.destinationIndex_, tn.neigh2_.destinationIndex_);
        //
        //         if (tn.neigh1_.indexInSolution_ > tn.neigh2_.indexInSolution_){
        //             tn.minInsertPosition = tn.neigh1_.indexInSolution_;
        //         } else {
        //             tn.minInsertPosition = tn.neigh2_.indexInSolution_;
        //         }
        //     }
        //     if (diff13Equal1){
        //         tn.addCost_ = getAndLogDistance(tn.thisNodeDestinationIndex_, tn.neigh1_.destinationIndex_)
        //                     + getAndLogDistance(tn.neigh3_.destinationIndex_, tn.thisNodeDestinationIndex_)
        //                     - getAndLogDistance(tn.neigh1_.destinationIndex_, tn.neigh3_.destinationIndex_);
        //
        //         if (tn.neigh1_.indexInSolution_ > tn.neigh3_.indexInSolution_){
        //             tn.minInsertPosition = tn.neigh1_.indexInSolution_;
        //         } else {
        //             tn.minInsertPosition = tn.neigh3_.indexInSolution_;
        //         }
        //     }
        //     if (diff23Equal1){
        //         tn.addCost_ = getAndLogDistance(tn.thisNodeDestinationIndex_, tn.neigh2_.destinationIndex_)
        //                     + getAndLogDistance(tn.neigh2_.destinationIndex_, tn.thisNodeDestinationIndex_)
        //                     - getAndLogDistance(tn.neigh2_.destinationIndex_, tn.neigh3_.destinationIndex_);
        //
        //         if (tn.neigh2_.indexInSolution_ > tn.neigh3_.indexInSolution_){
        //             tn.minInsertPosition = tn.neigh2_.indexInSolution_;
        //         } else {
        //             tn.minInsertPosition = tn.neigh3_.indexInSolution_;
        //         }
        //     }
        // }

        if (numberOfDifferencesEqualTo1 >= 1){
            // mindestens zwei Mal 1: alle
            tn.addCost_ = std::numeric_limits<double>::max();
            double minCost = std::numeric_limits<double>::max();
            if (diff12Equal1){
                if ((tn.neigh1_.indexInSolution_ > tn.neigh2_.indexInSolution_ && !diff12_reverse2) || diff12_reverse1){
                    double currentCost =  getAndLogDistance(tn.neigh2_.destinationIndex_, tn.thisNodeDestinationIndex_)
                    + getAndLogDistance(tn.thisNodeDestinationIndex_, tn.neigh1_.destinationIndex_)
                    - getAndLogDistance(tn.neigh2_.destinationIndex_, tn.neigh1_.destinationIndex_);

                    if (currentCost < minCost){

                        tn.addCost_ = currentCost;
                        minCost = currentCost;
                        tn.minInsertPosition = tn.neigh1_.indexInSolution_;
                        if (diff12_reverse1 || diff12_reverse2){
                            tn.minInsertPosition = solutionIndices.size();
                        }
                        // Rcpp::Rcout << "A update minCost / minInsertPos: " << tn.addCost_ << " / " << tn.minInsertPosition << std::endl;
                    }
                } else {
                    double currentCost =  getAndLogDistance(tn.neigh1_.destinationIndex_, tn.thisNodeDestinationIndex_)
                    + getAndLogDistance(tn.thisNodeDestinationIndex_, tn.neigh2_.destinationIndex_)
                    - getAndLogDistance(tn.neigh1_.destinationIndex_, tn.neigh2_.destinationIndex_);

                    if (currentCost < minCost){
                        tn.addCost_ = currentCost;
                        minCost = currentCost;
                        tn.minInsertPosition = tn.neigh2_.indexInSolution_;
                        if (diff12_reverse2 || diff12_reverse1){
                            tn.minInsertPosition = solutionIndices.size();
                        }
                        // Rcpp::Rcout << "B update minCost / minInsertPos: " << tn.addCost_ << " / " << tn.minInsertPosition << std::endl;
                    }
                }
            }
            if (diff13Equal1){
                if ((tn.neigh1_.indexInSolution_ > tn.neigh3_.indexInSolution_ && !diff13_reverse3) || diff13_reverse1){
                    double currentCost =  getAndLogDistance(tn.neigh3_.destinationIndex_, tn.thisNodeDestinationIndex_)
                    + getAndLogDistance(tn.thisNodeDestinationIndex_, tn.neigh1_.destinationIndex_)
                    - getAndLogDistance(tn.neigh3_.destinationIndex_, tn.neigh1_.destinationIndex_);

                    if (currentCost < minCost){
                        tn.addCost_ = currentCost;
                        minCost = currentCost;
                        tn.minInsertPosition = tn.neigh1_.indexInSolution_;
                        if (diff13_reverse1 || diff13_reverse3){
                            tn.minInsertPosition = solutionIndices.size();
                        }
                        // Rcpp::Rcout << "C update minCost / minInsertPos: " << tn.addCost_ << " / " << tn.minInsertPosition << std::endl;
                    }

                } else {
                    double currentCost =  getAndLogDistance(tn.neigh1_.destinationIndex_, tn.thisNodeDestinationIndex_)
                    + getAndLogDistance(tn.thisNodeDestinationIndex_, tn.neigh3_.destinationIndex_)
                    - getAndLogDistance(tn.neigh1_.destinationIndex_, tn.neigh3_.destinationIndex_);

                    if (currentCost < minCost){
                        tn.addCost_ = currentCost;
                        minCost = currentCost;
                        tn.minInsertPosition = tn.neigh3_.indexInSolution_;
                        if (diff13_reverse1 || diff13_reverse3){
                            tn.minInsertPosition = solutionIndices.size();
                        }
                        // Rcpp::Rcout << "D update minCost / minInsertPos: " << tn.addCost_ << " / " << tn.minInsertPosition << std::endl;
                    }
                }
            }
            if (diff23Equal1){
                if ((tn.neigh2_.indexInSolution_ > tn.neigh3_.indexInSolution_ && !diff23_reverse3) || diff23_reverse2){
                    double currentCost =  getAndLogDistance(tn.neigh3_.destinationIndex_, tn.thisNodeDestinationIndex_)
                    + getAndLogDistance(tn.thisNodeDestinationIndex_, tn.neigh2_.destinationIndex_)
                    - getAndLogDistance(tn.neigh3_.destinationIndex_, tn.neigh2_.destinationIndex_);

                    if (currentCost < minCost){
                        tn.addCost_ = currentCost;
                        minCost = currentCost;
                        tn.minInsertPosition = tn.neigh2_.indexInSolution_;
                        if (diff23_reverse2){
                            tn.minInsertPosition = solutionIndices.size();
                        }
                        // Rcpp::Rcout << "E update minCost / minInsertPos: " << tn.addCost_ << " / " << tn.minInsertPosition << std::endl;
                    }
                } else {
                    double currentCost =  getAndLogDistance(tn.neigh2_.destinationIndex_, tn.thisNodeDestinationIndex_)
                    + getAndLogDistance(tn.thisNodeDestinationIndex_, tn.neigh3_.destinationIndex_)
                    - getAndLogDistance(tn.neigh2_.destinationIndex_, tn.neigh3_.destinationIndex_);

                    if (currentCost < minCost){
                        tn.addCost_ = currentCost;
                        minCost = currentCost;
                        tn.minInsertPosition = tn.neigh3_.indexInSolution_;
                        if (diff23_reverse3){
                            tn.minInsertPosition = solutionIndices.size();
                        }
                        // Rcpp::Rcout << "F update minCost / minInsertPos: " << tn.addCost_ << " / " << tn.minInsertPosition << std::endl;
                    }
                }
            }
            // revert changes to neigh1,2,3
            // tn.neigh1_ = backUpTN.neigh1_;
            // tn.neigh2_ = backUpTN.neigh2_;
            // tn.neigh3_ = backUpTN.neigh3_;
        }
    }

    void dropOperator(std::vector<MyGraph::Node>& solution){
        double currentLength = evaluateSolutionMatrix(problemData_, solution).length_;

        while (currentLength > problemData_.budget_){
            std::vector<double> dropValues;

            // first, find the indices in the distance matrix corresponding to the nodes in the solution
            std::vector<int> indicesVector = nodeVectorToIndexVector(problemData_, solution);
            for (int i = 0; i < (int) solution.size(); i++){
                int prev = (i==0) ? 0 : indicesVector[i-1];
                int curr = indicesVector[i];
                int next = (i == (int) solution.size()-1) ? 0 : indicesVector[i+1];
                double lengthSaved = (getAndLogDistance(prev,curr) + getAndLogDistance(curr, next) - getAndLogDistance(prev,next));

                double dropValue = problemData_.nodeMap_[solution[i]].value_ / lengthSaved;

                additionalLogData_.numberOfTestedBitVectors_++;

                dropValues.push_back(dropValue);
            }

            // find minimal drop value (is sorting really necessary?)
            double minDrop = std::numeric_limits<double>::max();
            int minIndex = -1;
            for (int i = 0; i < (int) dropValues.size(); i++){
                if (dropValues[i] < minDrop){
                    minDrop = dropValues[i];
                    minIndex = i;
                }
            }
            
            if (minIndex == -1) {
                /* 
                 * This (rare) case occurs when for each nodes in the solution no decrease
                 * in path length occurs when it is removed.
                 * In this case: Remove a random node. It is possible that the path
                 * length decreases in the subsequent removals.
                 */
                Rcpp::Rcout << "Info: minIndex=-1 in dropOperator. Remove a different node...\n";
                printNodeIdsOfVector(problemData_, solution);
                Rcpp::Rcout << evaluateSolutionMatrix(problemData_, solution) << "\n";
                Rcpp::Rcout << "solutionSize: " << solution.size() << "\n";
                Rcpp::Rcout << "dropValues.size: " << dropValues.size() << "\n";
                for (int i = 0; i < (int) dropValues.size(); i++){
                    Rcpp::Rcout << dropValues[i] << " / " << minDrop << "\n";
                }
                minIndex = getRandomNumber(0, dropValues.size()-1);
                // Rcpp::stop("minIndex ist -1???");
            }
            solution.erase(solution.begin() + minIndex);

            currentLength = evaluateSolutionMatrix(problemData_, solution).length_;
        }
        // Rcpp::Rcout << "postDrop2: " << additionalLogData_.numberOfShortestPathCalls_ << ", bestValue: " << bestSolutionQuality_.value_ << "\n";
        // waitForInput("DE-output", true);
    }

    void resetSolver() override{
        if (initialSolutionBackup_.empty()) {
            initialSolutionBackup_.assign(initialSolution_.begin(), initialSolution_.end());
            initialSolutionQuality_ =  evaluateSolutionMatrix(problemData_, initialSolutionBackup_);
        }
        if (initialSolutionQuality_.length_ > problemData_.budget_){
            initialSolution_.clear();
            bestSolution_.clear();
            bestSolutionQuality_ = ResultData(0,0);
            Rcpp::Rcout << "Restart algorithm with empty initial solution.\n";
        } else {
            initialSolution_.assign(initialSolutionBackup_.begin(), initialSolutionBackup_.end());
            bestSolution_.assign(initialSolution_.begin(), initialSolution_.end());
            bestSolutionQuality_ = initialSolutionQuality_;
            Rcpp::Rcout << "Restart algorithm with the given initial solution.\n";
        }
        Rcpp::Rcout << "Current log data: " << additionalLogData_ << "\n";
        run(); 
        
        std::vector<MyGraph::Node> mySolution = asCompleteSolution(bestSolution_);
        writeSolution(mySolution, !calledAsImprover_);
        Rcpp::Rcout << "Log data at the end:" << additionalLogData_ << "\n";
        Rcpp::stop("Algorithm used restarts, but the termination criterion is now satisfied");
    }
    // double getAndLogDistance(int i, int j){
    //     additionalLogData_.numberOfShortestPathCalls_++;
    //     //Rcpp::Rcout << "(i,j): " << i << ", " << j << " - " << distanceMatrix_(i,j) << std::endl;
    //     return distanceMatrix_(i,j);
    // }

    // ResultData evaluateSolutionMatrix(ProblemData& problemData,
    //                                   const std::vector<MyGraph::Node>& solution,
    //                                   std::string targetCriterion = "value",
    //                                   bool forceLogging = false,
    //                                   bool abortOnInvalidity = true){
    //     return evaluateSolution(problemData_, solution, targetCriterion, forceLogging, abortOnInvalidity, &distanceMatrix_);
    // }

};


//' @export
// [[Rcpp::export]]
void callEa4OpSolver(const Rcpp::DataFrame& nodeDf,
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

    // for (MyGraph::NodeIt n(graph); n != lemon::INVALID; ++n){
    //     Rcpp::Rcout << nodeMap[n].coordinates_.x << "," << nodeMap[n].coordinates_.y << " - " << nodeMap[n].id_ << ", "
    //                 << nodeMap[n].type_ << ", " << nodeMap[n].value_ << std::endl;
    // }
    // for (MyGraph::ArcIt a(graph); a != lemon::INVALID; ++a){
    //     Rcpp::Rcout << arcMap[a].fromID_ << "," << arcMap[a].toID_ << " - " << arcMap[a].length_ << std::endl;
    // }




    ProblemData pd(graph, nodeMap, arcMap, startNode, destinations, budget, problemName);
    if (pathToChanges != "") {
        readDynamicChangesFromFile(pd, pathToChanges);
    }

    Rcpp::Environment myEnvironment(MYPACKAGE);
    Ea4OpSolver ea(pd, problemName, runNumber, myEnvironment, "value", fileSuffix, pathToDistanceMatrix);
    ea.run();


    /////////////////////////////

    std::vector<MyGraph::Node> mySolution = ea.asCompleteSolution(ea.bestSolution_, true, true);
    ea.writeSolution(mySolution, true);

}



//' @export
// [[Rcpp::export]]
void callEa4OpImprover(const Rcpp::DataFrame& nodeDf,
                       const Rcpp::DataFrame& arcDf,
                       const Rcpp::DataFrame& problemDf,
                       double budget,
                       std::string problemName,
                       unsigned int runNumber,
                       std::string pathToInitialSolution,
                       double dominanceOfInitialSolution = 0.5,
                       std::string fileSuffix = "",
                       std::string pathToChanges = "",
                       std::string pathToDistanceMatrix = "",
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
     * EA4OP does not have an initial solution, but a initial population.
     * How to incorporate the initial solution?
     * Idea: Initialize half of the population using mutated variants of the initial solution.
     */

    Rcpp::Environment myEnvironment(MYPACKAGE);
    Ea4OpSolver ea(pd, problemName, runNumber, myEnvironment, "value", fileSuffix, pathToDistanceMatrix);

    ea.typeOfHandling_ = budgetChangeHandlingMode;
    if (budgetChangeHandlingMode == Solver::HANDLING_TABLE) {
        ea.initializeSampledTable(minBudgetToHandle, maxBudgetToHandle, 
                                   budgetChangeTableSize);
    }
    
    ea.readInitialSolutionFromFile(pathToInitialSolution);
    ea.calledAsImprover_ = true;
    // waitForInput("beforeRun", DEBUG_ENABLED);
    ea.dominanceOfInitialSolution_ = dominanceOfInitialSolution;

    // if (pathToDistanceMatrix != "") {
    //     ea.readDi
    // } else {
    //
    // }
    ea.run();

    /////////////////////////////

    std::vector<MyGraph::Node> mySolution = ea.asCompleteSolution(ea.bestSolution_, true, false);
    ea.writeSolution(mySolution, false);

}
