//
//  GeneticAlgorithm.hpp
//  AKML Project
//
//  Created by Aldric Labarthe on 13/09/2023.
//
#ifndef GeneticAlgorithm_hpp
#define GeneticAlgorithm_hpp

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <functional>
#include <array>
#include <algorithm>
#include <random>
#include <string>
#include <utility>

#include "AKML_consts.hpp"
#include "UtilityLoops.hpp"
#include "NeuralNetwork.hpp"
#include "NeuralLayer.hpp"
#include "Matrices.hpp"

namespace akml {

class GeneticAlgorithmMethods {
public:

    /* @deprecated rather use the pointer version for child. */
    template <typename NEURAL_NET_TYPE, typename NEURAL_LAYER_TYPE>
    static void mergeLayers(const std::size_t layer, NEURAL_NET_TYPE* parent1, NEURAL_NET_TYPE* parent2, NEURAL_NET_TYPE& child){
        NEURAL_LAYER_TYPE* parent1_layer = (NEURAL_LAYER_TYPE*)parent1->getLayer(layer);
        NEURAL_LAYER_TYPE* parent2_layer = (NEURAL_LAYER_TYPE*)parent2->getLayer(layer);
        NEURAL_LAYER_TYPE* child_layer = (NEURAL_LAYER_TYPE*)child.getLayer(layer);
        
        if (parent2_layer == parent1_layer){
            *(child_layer->getWeightsAccess()) = *(parent2_layer->getWeightsAccess());
            *(child_layer->getBiasesAccess()) = *(parent2_layer->getBiasesAccess());
            return;
        }
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0,0.15);
        
        // Biases :
        for (std::size_t row(0); row < child_layer->getNeuronNumber(); row++){
            if (row %2 == 0)
                child_layer->getBiasesAccess()->operator()(row+1, 1) = parent1_layer->getBiasesAccess()->operator()(row+1, 1);
            else
                child_layer->getBiasesAccess()->operator()(row+1, 1) = parent2_layer->getBiasesAccess()->operator()(row+1, 1);
            child_layer->getBiasesAccess()->operator()(row+1, 1) = child_layer->getBiasesAccess()->operator()(row+1, 1) * (1+distribution(gen));
            
        }
        
        // weights :
        for (std::size_t row(0); row < child_layer->getNeuronNumber(); row++){
            for (std::size_t col(0); col < child_layer->getPreviousNeuronNumber(); col++){
                if ((row+col) %2 == 0)
                    child_layer->getWeightsAccess()->operator()(row+1, col+1) = parent1_layer->getWeightsAccess()->operator()(row+1, col+1);
                else
                    child_layer->getWeightsAccess()->operator()(row+1, col+1) = parent2_layer->getWeightsAccess()->operator()(row+1, col+1);
                
                child_layer->getWeightsAccess()->operator()(row+1, col+1) = child_layer->getWeightsAccess()->operator()(row+1, col+1) * (1+distribution(gen));
            }
        }
    }
    
    template <typename NEURAL_NET_TYPE, typename NEURAL_LAYER_TYPE>
    static void mergeLayers(const std::size_t layer, NEURAL_NET_TYPE* parent1, NEURAL_NET_TYPE* parent2, NEURAL_NET_TYPE* child){
        NEURAL_LAYER_TYPE* parent1_layer = (NEURAL_LAYER_TYPE*)parent1->getLayer(layer);
        NEURAL_LAYER_TYPE* parent2_layer = (NEURAL_LAYER_TYPE*)parent2->getLayer(layer);
        NEURAL_LAYER_TYPE* child_layer = (NEURAL_LAYER_TYPE*)child->getLayer(layer);
        
        if (parent2_layer == parent1_layer){
            *(child_layer->getWeightsAccess()) = *(parent2_layer->getWeightsAccess());
            *(child_layer->getBiasesAccess()) = *(parent2_layer->getBiasesAccess());
            return;
        }
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<double> distribution(0.0,0.15);
        
        // Biases :
        for (std::size_t row(0); row < child_layer->getNeuronNumber(); row++){
            if (row %2 == 0)
                child_layer->getBiasesAccess()->operator()(row+1, 1) = parent1_layer->getBiasesAccess()->operator()(row+1, 1);
            else
                child_layer->getBiasesAccess()->operator()(row+1, 1) = parent2_layer->getBiasesAccess()->operator()(row+1, 1);
            child_layer->getBiasesAccess()->operator()(row+1, 1) = child_layer->getBiasesAccess()->operator()(row+1, 1) * (1+distribution(gen));
            
        }
        
        // weights :
        for (std::size_t row(0); row < child_layer->getNeuronNumber(); row++){
            for (std::size_t col(0); col < child_layer->getPreviousNeuronNumber(); col++){
                if ((row+col) %2 == 0)
                    child_layer->getWeightsAccess()->operator()(row+1, col+1) = parent1_layer->getWeightsAccess()->operator()(row+1, col+1);
                else
                    child_layer->getWeightsAccess()->operator()(row+1, col+1) = parent2_layer->getWeightsAccess()->operator()(row+1, col+1);
                
                child_layer->getWeightsAccess()->operator()(row+1, col+1) = child_layer->getWeightsAccess()->operator()(row+1, col+1) * (1+distribution(gen));
            }
        }
    }
};

template <size_t NBLAYERS, size_t INPUTNUMBER, size_t OUTPUTNUMBER, size_t POP_SIZE=50>
class BaseGeneticAlgorithm : public GeneticAlgorithmMethods {
protected:
    std::function<void(NeuralNetwork<NBLAYERS>&)> NNInstructions;
    std::array<NeuralNetwork<NBLAYERS>*, POP_SIZE> networksPopulation;
    std::function<StaticMatrix <float, OUTPUTNUMBER, 1>(StaticMatrix <float, OUTPUTNUMBER, 1>)>* postactivation_process = nullptr;
    std::function<void(NeuralNetwork<NBLAYERS>*, NeuralNetwork<NBLAYERS>*)>* debug_printer = nullptr;
    
    template<int i, int j>
    struct mergeLayer_functor {
        static void run(akml::NeuralNetwork<NBLAYERS>* child, akml::NeuralNetwork<NBLAYERS>* parent1, akml::NeuralNetwork<NBLAYERS>* parent2) {
            GeneticAlgorithmMethods::mergeLayers<akml::NeuralNetwork<NBLAYERS>, akml::NeuralLayer<akml::NN_STRUCTURE[j], akml::NN_STRUCTURE[i]>>(i+2, parent1, parent2, child);
        }
    };
    
    double newcomers_rate = (double)1 / 12;
    int newcomers_bound;
    double reproduction_rate = (double)7 / 8 - newcomers_rate;
    int reproduction_bound;
    double conservation_rate = 1 - newcomers_rate - reproduction_rate;
    int conservation_bound = POP_SIZE;
    
public:
    static std::function<StaticMatrix <float, OUTPUTNUMBER, 1>(StaticMatrix <float, OUTPUTNUMBER, 1>)> ACTIVATE_ROUND;
    
    static inline std::function<void(NeuralNetwork<NBLAYERS>&, NeuralNetwork<NBLAYERS>*, NeuralNetwork<NBLAYERS>*)> DEFAULT_MERGING_INSTRUCTIONS = [](akml::NeuralNetwork<NBLAYERS>& child, akml::NeuralNetwork<NBLAYERS>* parent1, akml::NeuralNetwork<NBLAYERS>* parent2) {
        if (parent1 == parent2)
            akml::GeneticAlgorithmMethods::mergeLayers<akml::NeuralNetwork<NBLAYERS>, akml::NeuralLayer<akml::NN_STRUCTURE[0], 1>>(1, parent1, parent2, child);
        akml::for_<0, NBLAYERS-1>::template run<mergeLayer_functor, NeuralNetwork<NBLAYERS>*, NeuralNetwork<NBLAYERS>*, NeuralNetwork<NBLAYERS>*>(&child, parent1, parent2);
    };
    
    inline BaseGeneticAlgorithm(const std::function<void(NeuralNetwork<NBLAYERS>&)> instructions = NeuralNetwork<NBLAYERS>::DEFAULT_INIT_INSTRUCTIONS) : NNInstructions(instructions){
        for (std::size_t i(0); i < POP_SIZE; i++){
            networksPopulation[i] = new NeuralNetwork<NBLAYERS> (std::string("RANDOM_INIT"));
            NNInstructions(*networksPopulation[i]);
        }
        
        newcomers_bound = std::round(newcomers_rate * POP_SIZE);
        reproduction_bound = std::round((reproduction_rate + newcomers_rate) * POP_SIZE);
    }
    
    inline ~BaseGeneticAlgorithm() {
        for (std::size_t i(0); i < POP_SIZE; i++){
            if (networksPopulation[i] != nullptr)
                delete networksPopulation[i];
        }
    };
    
    inline void setNewGenerationRates(double newcomers_r, double reproduction_r){
        reproduction_rate = reproduction_r;
        newcomers_rate = newcomers_r;
        conservation_rate = 1 - newcomers_rate - reproduction_rate;
        newcomers_bound = std::max(std::round(newcomers_rate * POP_SIZE), (double)1);
        reproduction_bound = std::round((reproduction_rate + newcomers_rate) * POP_SIZE);
    }
    
    inline void setPostActivationProcess(std::function<StaticMatrix <float, OUTPUTNUMBER, 1>(StaticMatrix <float, OUTPUTNUMBER, 1>)> func){
        postactivation_process = &func;
    }
    
    inline void setDebugFunction(std::function<void(NeuralNetwork<NBLAYERS>*, NeuralNetwork<NBLAYERS>*)> func){
        debug_printer = &func;
    }
    
    inline std::array<NeuralNetwork<NBLAYERS>*, POP_SIZE> getNetworksPopulation() {
        return networksPopulation;
    }
    
    inline std::array< StaticMatrix <float, OUTPUTNUMBER, 1>, POP_SIZE> evaluateNN(StaticMatrix <float, INPUTNUMBER, 1>& input){
        std::array< StaticMatrix <float, OUTPUTNUMBER, 1>, POP_SIZE> localoutputs;
        for (std::size_t netid(0); netid < POP_SIZE; netid++){
            localoutputs[netid] = *(networksPopulation[netid]->template process<INPUTNUMBER, OUTPUTNUMBER>(input));
            if (this->postactivation_process != nullptr)
                localoutputs[netid] = postactivation_process->operator()(localoutputs[netid]);
        }
        return localoutputs;
    };
    
    inline StaticMatrix <float, OUTPUTNUMBER, 1> evaluateOneNN(std::size_t netid, StaticMatrix <float, INPUTNUMBER, 1>& input){
        StaticMatrix <float, OUTPUTNUMBER, 1> localoutput = *(networksPopulation[netid]->template process<INPUTNUMBER, OUTPUTNUMBER>(input));
        if (this->postactivation_process != nullptr)
            localoutput = postactivation_process->operator()(localoutput);
        return localoutput;
    };

    
    /*
     MSE should be sorted from worst to best.
     */
    inline void generateNewGeneration(const std::array< std::pair<float, NeuralNetwork<NBLAYERS>*>, POP_SIZE>& MSE, std::size_t iteration,
                                      const std::function<void(NeuralNetwork<NBLAYERS>&, NeuralNetwork<NBLAYERS>*, NeuralNetwork<NBLAYERS>*)> merging_instructions = BaseGeneticAlgorithm<NBLAYERS, INPUTNUMBER, OUTPUTNUMBER, POP_SIZE>::DEFAULT_MERGING_INSTRUCTIONS) {
        std::array<NeuralNetwork<NBLAYERS>*, POP_SIZE> newNetworksPopulation;
        
        float mean_MSE(0);
        for (std::size_t i(0); i < newcomers_bound; i++){
            newNetworksPopulation[i] = new NeuralNetwork<NBLAYERS> (std::string("RANDOM_") + std::to_string(iteration));
            NNInstructions(*newNetworksPopulation[i]);
            mean_MSE += MSE[i].first;
        }
        for (std::size_t i(newcomers_bound); i < reproduction_bound; i++){
            newNetworksPopulation[i] = new NeuralNetwork<NBLAYERS> (std::string("CHILD_") + std::to_string(iteration));
            NNInstructions(*newNetworksPopulation[i]);
            merging_instructions(*newNetworksPopulation[i],
                                 MSE[POP_SIZE-i-1].second,
                                 MSE[POP_SIZE-i].second);
            mean_MSE += MSE[i].first;
        }
        for (std::size_t i(reproduction_bound); i < POP_SIZE; i++){
            newNetworksPopulation[i] = new NeuralNetwork<NBLAYERS> (MSE[i].second->getCustomOriginField());
            NNInstructions(*newNetworksPopulation[i]);
            merging_instructions(*newNetworksPopulation[i], MSE[i].second, MSE[i].second);
            mean_MSE += MSE[i].first;
            
            if (debug_printer != nullptr)
                debug_printer->operator()(MSE[i].second, newNetworksPopulation[i]);
        }
        mean_MSE = mean_MSE/POP_SIZE;
        std::cout << "Best MSE = " << MSE[POP_SIZE-1].first << " - Mean=" << mean_MSE << " - Origin: " << MSE[POP_SIZE-1].second->getCustomOriginField() << std::endl;
        for (std::size_t i(0); i < POP_SIZE; i++){
            delete networksPopulation[i];
            networksPopulation[i] = newNetworksPopulation[i];
        }
    }
};

template <size_t NBLAYERS, size_t INPUTNUMBER, size_t OUTPUTNUMBER, size_t POP_SIZE=50>
class ManualGeneticAlgorithm : public BaseGeneticAlgorithm<NBLAYERS, INPUTNUMBER, POP_SIZE>{
    
};

template <size_t NBLAYERS, size_t INPUTNUMBER, size_t OUTPUTNUMBER, size_t TRAINING_LENGTH, size_t POP_SIZE=50>
class GeneticAlgorithm : public BaseGeneticAlgorithm<NBLAYERS, INPUTNUMBER, OUTPUTNUMBER, POP_SIZE> {
protected:
    std::array<StaticMatrix <float, INPUTNUMBER, 1>, TRAINING_LENGTH> inputs;
    std::array<StaticMatrix <float, OUTPUTNUMBER, 1>, TRAINING_LENGTH> outputs;
    
    
public:
    static std::function<float(std::array<StaticMatrix <float, OUTPUTNUMBER, 1>, TRAINING_LENGTH>&, std::array< StaticMatrix <float, OUTPUTNUMBER, 1> , TRAINING_LENGTH>&)> MSE;
    static std::function<float(std::array<StaticMatrix <float, OUTPUTNUMBER, 1>, TRAINING_LENGTH>&, std::array< StaticMatrix <float, OUTPUTNUMBER, 1> , TRAINING_LENGTH>&)> ERRORS_COUNT;
    
    
    inline GeneticAlgorithm(const std::array<StaticMatrix <float, INPUTNUMBER, 1>, TRAINING_LENGTH> in, const std::array<StaticMatrix <float, OUTPUTNUMBER, 1>, TRAINING_LENGTH> out, const std::function<void(NeuralNetwork<NBLAYERS>&)> instructions = NeuralNetwork<NBLAYERS>::DEFAULT_INIT_INSTRUCTIONS) : BaseGeneticAlgorithm<NBLAYERS, INPUTNUMBER, OUTPUTNUMBER, POP_SIZE>(instructions), inputs(in), outputs(out) {
    };
    
    /*
     Choices about generation of new population :
     Flop 1/12 is replaced by new networks generated randomly (introduction of immigration)
     Top (7/8-1/12) is selected and multiplied so that their offspring represents less than 2/3 of new generation
     The remainder is generated using former top players
     */
    
    inline NeuralNetwork<NBLAYERS>* trainNetworks(const int iterations,
                                                  const std::function<void(NeuralNetwork<NBLAYERS>&, NeuralNetwork<NBLAYERS>*, NeuralNetwork<NBLAYERS>*)> merging_instructions = BaseGeneticAlgorithm<NBLAYERS, INPUTNUMBER, OUTPUTNUMBER, POP_SIZE>::DEFAULT_MERGING_INSTRUCTIONS,
                                                  const std::function<float(std::array<StaticMatrix <float, OUTPUTNUMBER, 1>, TRAINING_LENGTH>&, std::array< StaticMatrix <float, OUTPUTNUMBER, 1> , TRAINING_LENGTH>&)> evalfunc = GeneticAlgorithm<NBLAYERS, INPUTNUMBER, OUTPUTNUMBER, TRAINING_LENGTH, POP_SIZE>::MSE) {
        for (std::size_t iteration(1); iteration <= iterations; iteration++){
            std::array< std::array< StaticMatrix <float, OUTPUTNUMBER, 1> , TRAINING_LENGTH>, POP_SIZE> localoutputs;
            std::array< std::pair<float, NeuralNetwork<NBLAYERS>*>, POP_SIZE> MSE;
            for (std::size_t netid(0); netid < POP_SIZE; netid++){
                float MSE_i = 0;
                for (std::size_t inputid(0); inputid < TRAINING_LENGTH; inputid++){
                    localoutputs[netid][inputid] = *(this->networksPopulation[netid]->template process<INPUTNUMBER, OUTPUTNUMBER>(inputs[inputid]));
                    if (this->postactivation_process != nullptr)
                        localoutputs[netid][inputid] = this->postactivation_process->operator()(localoutputs[netid][inputid]);
                }
                MSE_i = evalfunc(outputs, localoutputs[netid]);
                MSE[netid] = {MSE_i, this->networksPopulation[netid]};
            }
            std::sort(MSE.begin(), MSE.end(), [](const std::pair<float, NeuralNetwork<NBLAYERS>*> &x, const std::pair<float, NeuralNetwork<NBLAYERS>*> &y){
                return x.first > y.first;
            });
            
            this->generateNewGeneration(MSE, iteration, merging_instructions);
            
            if (MSE[POP_SIZE-1].first == 0.f){
                std::cout << "Best NN trained in " << iteration << " iterations (max allowed: " << iterations << ")\n";
                break;
            }
        }
        
        return this->networksPopulation[POP_SIZE-1];
        
    }
    
};

}

#endif /* GeneticAlgorithm_hpp */
