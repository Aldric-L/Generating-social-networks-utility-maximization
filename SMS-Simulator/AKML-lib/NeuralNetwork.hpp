//
//  NeuralNetwork.hpp
//  AKML Project
//
//  Created by Aldric Labarthe on 07/09/2023.
//

#ifndef NeuralNetwork_hpp
#define NeuralNetwork_hpp

#include <fstream>

#include "AKML_consts.hpp"
#include "UtilityLoops.hpp"
#include "Matrices.hpp"
#include "NeuralLayer.hpp"

namespace akml {

template <const std::size_t NBLAYERS=AKML_NEURAL_LAYER_NB>
class NeuralNetwork {
protected:
    const std::string customOriginField;
    std::array<AbstractNeuralLayer*, NBLAYERS> layers;
    
    template<std::size_t i, std::size_t j>
    struct setLayer_functor {
        static void run(akml::NeuralNetwork<NBLAYERS>& neuralnet, std::function<float(float)>& activation_func) {
            neuralnet.setLayer<akml::NN_STRUCTURE[j], akml::NN_STRUCTURE[i]>(i+2, activation_func);
        }
    };
    
    template<std::size_t i, std::size_t j>
    struct exportLayer_functor {
        static void run(std::array<AbstractNeuralLayer*, NBLAYERS>* layers, std::string* file_content) {
            NeuralLayer<akml::NN_STRUCTURE[j], akml::NN_STRUCTURE[i]>* layer = (NeuralLayer<akml::NN_STRUCTURE[j], akml::NN_STRUCTURE[i]>*) (*layers)[i+1];
            layer->saveLayer(file_content, i+1);
        }
    };
    
    template<std::size_t i, std::size_t j>
    struct importLayer_functor {
        static void run(std::array<AbstractNeuralLayer*, NBLAYERS>* layers, std::vector<std::vector<std::vector<std::string>>>* weights,
                        std::vector<std::vector<std::vector<std::string>>>* biases) {
            NeuralLayer<akml::NN_STRUCTURE[j], akml::NN_STRUCTURE[i]>* layer = (NeuralLayer<akml::NN_STRUCTURE[j], akml::NN_STRUCTURE[i]>*) (*layers)[i+1];
            // Biases :
            for (std::size_t row(0); row < layer->getNeuronNumber(); row++){
                layer->getBiasesAccess()->operator()(row+1, 1) = stof((*biases)[i+1][row][0]);
            }
            
            // weights :
            for (std::size_t row(0); row < layer->getNeuronNumber(); row++){
                for (std::size_t col(0); col < layer->getPreviousNeuronNumber(); col++){
                    layer->getWeightsAccess()->operator()(row+1, col+1) = stof((*weights)[i+1][row][col]);
                }
            }
            
        }
    };
    
public:
    
    static inline std::function<float(float)> SIGMOID = [](float x) {return 1/(1+exp(-x));};
    static inline std::function<float(float)> RELU = [](float x) {return std::max(0.f, x); };

    static inline std::function<float(float)> NO_ACTION = [](float x) {return x;};
    
    inline NeuralNetwork(std::string customOriginField="") : customOriginField(customOriginField) {
        for (std::size_t i(0); i < NBLAYERS; i++){
            layers[i] = nullptr;
        }
    };
    
    inline ~NeuralNetwork() {
        for (std::size_t i(0); i < NBLAYERS; i++){
            if (layers[i] != nullptr)
                delete layers[i];
        }
    };
    
    inline std::string getCustomOriginField(){ return customOriginField; }
    
    // This function only works with default initializer
    inline void saveNetwork(std::string filename="NNSave.akml"){
        std::ofstream file;
        if (filename != "cout"){
            file.open (filename);
        }
        
        std::string* file_content = new std::string;
        *file_content += "AKML_NN " + std::to_string(NBLAYERS) + "\n";
        
        *file_content += "AKML_NN_STRUCTURE ";
        for (std::size_t i(0); i < NBLAYERS; i++){
            if (layers[i] != nullptr)
                *file_content += std::to_string(layers[i]->neuronNumber) + " ";
        }
        *file_content += "\n";
        
        NeuralLayer<akml::NN_STRUCTURE[0], 1>* layer = (NeuralLayer<akml::NN_STRUCTURE[0], 1>*) (layers)[0];
        layer->saveLayer(file_content, 0);
    
        akml::for_<0, NBLAYERS-1>::template run<exportLayer_functor, std::array<AbstractNeuralLayer*, NBLAYERS>*, std::string*>(&layers, file_content);
        file << *file_content;
        if (filename != "cout"){
            file.close();
        }else{
            std::cout << "Exporting NeuralNetwork (mode COUT) \n";
            std::cout << *file_content;
        }
        
    }
    
    inline void importNetwork(std::string filename="NNSave.akml"){
        std::ifstream file;
        file.open (filename);
        std::string line;
        if (file.is_open()) {
            std::vector<std::vector<std::string>> buffer;
            std::vector<std::vector<std::vector<std::string>>> weights;
            std::vector<std::vector<std::vector<std::string>>> biases;

            while (file) {
                std::getline (file, line);
                if (line.length() > 0){
                    std::string word_buffer;
                    std::vector<std::string> line_buffer;
                    for (int incr(0); incr < line.length(); incr++){
                        if (line[incr] == ' ' && word_buffer.length() != 0){
                            line_buffer.push_back(word_buffer);
                            word_buffer = "";
                        }else{
                            word_buffer.push_back(line[incr]);
                        }
                    }
                    if (word_buffer.length() != 0){
                        line_buffer.push_back(word_buffer);
                        word_buffer = "";
                    }
                    if (line_buffer[0] == "AKML_NN"){
                        if (stoi(line_buffer[1]) != NBLAYERS)
                            throw std::invalid_argument("Unable to import this network : the size is invalid.");
                        else
                            line_buffer.clear();
                    }else if (line_buffer[0] == "AKML_NN_STRUCTURE"){
                        for (int lay(1); lay < line_buffer.size(); lay++){
                            if (stoi(line_buffer[lay]) != akml::NN_STRUCTURE[lay-1])
                                throw std::invalid_argument("Unable to import this network : the layers are invalid.");
                        }
                        line_buffer.clear();
                    }else if (line_buffer[0] == "AKML_LAYER_END_BIASES"){
                        biases.push_back(buffer);
                        buffer.clear();
                        line_buffer.clear();
                    }else if (line_buffer[0] == "AKML_LAYER_END_WEIGHTS"){
                        weights.push_back(buffer);
                        buffer.clear();
                        line_buffer.clear();
                    }else if (line_buffer[0] == "AKML_LAYER_WEIGHTS" || line_buffer[0] == "AKML_LAYER_BIASES"){
                        line_buffer.clear();
                    }else {
                        buffer.push_back(line_buffer);
                        line_buffer.clear();
                    }

                }
            }
            if (weights.size() != biases.size() || biases.size() == 0 || biases.size() != NBLAYERS)
                throw std::invalid_argument("Unable to import this network : the size is invalid (more biases than weights).");
            
            NeuralLayer<akml::NN_STRUCTURE[0], 1>* layer = (NeuralLayer<akml::NN_STRUCTURE[0], 1>*) (layers)[0];
            // Biases :
            for (std::size_t row(0); row < layer->getNeuronNumber(); row++){
                layer->getBiasesAccess()->operator()(row+1, 1) = stof((biases)[0][row][0]);
            }
            
            // weights :
            for (std::size_t row(0); row < layer->getNeuronNumber(); row++){
                for (std::size_t col(0); col < layer->getPreviousNeuronNumber(); col++){
                    layer->getWeightsAccess()->operator()(row+1, col+1) = stof((weights)[0][row][col]);
                }
            }
        
            akml::for_<0, NBLAYERS-1>::template run<importLayer_functor, std::array<AbstractNeuralLayer*, NBLAYERS>*, std::vector<std::vector<std::vector<std::string>>>*, std::vector<std::vector<std::vector<std::string>>>*>(&layers, &weights, &biases);
            
            //assert(1==2);
        }
        file.close();
    }
    
    template <std::size_t INPUTNUMBER>
    inline NeuralLayer<INPUTNUMBER, 1>* setFirstLayer(){
        layers[0] = new NeuralLayer<INPUTNUMBER, 1>(0);
        layers[0]->setFirstRow(true);
        layers[0]->setActivationFunction(NeuralNetwork<NBLAYERS>::NO_ACTION);
        return (NeuralLayer<INPUTNUMBER, 1>*)layers[0];
    };
    
    template <std::size_t NEURON_NUMBER, std::size_t PREVIOUS_NEURON_NUMBER>
    inline NeuralLayer<NEURON_NUMBER, PREVIOUS_NEURON_NUMBER>* setLayer(std::size_t layer=0, const std::function<float(float)> activation_function=NeuralNetwork<NBLAYERS>::SIGMOID){
        if (layer == 0){
            for (std::size_t i(0); i < NBLAYERS; i++){
                if (layers[i] == nullptr)
                    layer = i + 1;
            }
        }
        
        layers[layer-1] = new NeuralLayer<NEURON_NUMBER, PREVIOUS_NEURON_NUMBER>(layer-1);
        layers[layer-1]->setActivationFunction(activation_function);
        
        StaticMatrix<float, NEURON_NUMBER, 1> new_biases;
        new_biases.transform(StaticMatrix<float, NEURON_NUMBER, 1>::RANDOM_TRANSFORM);
        layers[layer-1]->setBiases(&new_biases);
        StaticMatrix<float, NEURON_NUMBER, PREVIOUS_NEURON_NUMBER> new_weights;
        new_weights.transform(StaticMatrix<float, NEURON_NUMBER, 1>::RANDOM_TRANSFORM);
        layers[layer-1]->setWeights(&new_weights);
        return (NeuralLayer<NEURON_NUMBER, PREVIOUS_NEURON_NUMBER>*)layers[layer-1];
    };
    
    inline AbstractNeuralLayer* getLayer(std::size_t layer){
        return layers[layer-1];
    };
    
    template <std::size_t INPUTNUMBER=NEURAL_INPUT_NB, std::size_t OUTPUTNUMBER=NEURAL_OUTPUT_NB>
    inline StaticMatrix<float, OUTPUTNUMBER, 1>* process(StaticMatrix<float, INPUTNUMBER, 1> &input){
        layers[0]->setInput(&input);
        for (std::size_t i(1); i < NBLAYERS; i++){
            if (layers[i] == nullptr)
                throw std::exception();
            layers[i]->setPreviousActivationLayer(layers[i-1]->getActivationLayer());
        }
        
        return  (StaticMatrix<float, OUTPUTNUMBER, 1>*)layers[NBLAYERS-1]->getActivationLayer();
    }
    
    static inline std::function<void(akml::NeuralNetwork<NBLAYERS>&)> DEFAULT_INIT_INSTRUCTIONS = [](NeuralNetwork<NBLAYERS>& net) {
        if (NBLAYERS != akml::NEURAL_LAYER_NB)
            throw std::invalid_argument("You should not use a default initializer for a neural network which is not standard in layers number.");
        
        net.setFirstLayer<akml::NN_STRUCTURE[0]>();
        akml::for_<0, NBLAYERS-1>::template run<setLayer_functor,akml::NeuralNetwork<NBLAYERS>&, std::function<float(float)>&>(net, akml::NeuralNetwork<NBLAYERS>::SIGMOID);
    };
    
    static inline std::function<void(akml::NeuralNetwork<NBLAYERS>&)> DEFAULT_INITRELU_INSTRUCTIONS = [](NeuralNetwork<NBLAYERS>& net) {
        if (NBLAYERS != akml::NEURAL_LAYER_NB)
            throw std::invalid_argument("You should not use a default initializer for a neural network which is not standard in layers number.");
        
        net.setFirstLayer<akml::NN_STRUCTURE[0]>();
        if (NBLAYERS < 3){
            akml::NeuralNetwork<NBLAYERS>::DEFAULT_INIT_INSTRUCTIONS;
            return;
        }
            

        akml::for_<0, NBLAYERS-2>::template run<setLayer_functor,akml::NeuralNetwork<NBLAYERS>&, std::function<float(float)>&>(net, akml::NeuralNetwork<NBLAYERS>::RELU);
        akml::for_<NBLAYERS-2, NBLAYERS-1>::template run<setLayer_functor,akml::NeuralNetwork<NBLAYERS>&, std::function<float(float)>&>(net, akml::NeuralNetwork<NBLAYERS>::SIGMOID);
    };
};

}

#endif /* NeuralNetwork_hpp */
