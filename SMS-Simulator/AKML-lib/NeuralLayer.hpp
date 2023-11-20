//
//  NeuralLayer.hpp
//  AKML Project
//
//  Created by Aldric Labarthe on 09/09/2023.
//

#ifndef NeuralLayer_hpp
#define NeuralLayer_hpp

#include "Matrix.hpp"

namespace akml {

class AbstractNeuralLayer {
protected:
    std::function<float(float)> activationFunction;
    bool isFirstRow = false;
    
public:
    const std::size_t neuronNumber, previousNeuronNumber, layerId;
    
    AbstractNeuralLayer(const std::size_t neuronNumber, const std::size_t previousNeuronNumber, const std::size_t layerId) : neuronNumber(neuronNumber), previousNeuronNumber(previousNeuronNumber), layerId(layerId) {}
    
    inline std::size_t getNeuronNumber(){ return neuronNumber; }
    inline std::size_t getPreviousNeuronNumber(){ return previousNeuronNumber; }
    
    inline void setActivationFunction(std::function<float(float)> actfunc){ activationFunction = actfunc; }
    
    inline void setFirstRow(const bool first_row){ isFirstRow = first_row; }
    
    inline virtual void setInput(MatrixPrototype<float>* arg) { return; };
    inline virtual void setBiases(MatrixPrototype<float>* new_biases) { return; };
    inline virtual void setWeights(MatrixPrototype<float>* new_weights) { return; };
    inline virtual void setPreviousActivationLayer(MatrixPrototype<float>* argument) { return; };
    inline virtual MatrixPrototype<float>* getActivationLayer() { MatrixPrototype<float>* rtrn = new MatrixPrototype<float>(1, 1); return rtrn; };
};


template <std::size_t NEURON_NUMBER, std::size_t PREVIOUS_NEURON_NUMBER>
class NeuralLayer : public AbstractNeuralLayer {
private:
    Matrix<float, NEURON_NUMBER, 1> ownActivationLayer;
    Matrix<float, PREVIOUS_NEURON_NUMBER, 1>* previousActivationLayer;
    Matrix<float, NEURON_NUMBER, PREVIOUS_NEURON_NUMBER> weights;
    Matrix<float, NEURON_NUMBER, 1> biases;
    
public:
    
    NeuralLayer(const std::size_t layerId) : AbstractNeuralLayer(NEURON_NUMBER, PREVIOUS_NEURON_NUMBER, layerId) {
        weights.transform([](float x, std::size_t row, std::size_t column) {return 1;});
        previousActivationLayer = nullptr;
    }
    
    inline ~NeuralLayer() {
        if (previousActivationLayer != nullptr)
            delete previousActivationLayer;
    };
    
    inline void setInput(MatrixPrototype<float>* arg) {
        if (!isFirstRow)
            throw std::exception();
        
        ownActivationLayer = *((Matrix<float, NEURON_NUMBER, 1>*)arg);
        if (previousActivationLayer == nullptr)
            delete previousActivationLayer;
        // This is oddly a memory leak... but thanksfully it is not necessary to define previous activation Layer for the first layer
        // It was only more "humanly-coherent"
        //previousActivationLayer = new Matrix<float, PREVIOUS_NEURON_NUMBER, 1>;
        //previousActivationLayer->operator()(1,1) = 1;
    }
    
    inline void setBiases(MatrixPrototype<float>* arg) {
        biases = *((Matrix<float, NEURON_NUMBER, 1>*)arg);
    }
    
    inline Matrix<float, NEURON_NUMBER, 1>* getBiasesAccess(){
        Matrix<float, NEURON_NUMBER, 1>* b_point (&biases);
        return b_point;
    }
    
    inline void setWeights(MatrixPrototype<float>* arg) {
        Matrix<float, NEURON_NUMBER, PREVIOUS_NEURON_NUMBER>* new_weights(0);
        new_weights = (Matrix<float, NEURON_NUMBER, PREVIOUS_NEURON_NUMBER>*)arg;
        weights = *new_weights;
    }
    
    inline Matrix<float, NEURON_NUMBER, PREVIOUS_NEURON_NUMBER>* getWeightsAccess (){
        Matrix<float, NEURON_NUMBER, PREVIOUS_NEURON_NUMBER>* w_point (&weights);
        return w_point;
    }
    
    inline Matrix<float, PREVIOUS_NEURON_NUMBER, 1>* getPreviousActivationLayer (){
        return previousActivationLayer;
    }
    
    inline void setPreviousActivationLayer(MatrixPrototype<float>* arg) {
        Matrix<float, PREVIOUS_NEURON_NUMBER, 1>* prev(0);
        prev = (Matrix<float, PREVIOUS_NEURON_NUMBER, 1>*)arg;
        previousActivationLayer = prev;
    }
    
    inline MatrixPrototype<float>* getActivationLayer(){
        if (!isFirstRow){
            ownActivationLayer = Matrix<float, NEURON_NUMBER, 1>::product(weights, *previousActivationLayer);
            ownActivationLayer += biases;
            ownActivationLayer.transform(activationFunction);
        }
        return &ownActivationLayer;
    }
    
    inline void saveLayer(std::string* buffer, int id){
        *buffer += ("AKML_LAYER_BIASES " + std::to_string(id) + "\n");
        // Biases :
        for (std::size_t row(0); row < this->getNeuronNumber(); row++){
            *buffer += std::to_string(this->getBiasesAccess()->operator()(row+1, 1)) + "\n";
        }
        
        *buffer += "AKML_LAYER_END_BIASES \n";
        *buffer += "AKML_LAYER_WEIGHTS " + std::to_string(id) + "\n";
        // weights :
        for (std::size_t row(0); row < this->getNeuronNumber(); row++){
            for (std::size_t col(0); col < this->getPreviousNeuronNumber(); col++){
                *buffer += std::to_string(this->getWeightsAccess()->operator()(row+1, col+1)) + " ";
            }
            *buffer += "\n";
        }
        *buffer += "AKML_LAYER_END_WEIGHTS \n";
    }
    
};

}

#endif /* NeuralLayer_hpp */
