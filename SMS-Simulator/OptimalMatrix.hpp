//
//  OptimalMatrix.hpp
//  SocialMatrixSimulation
//
//  Created by Aldric Labarthe on 10/06/2024.
//

#ifndef OptimalMatrix_hpp
#define OptimalMatrix_hpp

#include "Constants.hpp"
#include "AKML-lib/Matrices.hpp"
#include "AKML-lib/NeuralNetwork.hpp"
#include "Individual.hpp"

class OptimalMatrix {
    //just a convenient alias
    const std::function<float(const float)>& sig = akml::ActivationFunctions::SIGMOID.function;
    const std::function<float(const float)>& sigPrime = akml::ActivationFunctions::SIGMOID.derivative;
    
    std::vector<akml::DynamicMatrix<float>> PSProdBuffer;
    
    akml::DynamicMatrix<float> buildGlobalGradient(akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> adjacencyMatrix, const akml::Matrix<Individual*, GRAPH_SIZE, 1>& individuals, double tolerance= 1e-6);
    void followGradient(const akml::DynamicMatrix<float>& grad, akml::DynamicMatrix<float>& adjacencyMatrix, const akml::NeuralNetwork::GRADIENT_METHODS method=akml::NeuralNetwork::GRADIENT_METHODS::GRADIENT_ASCENT);
    
public:
    static inline unsigned short int MAX_PRECISION = 5;
    static inline unsigned long int MAX_ITERATIONS = 1000000;
    
    void setCompatibilityMatrix(const akml::DynamicMatrix<float>& compatibilityMatrix);
    
    std::pair<std::size_t, akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE>> compute(akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> adjacencyMatrix, const akml::Matrix<Individual*, GRAPH_SIZE, 1>& individuals, const std::size_t max_epochs = MAX_ITERATIONS, const double lr_moment1 = 0.9,const double lr_moment2 = 0.999, const double step_size = 0.001, const double epsilon=1e-8);
    
    akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> exportAffinityBuffer() const;
    
    akml::DynamicMatrix<float> computeObjectiveFunction(const akml::DynamicMatrix<float>& adjacencyMatrix, const akml::Matrix<Individual*, GRAPH_SIZE, 1>& individuals) const;
    

};


#endif /* OptimalMatrix_hpp */
