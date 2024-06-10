//
//  OptimalMatrix.cpp
//  SocialMatrixSimulation
//
//  Created by Aldric Labarthe on 10/06/2024.
//

#include "OptimalMatrix.hpp"

akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> OptimalMatrix::compute(akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> adjacencyMatrix, const akml::Matrix<Individual*, GRAPH_SIZE, 1>& individuals) {
    
    // We first have to build the P_S Matrix
    std::vector<akml::Matrix<float, P_DIMENSION, 1>> P_S_temp(individuals.getNRows());
    for (std::size_t i(0); i < GRAPH_SIZE-1; i++)
        P_S_temp[i] = individuals[{i,0}]->getP();
    akml::DynamicMatrix<float> P_S (P_S_temp);
    
    // We bufferize all compatibility scores that are constants at each iteration
    for (std::size_t indiv(0); indiv < individuals.getNRows(); indiv++){
        akml::DynamicMatrix<float> P_prod (akml::matrix_product(akml::transpose(P_S), individuals[{indiv, 0}]->getP()));
        P_prod = P_prod * (1/(float)P_DIMENSION);
        PSProdBuffer.push_back(std::move(P_prod));
    }
    
    //std::cout << "Before we start we are at: \n";
    //std::cout << computeObjectiveFunction(adjacencyMatrix, individuals);
    
    // #############################################################################
    // ############## ADAM Optimizer and gradient ascent
    
    const std::size_t max_epochs = 100000;
    double lr_moment1 = 0.9;
    double lr_moment2 = 0.999;
    double step_size = 0.001;
    double tolerance = 1e-6;
    double epsilon=1e-8;
    std::size_t epochs(1);

    akml::DynamicMatrix<float> grad (GRAPH_SIZE*GRAPH_SIZE, 1);
    for (std::size_t i(0); i < grad.getNRows(); i++)
        grad[{i, 0}] = 0.f;
    akml::DynamicMatrix<float> firstMoment (grad);
    akml::DynamicMatrix<float> secondMoment(firstMoment);
    akml::DynamicMatrix<float> correctedFirstMoment(firstMoment);
    akml::DynamicMatrix<float> correctedSecondMoment(firstMoment);
                                                                            
    do {
        grad = buildGlobalGradient(adjacencyMatrix, individuals);
        firstMoment = lr_moment1 * firstMoment + (1.f - lr_moment1)*grad;
        secondMoment = lr_moment2 * secondMoment + (1.f - lr_moment2)*akml::hadamard_product(grad, grad);
            
        correctedFirstMoment = firstMoment * (1/(1-std::pow(lr_moment1,epochs)));
        correctedSecondMoment = secondMoment * (1/(1-std::pow(lr_moment2,epochs)));
            
        grad = step_size * akml::hadamard_division(correctedFirstMoment,
                                                akml::transform(correctedSecondMoment, (std::function<float(float)>) ([epsilon](float x) {return std::sqrt(x) + epsilon; })));
        
        followGradient(grad, adjacencyMatrix);
        epochs++;
    }while (epochs < max_epochs);
    adjacencyMatrix.transform(sig);
    //std::cout << "AlphaCorrected:\n" << correctedAdjacencyMatrix;
    //std::cout << "\nResult:\n" << computeObjectiveFunction(adjacencyMatrix, individuals);
    return adjacencyMatrix;
    
}

akml::DynamicMatrix<float> OptimalMatrix::buildGlobalGradient(akml::Matrix<float, GRAPH_SIZE, GRAPH_SIZE> adjacencyMatrix, const akml::Matrix<Individual*, GRAPH_SIZE, 1>& individuals) {
    akml::DynamicMatrix<float> finalGradient (GRAPH_SIZE*GRAPH_SIZE, 1);
    
    // We apply directly sigmoid on our local version of adjacency matrix to stick in the (0,1) interval
    adjacencyMatrix.transform(sig);
    
    for (std::size_t indiv(0); indiv < individuals.getNRows(); indiv++){
        akml::DynamicMatrix<float> localAlpha = akml::getRowAsColumn(adjacencyMatrix, indiv+1);
        akml::DynamicMatrix<float> localGrad = individuals[{indiv, 0}]->getUtilityFunction()->derivative(PSProdBuffer.at(indiv), localAlpha);
        
        localGrad = akml::hadamard_product(localGrad, akml::transform(localAlpha, sigPrime));
        
        std::copy(localGrad.getStorage(), localGrad.getStorageEnd(), finalGradient.getStorage() + GRAPH_SIZE*indiv);
    }
    return finalGradient;
}


void OptimalMatrix::followGradient(const akml::DynamicMatrix<float>& grad, akml::DynamicMatrix<float>& adjacencyMatrix, const akml::NeuralNetwork::GRADIENT_METHODS method) {
        for (std::size_t i(0); i < adjacencyMatrix.getNRows(); i++){
            for (std::size_t j(0); j < adjacencyMatrix.getNColumns(); j++)
                adjacencyMatrix[{i, j}] = adjacencyMatrix[{i, j}] + method * grad[{i*adjacencyMatrix.getNColumns()+j,0}];
        }
}

akml::DynamicMatrix<float> OptimalMatrix::computeObjectiveFunction(const akml::DynamicMatrix<float>& adjacencyMatrix, const akml::Matrix<Individual*, GRAPH_SIZE, 1>& individuals) {
    akml::DynamicMatrix<float> utilities (GRAPH_SIZE, 1);
    for (std::size_t indiv(0); indiv < individuals.getNRows(); indiv++)
        utilities[{indiv, 0}] = individuals[{indiv, 0}]->getUtilityFunction()->function(PSProdBuffer.at(indiv),  akml::getRowAsColumn(adjacencyMatrix, indiv+1));
    return utilities;
}
