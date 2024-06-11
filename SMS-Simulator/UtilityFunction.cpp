//
//  UtilityFunction.cpp
//  SocialMatrixSimulation
//
//  Created by Aldric Labarthe on 10/06/2024.
//

#include "UtilityFunction.hpp"

float ALKYUtility::function(const akml::DynamicMatrix<float>& affinityVector, const akml::DynamicMatrix<float>& weightsVector) const {
    float LHS = akml::inner_product(weightsVector, affinityVector)*kappa;
    
    float RHS1 = 0;
    float RHS2 = 0;
    for (std::size_t i(0); i < weightsVector.getNRows(); i++){
        if (weightsVector[{i, 0}] > 0){
            RHS1 += std::pow(weightsVector[{i, 0}], gamma) / (1.f - std::pow(weightsVector[{i, 0}], gamma));
            RHS2 += weightsVector[{i, 0}];
        }
    }
    RHS2 = std::pow(RHS2, delta);
    return LHS - RHS1 - RHS2;
}

float ALKYUtility::localDerivative(const float& affinityScore, const float& weight, const akml::DynamicMatrix<float>& weightsVector) const {
    return affinityScore * kappa - ((gamma * std::pow(weight, gamma-1)) / std::pow(1 - std::pow(weight, gamma), 2)) - delta * std::pow(akml::sum_column(weightsVector), delta-1);
}

akml::DynamicMatrix<float> ALKYUtility::derivative(const akml::DynamicMatrix<float>& affinityVector, const akml::DynamicMatrix<float>& weightsVector) const {
    akml::DynamicMatrix<float> gradient(weightsVector.getNRows(), 1);
    
    if (affinityVector.getNRows() != gradient.getNRows())
        throw std::runtime_error("Matrix dimension failure while computing utility");

    float sumWeights = delta * std::pow(akml::sum_column(weightsVector), delta-1);
    
    for (std::size_t i(0); i < gradient.getNRows(); i++)
        gradient[{i, 0}] = (weightsVector[{i, 0}] > 0) ? affinityVector[{i, 0}] * kappa - sumWeights - (gamma * std::pow(weightsVector[{i, 0}], gamma-1) / std::pow(1.f - std::pow(weightsVector[{i, 0}], gamma), 2)) : 0;
    
    return gradient;
}
