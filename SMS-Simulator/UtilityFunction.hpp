//
//  UtilityFunction.hpp
//  SocialMatrixSimulation
//
//  Created by Aldric Labarthe on 10/06/2024.
//

#ifndef UtilityFunction_hpp
#define UtilityFunction_hpp

#include "Constants.hpp"

class UtilityFunction {
public:
    virtual float function(const akml::DynamicMatrix<float>& affinityVector, const akml::DynamicMatrix<float>& weightsVector) const = 0;
    virtual float localDerivative(const float& affinityScore, const float& weight, const akml::DynamicMatrix<float>& weightsVector) const = 0;
    
    virtual akml::DynamicMatrix<float> derivative(const akml::DynamicMatrix<float>& affinityVector, const akml::DynamicMatrix<float>& weightsVector) const = 0;
    virtual ~UtilityFunction() = default;
};

class KZALUtility : public UtilityFunction {
public:
    const float gamma;
    const float kappa;
    const float delta;
    
    inline KZALUtility(float gamma, float kappa, float delta) : gamma(gamma), kappa(kappa), delta(delta) {};
    
    float function(const akml::DynamicMatrix<float>& affinityVector, const akml::DynamicMatrix<float>& weightsVector) const;
    float localDerivative(const float& affinityScore, const float& weight, const akml::DynamicMatrix<float>& weightsVector) const;
    akml::DynamicMatrix<float> derivative(const akml::DynamicMatrix<float>& affinityVector, const akml::DynamicMatrix<float>& weightsVector) const;
};
#endif /* UtilityFunction_hpp */
