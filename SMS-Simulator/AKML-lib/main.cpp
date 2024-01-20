
//                                                                 AA     LL
//                 pzr--}fC/][dn1wYnfJ                            AAAA    LL
//              c{?_x[v){Ur}][ft(vX(|wzrmOcQ                     AA |AA   LL
//            #][_j(L}u-[[{x[{)q|)fztfXxUQuJkUm                 AAAAAAA   LL
//        .])-1-U)t}-v})}}]L/(}hZjjmtXdYCZuCkJCnLC|Xx          |AA   AAA  LLLLLLLL
//      I{(}11_-+/-q{p}f[1Oj{){QnfL|J00kC*#aqh0zQnjvwz/fCj^
//    "tf1f(fxtujXzt(jcf{{vX/t|bOxJ8hmJ0zwLwJCqzCxYLrYXXYcxx?
//   >()vxcJpUaXMUfdtcuj{]XrtrkkQC0#WhwkhokwCCzChbdcQB@BuO|vxn[
//  "f(|xcOhW&OB@B*BMvJ//(cCzXmqzhB%Mo*qrbUZLbma0&paB@B&X*Zw)xruLx/"
//  rv/    a%Mw88(  ./jujn0nvu   mh%%B%B%ad*Qk0Y0#B@@@@@@*puQqdt/ntCujx/)
// .xu  O  oh88c  .`.1ptL|OXY   @@@@c   #mJhd_        BBqhoLYkd#ju/j|rtnfnrjxjx~'
// `YuO  doBa`    nOCOXmQhmw    @@@B"    #ozbJk,         {?{jb*QbO#W*8*#mha0LxfIl`
//  cCZwdk&,      jzzQcak*u.    %%%%B    WoddJmO_            ddddbM*p/U0dbM*p/Ueudk
//   hMW          J[fY:          wBB8    M*p/U0{                 wLwJCqzCxYwLwJCqzfdfd
//                !qv              !qv   `}czjf:                    wLwJCqzCxYwLwJCddddfd
//

// This file provides an example of the usage of AKML for simple problems like the famous XOR
#define XOR_DIM 3

// AKML needs to first have the default layout of NeuralNetworks defiened
#if XOR_DIM == 3
    #define AKML_NEURAL_LAYER_NB 4
    #define AKML_NN_STRUCTURE { 3, 6, 3, 1 }
#endif

#if XOR_DIM == 2
    #define AKML_NEURAL_LAYER_NB 3
    #define AKML_NN_STRUCTURE { 2, 2, 1 }
#endif

// AKML.hpp contains the basics utilities needed to build a NeuralNetwork
#include "AKML.hpp"
// If you want to use a GA, you need to include it (it is not in the AKML.hpp)
#include "GeneticAlgorithm.cpp"

int main(int argc, const char * argv[]) {
    using SINGLETON = std::array <std::array <float, 1>, 1>;

    std::cout << "Hello, terrible World!\n";
    
    // NeuralNetwork initialization example
    //akml::NeuralNetwork neuralnet;
    /*
    Manual initialization example :
    neuralnet.setFirstLayer<3>();
    neuralnet.setLayer<6, 3>(2);
    neuralnet.setLayer<3, 6>(3);
    neuralnet.setLayer<1, 3>(4);*/
    //akml::NeuralNetwork<>::DEFAULT_INIT_INSTRUCTIONS(neuralnet);
    
    #if XOR_DIM == 3
        std::array<akml::StaticMatrix <float, 3, 1>, 8> inputs =
        {{ {{ {1,0,0} }}, {{ {0,0,1} }}, {{ {0,1,0} }},
           {{ {1,1,0} }}, {{ {1,0,1} }}, {{ {0,1,1} }},
           {{ {0,0,0} }}, {{ {1,1,1} }} }};
        
        std::array<akml::StaticMatrix <float, 1, 1>, 8> outputs =
        {{ (SINGLETON){{ {1} }}, (SINGLETON){{ {1} }}, (SINGLETON){{ {1} }},
           (SINGLETON){{ {0} }}, (SINGLETON){{ {0} }}, (SINGLETON){{ {0} }},
           (SINGLETON){{ {0} }}, (SINGLETON){{ {0} }} }};
        

        akml::GeneticAlgorithm<4, 3, 1, 8> ga (inputs, outputs);
    #endif
    
    #if XOR_DIM == 2
        std::array<akml::Matrix <float, 2, 1>, 4> inputs = {{ {{ {1,0} }}, {{ {0,0} }}, {{ {0,1} }}, {{ {1,1} }} }};
        
        std::array<akml::Matrix <float, 1, 1>, 4> outputs = {{ (SINGLETON){{ {1} }}, (SINGLETON){{ {0} }}, (SINGLETON){{ {1} }}, (SINGLETON){{ {0} }} }};
        
        akml::GeneticAlgorithm<3, 2, 1, 4> ga (inputs, outputs);
    #endif
    akml::NeuralNetwork<>* bestnet = ga.trainNetworks(5000);
    
    for (std::size_t inputid(0); inputid<inputs.size(); inputid++){
        std::cout << "Testing with " << std::endl;
        std::cout << inputs[inputid];
        std::cout << "Output :" << std::endl;
        std::cout << *bestnet->process<>(inputs[inputid]);
        std::cout << "Output expected :" << std::endl;
        std::cout << outputs[inputid];
    }
    bestnet->saveNetwork("cout");
    
    
    return 0;
}
