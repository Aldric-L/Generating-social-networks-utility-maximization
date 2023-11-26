//
//  Save.hpp
//  AKML Project
//
//  Created by Aldric Labarthe on 05/11/2023.
//

#ifndef Save_hpp
#define Save_hpp

#include "../AKML.hpp"
#include <tuple>

namespace akml {

    template <unsigned short int parameters_nb, typename... types>
    class Save {
        private:
            template<int i>
            struct exportLayer_functor {
                static std::string to_string(const std::string elem){ return elem; }
                template <typename T>
                static std::string to_string(const T elem){ return std::to_string(elem); }
                
                static void run(std::string& rtrn_to_edit, std::tuple<types...>& parameters) {
                    rtrn_to_edit += to_string(std::get<i>(parameters));
                    if (i != parameters_nb-1)
                        rtrn_to_edit += ", ";
                }
            };
        
        protected:
            std::array<std::string, parameters_nb> parameters_name;
            std::tuple<types...> parameters;
        
        public:
            static inline std::array<std::string, parameters_nb> default_parameters_name;
        
            inline Save(std::tuple<types...>& param) : parameters(param), parameters_name(default_parameters_name) {};
        
            inline Save(types... params) : parameters_name(default_parameters_name){
                parameters = std::make_tuple(params...);
            };
        
            inline Save() : parameters_name(default_parameters_name){};
        
            inline void setParameterNames(std::array<std::string, parameters_nb> pn){
                parameters_name = pn;
            }
        
            template <int param_id, typename type>
            inline void editParameter(type& param_val){
                std::get<param_id>(parameters) = param_val;
            }
        
            inline void editParameters(types... params){
                parameters = std::make_tuple(params...);
            }
        
            template <int param_id, typename type>
            inline type& getParameter(){
                return std::get<param_id>(parameters);
            }
            
            inline std::string printAsCSV() {
                std::string rtrn ("");
                akml::for_<0, parameters_nb>::template run<exportLayer_functor,std::string&, std::tuple<types...>&>(rtrn, parameters);
                return rtrn;
            }
        
            inline std::string printTitleAsCSV() {
                std::string rtrn ("");
                for (int i(0); i < parameters_nb; i++){
                    rtrn += parameters_name[i];
                    if (i != parameters_nb-1)
                        rtrn += ", ";
                }
                return rtrn;
            }
    };

}

#endif /* Save_hpp */
