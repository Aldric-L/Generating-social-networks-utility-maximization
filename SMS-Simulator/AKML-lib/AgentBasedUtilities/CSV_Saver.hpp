//
//  CSV-Saver.hpp
//  AKML Project
//
//  Created by Aldric Labarthe on 05/11/2023.
//

#ifndef CSV_Saver_hpp
#define CSV_Saver_hpp

#include <stdio.h>
#include <vector>
#include <fstream>
#include <iomanip>

#include "Save.hpp"

namespace akml {

    template <class SaveClass>
    class CSV_Saver {
    protected:
        std::vector<std::pair<std::size_t, SaveClass*>> memory;
        
    public:
        inline void addSave(SaveClass* iteration) {
            memory.push_back(std::make_pair((memory.size() == 0) ? 1 : memory.back().first+1, iteration));
        }
        
        inline void saveToCSV(const std::string filepathandname="data.csv") {
            std::ofstream file;
            file.open (filepathandname);
            file << "iteration," << memory[0].second->printTitleAsCSV() <<"\n";
            for (int it(0); it < memory.size(); it++){
                file << it << "," << memory[it].second->printAsCSV() <<"\n";
            }
            file.close();
        }
        
        inline ~CSV_Saver() {
            for (int it(0); it < memory.size(); it++){
                delete memory[it].second;
            }
        }
        
    };

}

#endif /* CSV_Saver_hpp */
