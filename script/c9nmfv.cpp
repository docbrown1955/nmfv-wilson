// This file is part of NMFV_WILSONS.
//
// NMFV_WILSONS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// NMFV_WILSONS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with NMFV_WILSONS. If not, see <https://www.gnu.org/licenses/>.

#include "c9_nmfv.h"
#include "lha.h"
#include "read_spectrum.h"
#include "lhaData.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <fstream>
//#include <filesystem>
#include "ckm.h"
#include "marty/looptools_interface.h"

using namespace c9_nmfv;

enum Mode {
    Mass, C7, C7p, C9, C10, DiM, DiMp
};

bool isFunction(std::string_view name, Mode mode)
{
    switch (mode) {
        case Mass: return name[0] == 'm';
        case C7:   return name[1] == '7' && name[2] != 'p';
        case C7p:  return name[1] == '7' && name[2] == 'p';
        case C9:   return name[1] == '9';
        case C10:  return name[1] == '1';
        case DiM:  return name[0] == 'D' && name[6] != 'p';
        case DiMp: return name[0] == 'D' && name[6] == 'p';
    }
    return false;
}
int main(int argc, char ** argv) {
  // In case of parallel computation on lists
  int pid;
  if (argc != 2) {
    std::cout << "This program takes exactly 1 argument : \n" << "int id_file\t \t int id of file number to treat from split data\n" ;
      exit(EXIT_FAILURE);
    }
    else{
      std::istringstream iss(argv[1]);
      if (iss >> pid) {
        std::cout << "Arg line parser successful, proceeding now.\n";
      }
      else{
        std::cout << "Unexpected behavior, check input : "<< argv[1] << "\n";
        exit(EXIT_FAILURE);
      }
    }
  
    char buf[100];
    int char_count_1; int char_count_2;
    char filename[200];

    std::string folder = "/data/boussejra/slha2_list";
    
    if( pid >=0 && pid < 10){
      char_count_1 = snprintf(buf,100,"/data/boussejra/cw_NMFV_mty_0%d/",pid);
      char_count_2 = snprintf(filename, 200, "/data/boussejra/slha2_list/slha2_list_0%d",pid);
    }
    else if(pid < 100){
      char_count_1 = snprintf(buf,100,"/data/boussejra/cw_NMFV_mty_%d/",pid);
      char_count_2 = snprintf(filename, 200, "/data/boussejra/slha2_list/slha2_list_%d",pid);
    }
    else{
      std::cout << "Error while creating list at " << folder << " with pid : " << pid << std::endl;
      exit(EXIT_FAILURE);
    }
    if(char_count_1 < 0 or char_count_1 >100){
      std::cout << "Error while writing datadir name, char_count_1 = %d " <<  char_count_1 << "\n";
      exit(EXIT_FAILURE);
    }
    if (char_count_2 <0 or char_count_2>200) {
      std::cout << "List name creation failed. Exiting...\n" ;
      exit(EXIT_FAILURE);
    }

    const std::string analysisName = buf ;

    std::cout << "Treating : " << filename << " in folder : " << analysisName << '\n';
    
    std::ofstream other_re(analysisName + "other_contrib_re.txt");
    other_re.precision(10);
    other_re << std::left;

    std::ofstream global_re(analysisName + "total_contrib_re.txt");
    global_re.precision(10);
    global_re << std::left;

    std::ofstream particular_re(analysisName + "all_contrib_re.txt");
    particular_re.precision(10);
    particular_re << std::left;

    std::ofstream global_im(analysisName + "total_contrib_im.txt");
    global_im.precision(10);
    global_im << std::left;

    std::ofstream particular_im(analysisName + "all_contrib_im.txt");
    particular_im.precision(10);
    particular_im << std::left;

    std::ofstream parametrization(analysisName + "parametrization.txt");
    for (const auto &f : f_G) {
        if (!isFunction(f.name, Mass))
            parametrization << f.name << '\t';
    }
    parametrization.close();

    param_t params;
    fillCKM(params);
    std::vector<std::string> fileNames;
    fileNames.reserve(static_cast<int>(1e6));
    std::cout << "Press enter to continue" << '\n';
    std::cin.get();
    std::ifstream fileNameIn(filename);
    std::string name;
    while (!(fileNameIn >> name).eof()) {
        std::cout << "New file : " << name << '\n';
        fileNames.push_back(name);
    }
    std::cout << fileNames.size() << " total files" << '\n';
    std::cout << "Press enter to continue" << '\n';
    std::cin.get();
    const size_t nFiles = fileNames.size();
    std::ofstream out(analysisName + "masses.txt");
    int bad_files=0;
    for (size_t i = 0; i != nFiles; ++i) {
        std::cout << "Iteration " << i+1 << " / " << nFiles << '\n';
        std::string fileName =  fileNames[i];
        if (auto fileData = std::ifstream(fileName); fileData) {
            auto data = mty::lha::Reader::readFile(fileName);
            if (!data.empty()) {
                //readSpectrum(params, data, 1, false);
              try{
                  readSpectrum(params, data);
              }catch(std::bad_optional_access const& exception) {
                  std::cout << "Error during handling file";
                  bad_files++;
                  continue;
              }
                for (auto value : data.getValues("VALUES")) {
                    other_re << value << '\t';
                }
                other_re << std::endl;
            }
            else {
                std::cerr << "No file at " << i << " !" << '\n';
                continue;
            }
        } else {
            std::cerr << "No file for " << fileName << "at position " << i << " !" << '\n';
            continue;
        }
        out << params.m_C1p << " " << params.m_N_1 << " " << params.m_sb_L
            << " " << params.m_sc_L << '\n';
        std::cout << "masses file saved at : " << analysisName +"masses.txt" << "\n";
        // continue;
        setMu((double)crealq(params.M_W));
        params.reg_prop = 1e-5;

        complex_t c7 {0}, c7p {0}, c9 {0}, c10 {0}, dim {0}, dimp {0};
        for (const auto &f : f_G) {
            if (isFunction(f.name, Mass))
                continue;
            if (isFunction(f.name, DiM) || isFunction(f.name, DiMp)) {
                setMu((double)crealq(params.M_W));
                params.s_12 = params.m_mu*params.m_mu;
            }
            else {
                setMu((double)crealq(params.M_W));
                params.s_12 = (params.m_b*params.m_b + params.m_s*params.m_s) / 2;
                params.s_34 = -params.m_mu*params.m_mu;
            }
            if(params.m_mu == 0)
              std::cout << "Catching error : params.m_mu="<< params.m_mu << "\n";
            
            std::cout << params.m_C1p << " " << params.m_N_1 << " " << params.m_sb_L
            << " " << params.m_sc_L << '\n';
        // continue;
            clearcache();
            auto value = f(params);
            clearcache();

            particular_re << std::setw(20) << (double)crealq(value);
            particular_im << std::setw(20) << (double)cimagq(value);
            if (isFunction(f.name, C7))   c7   += value;
            if (isFunction(f.name, C7p))  c7p  += value;
            if (isFunction(f.name, C9))   c9   += value;
            if (isFunction(f.name, C10))  c10  += value;
            if (isFunction(f.name, DiM))  dim  += value;
            if (isFunction(f.name, DiMp)) dimp += value;
        }
        particular_re << std::endl;
        particular_im << std::endl;

        global_re << std::setw(20) << (double)(crealq(c7));
        global_re << std::setw(20) << (double)(crealq(c7p));
        global_re << std::setw(20) << (double)(crealq(c9));
        global_re << std::setw(20) << (double)(crealq(c10));
        global_re << std::setw(20) << (double)(crealq(dim));
        global_re << std::setw(20) << (double)(crealq(dimp)) << std::endl;

        global_im << std::setw(20) << (double)(cimagq(c7));
        global_im << std::setw(20) << (double)(cimagq(c7p));
        global_im << std::setw(20) << (double)(cimagq(c9));
        global_im << std::setw(20) << (double)(cimagq(c10));
        global_im << std::setw(20) << (double)(cimagq(dim));
        global_im << std::setw(20) << (double)(cimagq(dimp)) << std::endl;
    }

    global_re.close();
    global_im.close();
    particular_re.close();
    particular_im.close();
    std::cout << "Program temrianted successfully. Bad files :" << bad_files << "/" << nFiles << "." << '\n';
    return 0;
}

