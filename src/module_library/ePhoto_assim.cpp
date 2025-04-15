#include <stdlib.h>
//#include <string>
//#include <map>
#include <sstream>
#include <vector>
#include <iterator>
#include <boost/algorithm/string_regex.hpp>
#include <boost/regex.hpp>
#include "cxxopts.hpp"
#include "globals.hpp"
#include "modules/trDynaPS.hpp"
#include "modules/CM.hpp"
#include "modules/EPS.hpp"
#include "modules/PR.hpp"
#include "drivers/drivers.hpp"
#include "Variables.hpp"
#include "ePhoto_assim.h"

const boost::regex token("\\s+");
using namespace ePhotosynthesis;
using namespace ePhotosynthesis::drivers;
using namespace ePhotosynthesis::modules;

ephoto_outputs assim_ephoto(double LeafT, double PAR, double Ci,double exp_id)
{
        bool record    = true;//turn this on to calculate penalties
        bool runBioCro = true;//turn this on to NOT substract the penalties from A
        double stoptime = 5000.0, begintime=0.0, stepsize=0.5;
        int  maxSubSteps = 2500;
//        double abstol = 9.9e-6, reltol = 1e-4;
//        double abstol = 9.9e-6, reltol = 1e-5;
        double abstol = 1e-6, reltol = 1e-6;
        std::string evn_input="InputEvn.txt";
        std::string atpcost="InputATPCost.txt";
        std::string enzymeFile="Einput_all_experiments/Einput7_"+std::to_string(static_cast<int>(exp_id))+".txt";
        std::map<std::string, std::string> inputs;

        readFile1(evn_input, inputs);
        readFile1(atpcost, inputs);

        Variables *theVars = new Variables();

        readFile2(enzymeFile, theVars->EnzymeAct);

        theVars->CO2_in = Ci;  //get value from BioCro 
        theVars->TestLi = PAR; //get value from BioCro 
        if (stoi(inputs.at("SucPath"), nullptr) > 0)
            CM::setTestSucPath(true);
        theVars->TestATPCost = stoi(inputs.at("ATPCost"), nullptr);
        theVars->record = record;
        theVars->runBioCro = runBioCro;
        theVars->useC3 = true;
        theVars->RUBISCOMETHOD = 2;
	PR::setRUBISCOTOTAL(3);
//this is the scaling factor for some enzymes. check EPS_Drive.cpp to see
//which enzymes are being scaled
//        theVars->sensitivity_sf = enzyme_sf;

        Driver *maindriver;
       
//        std::cout<< "leafT is "<<LeafT<<std::endl;
//        std::cout<< "PAR is "<<theVars->TestLi<<std::endl;
//        std::cout<< "co2 is "<<theVars->CO2_in<<std::endl;
//        std::cout<< "alpha2 is"<<theVars->alpha2<<std::endl;
//        std::cout<< "JMAX is"<<theVars->EnzymeAct.at("Jmax")<<std::endl;

        maindriver = new EPSDriver(theVars, begintime, stepsize, stoptime, maxSubSteps, abstol, reltol, 1, 1, LeafT);

        std::vector<double> ResultRate = maindriver->run();

       if (theVars != nullptr) {
           maindriver->inputVars= nullptr;
           delete theVars;
       }   
       delete maindriver;
       
       ephoto_outputs result; 
       result.A       = ResultRate[0];
       result.Vc      = ResultRate[1];
       result.PR      = ResultRate[2]; //photorepiration
       result.penalty = ResultRate[3];
       return result;
}

void readFile1(const std::string &filename, std::map<std::string, std::string> &mapper) {
    std::vector<std::string> tempVec;
    std::string input;
    std::ifstream inputfile(filename);
    if(inputfile.fail()) {
        std::cout << "Could not open " << filename << " for reading" << std::endl;
        exit(EXIT_FAILURE);
    }
    while (getline(inputfile, input)) {
        if (input.empty())
            return;
        boost::algorithm::split_regex(tempVec, input, token);
        mapper.insert(std::pair<std::string, std::string>(tempVec[0], tempVec[1]));
    }
}

void readFile2(const std::string &filename, std::map<std::string, double> &mapper) {
    std::vector<std::string> tempVec;
    std::string input;
    std::ifstream inputfile(filename);
    if(inputfile.fail()) {
        std::cout << "Could not open " << filename << " for reading" << std::endl;
        exit(EXIT_FAILURE);
    }
    int count = 0;
    while (getline(inputfile, input)) {
        if (input.empty())
            return;
        boost::algorithm::split_regex(tempVec, input, token);
        double d;
        std::stringstream ss(tempVec[1]);
        ss >> d;
        if (count < 27)
            d /= 30.;
        count++;
        mapper.insert(std::pair<std::string, double>(tempVec[0], d));
    }
}
