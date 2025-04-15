#ifndef EPHOTO_ASSIM_H 
#define EPHOTO_ASSIM_H 
#include <string>
#include <map>

struct ephoto_outputs {
    double A;
    double Vc;
    double PR;
    double penalty;
};
ephoto_outputs assim_ephoto(double LeafT,double PAR, double Ci, double exp_id);

void readFile1(const std::string &filename, std::map<std::string,std::string> &mapper);
void readFile2(const std::string &filename, std::map<std::string, double> &mapper); 

#endif
