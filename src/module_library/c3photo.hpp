#ifndef C3PHOTO_H
#define C3PHOTO_H

#include <map>
#include <string>
struct c3_str {
    double Assim;
    double Gs;
    double Ci;
    double GrossAssim;
    double iterTimes;
    double penalty;
};

struct c3_str c3photoC(
    double const Qp,
    double const Tleaf,
    double const RH,
    double const Vcmax0,
    double const Jmax0,
    double const TPU_rate_max,
    double const Rd0,
    double const bb0,
    double const bb1,
    double const Gs_min,
    double Ca,
    double const AP,
    double const O2,
    double const thet,
    double const StomWS,
    int const water_stress_approach,
    double const electrons_per_carboxylation,
    double const electrons_per_oxygenation,
    double const exp_id);

struct c3_str c3photoC_FvCB(
    double const Qp,
    double const Tleaf,
    double const RH,
    double const Vcmax0,
    double const Jmax0,
    double const TPU_rate_max,
    double const Rd0,
    double const bb0,
    double const bb1,
    double const Gs_min,
    double Ca,
    double const AP,
    double const O2,
    double const thet,
    double const StomWS,
    int const water_stress_approach,
    double const electrons_per_carboxylation,
    double const electrons_per_oxygenation);

double solc(double LeafT);
double solo(double LeafT);

struct ephoto_Result {
    double A;
    double penalty;
};
ephoto_Result assim_ephoto(double LeafT,double PAR, double Ci, double exp_id);

void readFile1(const std::string &filename, std::map<std::string,std::string> &mapper);
void readFile2(const std::string &filename, std::map<std::string, double> &mapper); 

#endif
