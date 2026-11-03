#ifndef EPHOTO_ASSIM_H 
#define EPHOTO_ASSIM_H 
#include <string>
#include <map>

enum class ephoto_status { success, timed_out, retries_exhausted, invalid_output };

inline const char* ephoto_status_message(ephoto_status status) {
    switch (status) {
        case ephoto_status::success: return "success";
        case ephoto_status::timed_out: return "timeout";
        case ephoto_status::retries_exhausted: return "retries_exhausted";
        default: return "invalid_output";
    }
}

struct ephoto_outputs {
    double A;
    double Vc;
    double PR;
    double penalty;
    ephoto_status status = ephoto_status::invalid_output;
};
ephoto_outputs assim_ephoto(double LeafT,double PAR, double Ci, double exp_id);

void readFile1(const std::string &filename, std::map<std::string,std::string> &mapper);
void readFile2(const std::string &filename, std::map<std::string, double> &mapper); 

#endif
