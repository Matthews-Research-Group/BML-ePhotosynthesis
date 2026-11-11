#ifndef EPHOTO_ASSIM_H 
#define EPHOTO_ASSIM_H 
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>

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
// One biochemical seed per leaf context, separate from the exact-output cache.
struct ephoto_state {
    std::vector<double> values;
    double temperature = 0.0, ppfd = 0.0, ci = 0.0, experiment = 0.0;

    bool nearby(double t, double light, double carbon, double profile) const {
        return !values.empty() && std::isfinite(t) && std::isfinite(light) &&
            std::isfinite(carbon) && profile == experiment &&
            std::abs(t - temperature) <= 2.0 &&
            std::abs(light - ppfd) <= std::max(10.0, 0.1 * std::abs(ppfd)) &&
            std::abs(carbon - ci) <= std::max(10.0, 0.1 * std::abs(ci));
    }
};
ephoto_outputs assim_ephoto(double LeafT, double PAR, double Ci, double exp_id,
                            ephoto_state* state = nullptr);

void readFile1(const std::string &filename, std::map<std::string,std::string> &mapper);
void readFile2(const std::string &filename, std::map<std::string, double> &mapper); 

#endif
