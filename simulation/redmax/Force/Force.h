#pragma once
#include "Common.h"
#include "Utils.h"

namespace redmax {

class Simulation;

class Force {
public:
    // Simulation
    Simulation* _sim;

    // time report
    class Timer {
        std::map<std::string, long long> _time;
    public:
        void add(std::string name, long long duration) {
            auto it = _time.find(name);
            if (it == _time.end()) {
                _time[name] = duration;
            } else {
                it->second += duration;
            }
        }

        void print_time_report() {
            for (auto it = _time.begin();it != _time.end();it++) {
                std::cerr << "|" << std::setw(20) << std::left << it->first << "|" << std::setw(7) << std::right << it->second / 1000 << "(ms) |" << std::endl;
            }
        }

        void clear() {
            _time.clear();
        }
    } _time;

    Force(Simulation* sim) { _sim = sim; }

    virtual void init() {}
    virtual void computeForce(VectorX& fm, VectorX& fr, bool verbose = false) = 0;
    virtual void computeForceWithDerivative(VectorX& fm, VectorX& fr, MatrixX& Km, MatrixX& Dm, MatrixX& Kr, MatrixX& Dr, bool verbose = false) = 0;
    virtual void computeForceWithDerivative(
        VectorX& fm, VectorX& fr, 
        MatrixX& Km, MatrixX& Dm, 
        MatrixX& Kr, MatrixX& Dr, 
        MatrixX& dfm_dp, MatrixX& dfr_dp, bool verbose = false) {};

    virtual void test_derivatives_runtime() {}
    virtual void test_design_derivatives_runtime() {}
    
    void reset_time_report() {
        _time.clear();
    }
    void print_time_report() {
        _time.print_time_report();
    }
};

}