#include <cmath>
#include <fstream>
#include <iostream>
#include <complex>
#include <RandomNumber.hpp>
#include <bitset>
#include <algorithm>
#include <numeric>
#include <vector>

using namespace std;
using namespace Snu::Cnrc;
const int NN = 8;

typedef complex<double> dcomplex;
using Size = unsigned int;

vector<dcomplex> normalize(vector<dcomplex>& state) {
    double normalizationFactor = 0;
    for(Size i=0; i<state.size(); ++i) normalizationFactor += norm(state[i]);
    for(Size i=0; i<state.size(); ++i) state[i] /= sqrt(normalizationFactor);
    return state;
}

vector<dcomplex> vectorCal(vector<dcomplex>& state, vector<dcomplex>& k, const double dt, const double num) {
    vector<dcomplex> cc(state.size());
    for(Size i=0; i<state.size(); ++i) cc[i] = state[i]+k[i]*dt/num;
    return cc;
}

// Decay Lindblad Operator
double calDecayDeltaP(Size site, const double dt, vector<dcomplex>& state, const Size N, const double gamma) {
    double deltaP = 0;
    for(Size i=0; i<pow(2, N); ++i) {
        if(bitset<NN>(i)[site] == 1) {
            deltaP += gamma*dt*(norm(state[i]));
        }
    }
    return deltaP;
}

vector<dcomplex> emissionDecay(Size site, vector<dcomplex>& state, const Size N) {
    vector<dcomplex> cc(pow(2, N));
    for(Size i=0; i<pow(2, N); ++i) {
        if(bitset<NN>(i)[site] == 1) {
            cc[bitset<NN>(i).flip(site).to_ulong()] += state[i];
        }
    }
    cc = normalize(cc);
    return cc;
}

// Runge-Kutta method for the differential value of wavefunction
vector<dcomplex> dcdt(vector<dcomplex> state, const Size N, const double omega, const double gamma) {
    dcomplex I(0, 1);
    vector<dcomplex> cc(pow(2,N));
    for(Size j=0; j<N; ++j) {
        for(Size i=0; i<pow(2, N); ++i) {
            Size jPrev = j == 0 ? N-1 : j-1;
            Size jNext = j == N-1 ? 0 : j+1;
            Size jFlip = bitset<NN>(i).flip(j).to_ulong();
            if(bitset<NN>(i)[j] == 1) {
                cc[i] -= (dcomplex)(gamma/2.)*state[i]; // decay
            }
            cc[jFlip] -= I*(dcomplex)(omega*(bitset<NN>(i)[jPrev]+bitset<NN>(i)[jNext]))*(state[i]); // Quantum part
        }
    }
    return cc;
}

int main(int argc, char *argv[]) {
    const Size number_of_time_step = stoul(argv[1]);
    const Size N = stoul(argv[2]);
    const double omega = stod(argv[3]);
    const double number_of_ensemble = stoul(argv[4]);
    const double gamma = 1;
    vector<dcomplex> state(pow(2, N));
    vector<dcomplex> k1(pow(2, N));
    vector<dcomplex> k2(pow(2, N));
    vector<dcomplex> k3(pow(2, N));
    vector<dcomplex> k4(pow(2, N));
    double dt = 0.01;

    
    RandomRealGenerator rnd(0.0, 1.0);
    vector<double> orderParameter(number_of_time_step);
    vector<double> deltaP(N);
    Size siteEvolved;
    for(Size ensembleNum=0; ensembleNum<number_of_ensemble; ++ensembleNum) {
        for(Size i=0; i<pow(2,N)-1; ++i) {
            state[i] = {0, 0};
        }
        state[pow(2,N)-1] = {1, 0};

        for(Size t=0; t<number_of_time_step; ++t) {
            for(Size site=0; site<N; ++site) deltaP[site] = calDecayDeltaP(site, dt, state, N, gamma);
            partial_sum(deltaP.begin(), deltaP.end(), deltaP.begin());
            double p = rnd();
            siteEvolved = distance(deltaP.begin(), lower_bound(deltaP.begin(), deltaP.end(), p));
            if(siteEvolved != N) {
                state = emissionDecay(siteEvolved, state, N);
            }
            else { // no emission takes place
                k1 = dcdt(state, N, omega, gamma);
                k2 = dcdt(vectorCal(state, k1, dt, 2), N, omega, gamma);
                k3 = dcdt(vectorCal(state, k2, dt, 2), N, omega, gamma);
                k4 = dcdt(vectorCal(state, k3, dt, 1), N, omega, gamma);
                for(Size i=0; i<pow(2, N); ++i) {
                    //state[i] = (state[i] + dt/6.*(k1[i]+2.*k2[i]+2.*k3[i]+k4[i])) / (dcomplex)(sqrt(1.-deltaP[N-1]));
                    state[i] = (state[i] + dt/6.*(k1[i]+2.*k2[i]+2.*k3[i]+k4[i]));
                }
                state = normalize(state);
            }  
            for(Size i=0; i<pow(2, N); ++i) {
                orderParameter[t] += (bitset<NN>(i).count()*norm(state[i]));
            }
        }
    }
    for(Size t=0; t<number_of_time_step; ++t) cout << t*dt*gamma << '\t' << orderParameter[t]/number_of_ensemble/N << endl;
    return 0;
}