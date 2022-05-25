#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <string>
#include <thread>
#include <future>
#include <sstream>
#include <iomanip>
#include <math.h>

#include "channels.hpp"
#include "convolutional_code.hpp"

#define MAX_SIMS 10000
#define MIN_EVENTS 10

#define Es 1

using namespace std;

/**
 * Coding simulation result
 */
struct Result {
    uint64_t n;
    uint64_t block_errors;
    uint64_t bit_errors;
};

/**
 * Convolutional code description
 */
struct Code {
    int n; // Channel symbols per information symbol
    int m; // Constraint length of code
    int l; // Block length
    vector<uint32_t> gs; // Generator polynomials
};

/**
 * Reverse bit order of polynomial
 */
uint32_t rev(uint32_t k, int m) {
    uint32_t r = 0;
    for (int i = 0; i < m+1; i++) {
        r = (r << 1) | ((k >> i) & 1);
    }
    return r;
}

vector<float> modulate_bpsk(vector<uint8_t> input) {
    vector<float> symbols = {-1.0, 1.0};
    vector<float> out;
    for (uint32_t i = 0; i < input.size(); i++) {
        out.push_back(symbols[input[i]]);
    }
    return out;
}

Result simulate_communication(AWGNConvolutionalCode cc, int inf_len, int N, int MAX_N, int min_events, float snr) {
    stringstream result;
    AWGN channel(snr, inf_len/(float)(cc.inv_rate * (inf_len + cc.m)));
    vector<uint8_t> inf(inf_len, 0); // Zero inf vector
    vector<uint8_t> codeword = cc.encode(inf);
    vector<float> modulated = modulate_bpsk(codeword);

    uint64_t i = 0;
    uint64_t block_errors = 0;
    uint64_t bit_errors = 0;

    for (; i < MAX_N && (i < N || block_errors < min_events || (i-block_errors) < min_events); i++) {
        int errors = 0;
        vector<float> received = channel.transmit(modulated);
        vector<uint8_t> decoded = cc.decode(received, snr);
        
        for (int j = 0; j < inf_len; j++) {
            if (decoded[j] != inf[j]) {
                errors++;
            }
        }
        if (errors > 0) block_errors++;
        bit_errors += errors;
    }
    return {i, block_errors, bit_errors};
}

void simulate_threaded_communication(Code c, int N, int threads) {
    AWGNConvolutionalCode cc = AWGNConvolutionalCode(c.n, c.m, c.gs);
    
    string filename = "results/" + to_string(c.n) + "_1_" + to_string(c.m) + "_l" + to_string(c.l) + ".csv";
    stringstream results;
    results << "len,snr,n,block_errors,bit_errors\n";
    float snr_init = 2.5;
    bool rising = true;
    float snr = snr_init;
    while (true) {
        vector<future<Result>> futures;
        for (int i = 0; i < threads; i++)
            futures.push_back(async(
                        simulate_communication, 
                        cc, 
                        c.l, 
                        N / threads, 
                        MAX_SIMS / threads, 
                        (MIN_EVENTS + threads - 1) / threads, 
                        snr
            ));

        Result totals = {0,0,0};
        for (auto &fut : futures) {
            Result r = fut.get();
            totals.n += r.n;
            totals.block_errors += r.block_errors;
            totals.bit_errors += r.bit_errors;
        }
        if (totals.block_errors >= MIN_EVENTS) {
            results 
                << c.l
                << "," << fixed << setprecision(1) << snr 
                << "," << totals.n 
                << "," << totals.block_errors 
                << "," << totals.bit_errors 
                << endl;
        }
        printf("(%d,%d,%d):\t%.1f,%ld,%ld,%ld\n", c.n, 1, c.m, snr, totals.n, totals.block_errors, totals.bit_errors);
        snr += rising ? 0.5 : -0.5;
        if (totals.block_errors < MIN_EVENTS) {
            rising = false;
            snr = snr_init - 0.5;
        } 
        else if (totals.block_errors > totals.n - MIN_EVENTS || snr < -10) {
            break;
        }
    }
    ofstream f;
    f.open(filename);
    f << results.str();
    f.close();
}

long double simulate_BSC_rova(BSCConvolutionalCode cc, Code c, int N, float bsc_p) {
    BSC channel(bsc_p);
    vector<uint8_t> inf(c.l,0);
    vector<uint8_t> codeword = cc.encode(inf);

    long double total = 0;
    for (int i = 0; i < N; i++) {
        vector<uint8_t> received = channel.transmit(codeword);
        total += cc.rova(received, bsc_p);
    }
    return total;
}

void simulate_threaded_BSC_min_entropy(Code c, int N, int threads) {
    BSCConvolutionalCode cc = BSCConvolutionalCode(c.n, c.m, c.gs);

    string filename = "results/min-entropy_BSC_" + to_string(c.n) + "_1_" + to_string(c.m) + "_l" + to_string(c.l) + ".csv";
    stringstream results;
    results << "p,H\n0.0,0.0\n";
    for (float p = 0.025; p < 0.501; p += 0.025) {
        vector<future<long double>> futures;
        for (int i = 0; i < threads; i++)
            futures.push_back(async(simulate_BSC_rova, cc, c, N / threads, p));

        long double total = 0;
        for (auto &fut : futures) {
            total += fut.get();
        }
        long double H_min = -log2(total / N);
        results 
            << fixed << setprecision(3) << p
            << "," << fixed << setprecision(10)  << H_min
            << endl;
        printf("(%d,%d,%d,%d):\t%.3f,%.10Lf\n", c.n, 1, c.m, c.l, p, H_min);
    }
    ofstream f;
    f.open(filename);
    f << results.str();
    f.close();
}

long double simulate_BSC_entropy(BSCConvolutionalCode cc, Code c, int N, float bsc_p) {
    BSC channel(bsc_p);
    vector<uint8_t> inf(c.l, 0); // Zero inf vector
    vector<uint8_t> codeword = cc.encode(inf);

    long double total = 0;
    for (int i = 0; i < N; i++) {
        vector<uint8_t> received = channel.transmit(codeword);
        total += cc.entropy(received, bsc_p);
    }
    return total;
}

void simulate_threaded_BSC_entropy(Code c, int N, int threads) {
    BSCConvolutionalCode cc = BSCConvolutionalCode(c.n, c.m, c.gs);

    string filename = "results/entropy_BSC_" + to_string(c.n) + "_1_" + to_string(c.m) + "_l" + to_string(c.l) + ".csv";
    stringstream results;
    results << "p,H\n0.0,0.0\n";
    for (float p = 0.025; p < 1; p += 0.025) {
        vector<future<long double>> futures;
        for (int i = 0; i < threads; i++)
            futures.push_back(async(simulate_BSC_entropy, cc, c, N / threads, p));

        long double total = 0;
        for (auto &fut : futures) {
            total += fut.get();
        }
        long double H = total / N;
        results 
            << fixed << setprecision(3) << p
            << "," << fixed << setprecision(10)  << H
            << endl;
        printf("(%d,%d,%d,%d):\t%.3f,%.10Lf\n", c.n, 1, c.m, c.l, p, H);
    }
    results << "1.0,0.0" << endl;
    ofstream f;
    f.open(filename);
    f << results.str();
    f.close();
}

long double simulate_AWGN_rova(AWGNConvolutionalCode cc, Code c, int N, float snr) {
    AWGN channel(snr, c.l/(float)(cc.inv_rate * (c.l + cc.m)));
    vector<uint8_t> inf(c.l, 0); // Zero inf vector
    vector<uint8_t> codeword = cc.encode(inf);
    vector<float> modulated = modulate_bpsk(codeword);

    long double total = 0;
    for (int i = 0; i < N; i++) {
        vector<float> received = channel.transmit(modulated);
        total += cc.rova(received, snr);
    }
    return total;
}

void simulate_threaded_AWGN_min_entropy(Code c, int N, int threads) {
    AWGNConvolutionalCode cc = AWGNConvolutionalCode(c.n, c.m, c.gs);
    
    string filename = "results/min-entropy_AWGN_" + to_string(c.n) + "_1_" + to_string(c.m) + "_l" + to_string(c.l) + ".csv";
    stringstream results;
    results << "snr,H_min\n";
    for (float snr = -20; snr <= 10; snr += 0.5) {
        vector<future<long double>> futures;
        for (int i = 0; i < threads; i++)
            futures.push_back(async(simulate_AWGN_rova, cc, c, N / threads, snr));

        long double total = 0;
        for (auto &fut : futures) {
            total += fut.get();
        }
        long double H_min = -log2(total / N);
        results 
            << fixed << setprecision(2) << snr
            << "," << fixed << setprecision(10)  << H_min
            << endl;
        printf("(%d,%d,%d,%d):\t%.3f,%.10Lf\n", c.n, 1, c.m, c.l, snr, H_min);
    }
    ofstream f;
    f.open(filename);
    f << results.str();
    f.close();
}

long double simulate_AWGN_entropy(AWGNConvolutionalCode cc, Code c, int N, float snr) {
    AWGN channel(snr, c.l/(float)(cc.inv_rate * (c.l + cc.m)));
    vector<uint8_t> inf(c.l, 0); // Zero inf vector
    vector<uint8_t> codeword = cc.encode(inf);
    vector<float> modulated = modulate_bpsk(codeword);

    long double total = 0;
    for (int i = 0; i < N; i++) {
        vector<float> received = channel.transmit(modulated);
        total += cc.entropy(received, snr);
    }
    return total;
}

void simulate_threaded_AWGN_entropy(Code c, int N, int threads) {
    AWGNConvolutionalCode cc = AWGNConvolutionalCode(c.n, c.m, c.gs);
    
    string filename = "results/entropy_AWGN_" + to_string(c.n) + "_1_" + to_string(c.m) + "_l" + to_string(c.l) + ".csv";
    stringstream results;
    results << "snr,H\n";
    for (float snr = -20; snr <= 10; snr += 0.5) {
        if (snr > 8)
            N = 10000000;
        vector<future<long double>> futures;
        for (int i = 0; i < threads; i++)
            futures.push_back(async(simulate_AWGN_entropy, cc, c, N / threads, snr));

        long double total = 0;
        for (auto &fut : futures) {
            total += fut.get();
        }
        long double H = total / N;
        results 
            << fixed << setprecision(2) << snr
            << "," << fixed << setprecision(10)  << H
            << endl;
        printf("(%d,%d,%d,%d):\t%.3f,%.10Lf\n", c.n, 1, c.m, c.l, snr, H);
    }
    ofstream f;
    f.open(filename);
    f << results.str();
    f.close();
}

float Q(float x) {
    return 0.5 * erfc(x / sqrt(2));
}

void simulate_threaded_AWGN_hard_decision_entropy(Code c, int N, int threads) {
    BSCConvolutionalCode cc = BSCConvolutionalCode(c.n, c.m, c.gs);

    float coderate = c.l/(float)(cc.inv_rate * (c.l + cc.m));

    string filename = "results/entropy_AWGN_hard_dec_" + to_string(c.n) + "_1_" + to_string(c.m) + "_l" + to_string(c.l) + ".csv";
    stringstream results;
    results << "snr,H\n";
    for (float snr = -20; snr <= 10; snr += 0.5) {
        float snr_lin = pow(10, snr / 10.0);
        float Eb = Es / coderate;
        float N0 = Eb / snr_lin;
        float p = Q(sqrt(2*Es / N0)); // Lin & Costello eq. 1.4
        vector<future<long double>> futures;
        for (int i = 0; i < threads; i++)
            futures.push_back(async(simulate_BSC_entropy, cc, c, N / threads, p));

        long double total = 0;
        for (auto &fut : futures) {
            total += fut.get();
        }
        long double H = total / N;
        results 
            << fixed << setprecision(3) << snr
            << "," << fixed << setprecision(10)  << H
            << endl;
        printf("(%d,%d,%d,%d):\t%.2f,%.10Lf\n", c.n, 1, c.m, c.l, snr, H);
    }
    ofstream f;
    f.open(filename);
    f << results.str();
    f.close();
}

int main() {
    /* Seed the PRNG */
    srand (static_cast <unsigned> (time(0)));

    /* Simulation parameters */
    int N = 10000;
    int threads = 10;

    /* Codes from Lin & Costello */
    vector<Code> codes = {
        //{2, 2, 16,  {rev(05, 2),rev(07, 2)}},
        //{2, 4, 16,  {rev(027, 4),rev(031, 4)}},
        //{2, 6, 16,  {rev(0117, 6),rev(0155, 6)}},
        //{2, 2, 32,  {rev(05, 2),rev(07, 2)}},
        //{2, 4, 32,  {rev(027, 4),rev(031, 4)}},
        //{2, 6, 32,  {rev(0117, 6),rev(0155, 6)}},
        //{2, 2, 64,  {rev(05, 2),rev(07, 2)}},
        //{2, 4, 64,  {rev(027, 4),rev(031, 4)}},
        //{2, 6, 64,  {rev(0117, 6),rev(0155, 6)}},
        {2, 2, 128,  {rev(05, 2),rev(07, 2)}},
        {2, 4, 128,  {rev(027, 4),rev(031, 4)}},
        {2, 6, 128,  {rev(0117, 6),rev(0155, 6)}},
        //{3, 2, 128,  {rev(05, 2),rev(07, 2),rev(07, 2)}},
        //{3, 4, 128,  {rev(025, 4),rev(033, 4),rev(037, 4)}},
        //{3, 6, 128,  {rev(0117, 6),rev(0127, 6),rev(0155, 6)}},
        //{1, 1, 2,  {rev(03, 1)}},
    };

    for (Code c : codes) {
        printf("Simulating code (%d, 1, %d)\n", c.n, c.m);
        //simulate_threaded_BSC_entropy(c, N, threads);
        //simulate_threaded_AWGN_entropy(c, N, threads);
        //simulate_threaded_communication(c, N, threads);
        simulate_threaded_BSC_min_entropy(c, N, threads);
        //simulate_threaded_AWGN_hard_decision_entropy(c, N, threads);
    }

    printf("All simulations completed\n");
}
