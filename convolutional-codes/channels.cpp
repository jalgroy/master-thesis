#include <vector>

#include "channels.hpp"

using namespace std;

/* BSC Constructor */
BSC::BSC(float p) {
    this->p = p;
    this->dist = bernoulli_distribution(p);
}

/* BSC transmit */
vector<uint8_t> BSC::transmit(vector<uint8_t> input) {
    vector<uint8_t> out;
    for (uint32_t i = 0; i < input.size(); i++) {
        out.push_back((dist(gen)) ? input[i] ^ 1 : input[i]);
    }
    return out;
}

/* AWGN Constructor */
AWGN::AWGN(float snr_dB, float coderate) {
    float Eb = 1.0;
    float snr = pow(10, snr_dB/10.0)*coderate;
    float N0 = Eb/snr;
    float std_dev = sqrt(N0/2);

    this->e2 = mt19937(this->rd());
    this->dist = normal_distribution<>(0, std_dev);
}

/* AWGN transmit */
vector<float> AWGN::transmit(vector<float> input) {
    vector<float> out;
    for (uint32_t i = 0; i < input.size(); i++) {
        out.push_back(input[i] + dist(e2));
    }
    return out;
}
