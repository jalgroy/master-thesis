#include <iostream>
#include <vector>
#include <random>

using namespace std;

class BSC {
    public:
        float p;
        bernoulli_distribution dist;
        default_random_engine gen;

        BSC(float p);

        /**
         * Transmit bit vector over BSC
         *
         * @param[in]   input   Vector of bits to be transmitted
         *
         * @return      Vector of bits received from channel
         */
        vector<uint8_t> transmit(vector<uint8_t> input);
};

class AWGN {
    public:
        random_device rd;
        mt19937 e2; // Mersenne Twister
        normal_distribution<> dist;

        AWGN(float snr_dB, float coderate);

        /**
         * Transmit vector of real-valued symbols over AWGN
         * 
         *
         * @param[in]   input   Vector of symbols to be transmitted
         *
         * @return      Vector of bits received from channel
         */
        vector<float> transmit(vector<float> input);
};
