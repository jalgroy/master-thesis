#ifndef CONV_CODE_H
#define CONV_CODE_H

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>

using namespace std;

template <typename T>
class ConvolutionalCode {
    private:
        /**
         * Edge transition probability
         *
         * @param[in]   ch_symbols      Channel symbols
         * @param[in]   edge_symbols    Edge symbols
         * @param[in]   channel_quality Channel quality (e.g. SNR or error probability, depending on the type of channel)
         *
         * @return      P(ch_symbols | edge_symbols)
         */
        virtual long double edge_p(vector<T> ch_symbols, vector<uint8_t> edge_symbols, float channel_quality) = 0;

        /**
         * Viterbi metric
         *
         * @param[in]   ch_symbols      Channel symbols
         * @param[in]   edge_symbols    Edge symbols
         * @param[in]   channel_quality Channel quality (e.g. SNR or error probability, depending on the type of channel)
         *
         * @return      -log(P(ch_symbols | edge_symbols))
         */
        long double metric(vector<T> ch_symbols, vector<uint8_t> edge_symbols, float channel_quality) {
            return -log2(edge_p(ch_symbols, edge_symbols, channel_quality));
        }

        /**
         * Binary entropy
         *
         * @param[in]   p       Probability
         *
         * @return      h(p)
         */
        long double h(long double p) {
            if (p == 0.0 || p == 1.0) return 0;
            return -p*log2(p)-(1-p)*log2(1-p);
        }
        
        /**
         * Get trellis output for a given state and an input symbol
         *
         * @param[in]   state   Trellis state
         * @param[in]   input   Input symbol
         *
         * @return      Output symbols
         */
        vector<uint8_t> get_output(uint32_t state, uint8_t input) {
            vector<uint8_t> output;
            for (uint32_t j = 0; j < inv_rate; j++) {
                uint8_t g = ((generators[j] >> m) & 1) & input;
                for (int k = m-1; k >= 0; k--) {
                    g ^= ((generators[j] >> k) & 1) & ((state >> (m-1-k)) & 1);
                }
                output.push_back(g);
            }
            return output;
        }

        /**
         * Trellis node for Viterbi
         */
        struct Node {
            /* State */
            uint32_t state;
            /* Metric */
            long double m;
            /* Backpointer */
            long double bp;
        };

        /**
         * Trellis node for entropy calculation
         */
        struct HNode {
            /* State */
            uint32_t state;
            /* Forward probability */
            long double p;
            /* Backward probability */
            long double q;
            /* Binary entropy */
            long double h;
            /* State probability pi */
            long double pi;
        };

        /**
         * Trellis node for ROVA
         */
        struct RNode {
            /* State */
            uint32_t state;
            /* Metric */
            long double m;
            /* Backpointer */
            long double bp;
            /* P(survivor path correct) */
            long double P;
            /* P(non-survivor path correct) */
            long double PN;
        };
    public:
        uint32_t inv_rate;
        uint32_t m;
        vector<uint32_t> generators;

        /**
         * Convolutional code constructor
         *
         * @param[in]   inv_rate        Inverse of the code rate. I.e. 3 if the code rate is 1/3
         * @param[in]   m               Number of state bits. Trellis will have 2^m states per depth.
         * @param[in]   generators      Vector of m generator polynomials for the encoder.
         */
        ConvolutionalCode(uint32_t inv_rate, uint32_t m, vector<uint32_t> generators) {
            this->inv_rate = inv_rate;
            this->m = m;
            this->generators = generators;
        }
        
        /**
         * Encode information vector using trellis.
         *
         * @param[in]   inf_vector      Information bit vector
         *
         * @return      Vector of encoded bits
         */
        vector<uint8_t> encode(vector<uint8_t> inf_vector) {
            /* Pad information vector with m zero bits */
            for (uint32_t i = 0; i < m; i++) inf_vector.push_back(0);
            
            /* Start at state zero */
            uint32_t state = 0;

            vector<uint8_t> output;
            for (uint32_t i = 0; i < inf_vector.size(); i++) {
                vector<uint8_t> out = get_output(state, inf_vector[i]);
                output.insert(output.end(), out.begin(), out.end());
                state = (state << 1) | (inf_vector[i] & 1);
            }
            return output;
        }

        /**
         * Decode vector from channel. Hard decisions.
         *
         * @param[in]   ch_vector       Vector from channel
         * @param[in]   channel_quality Channel quality (e.g. SNR or error probability, depending on the channel)
         *
         * @return      Decoded information vector
         */
        vector<uint8_t> decode(vector<T> ch_vector, float channel_quality) {
            /* Construct trellis */
            int n = ch_vector.size() / inv_rate;
            vector<Node> trellis[n+1];
            for (int i = 0; i <= n; i++) {
                for (int s = 0; s < 1 << min(i, (int)m); s++) {
                    Node n;
                    n.state = s;
                    if (i == 0) n.m = 0;
                    trellis[i].push_back(n);
                }
            }
            /* --- Viterbi decoding algorithm --- */
            for (int i = 0; i < n; i++) {
                for (uint32_t s = 0; s < trellis[i+1].size(); s++) {
                    float ms = INFINITY;
                    for (uint32_t j = 0; j < 2; j++) {
                        uint32_t sp = (s >> 1) | (j << (m-1)); // Previous states with edges to current state
                        if (sp >= trellis[i].size()) continue;
                        vector<T> ch_symbols(&ch_vector[i*inv_rate], &ch_vector[i*inv_rate+inv_rate]);
                        vector<uint8_t> edge_symbols = get_output(sp, s & 1);
                        /* DEBUG */
                        float m = trellis[i][sp].m + metric(ch_symbols, edge_symbols, channel_quality);
                        if (m <= ms) {
                            ms = m;
                            trellis[i+1][s].bp = sp;
                        }
                    }
                    trellis[i+1][s].m = ms;
                }
            }
            /* Backtrack to collect decoded vector */
            vector<uint8_t> out;
            uint32_t bp = trellis[n][0].bp;
            out.insert(out.begin(), 0);
            for (int i = n-1; i > 0; i--) {
                out.insert(out.begin(), (uint8_t)(trellis[i][bp].state & 1));
                bp = trellis[i][bp].bp;
            }
            /* Remove padding */
            for (uint32_t i = 0; i < m; i++) out.pop_back();

            return out;
        }
        
        /**
         * Calculate conditional entropy H(X | R = r)
         *
         * @param[in]   r               Received channel vector 
         * @param[in]   channel_quality Channel quality (e.g. SNR or error probability, depending on the type of channel)
         *
         * @return      H(X | R = r)
         */
        long double entropy(vector<T> r, float channel_quality) {
            /* --- Construct trellis --- */
            int n = r.size() / inv_rate;
            vector<HNode> trellis[n+1];
            for (int i = 0; i <= n; i++) {
                for (int s = 0; s < 1 << min(i, (int)m); s++) {
                    HNode node;
                    node.state = s;
                    if (i == 0) {
                        node.p = 1;
                        node.h = 0;
                    }
                    if (i == n) node.q = 1;
                    trellis[i].push_back(node);
                }
            }
            /* --- Forward pass --- */
            for (int i = 0; i < n; i++) {
                for (uint32_t s = 0; s < trellis[i+1].size(); s++) {
                    if (trellis[i].size() < trellis[i+1].size()) {
                        /* Less than two incoming edges per state */
                        uint32_t sp = (s >> 1); // Previous state with edge to current state

                        vector<T> rl(&r[i*inv_rate], &r[i*inv_rate+inv_rate]);
                        vector<uint8_t> edge_symbols = get_output(sp, s & 1);

                        trellis[i+1][s].p = trellis[i][sp].p * edge_p(rl, edge_symbols, channel_quality);
                        trellis[i+1][s].h = 0;
                    }
                    else {
                        long double p[] = {0,0};
                        for (uint32_t j = 0; j < 2; j++) {
                            uint32_t sp = (s >> 1) | (j << (m-1)); // Previous states with edges to current state
                            vector<T> rl(&r[i*inv_rate], &r[i*inv_rate+inv_rate]);
                            vector<uint8_t> edge_symbols = get_output(sp, s & 1);

                            p[j] = trellis[i][sp].p * edge_p(rl, edge_symbols, channel_quality);
                        }
                        trellis[i+1][s].p = p[0] + p[1];
                        trellis[i+1][s].h = h(p[0] / (p[0] + p[1]));
                    }
                }
            }
            /* --- Backward pass --- */
            for (int i = n-1; i >= 0; i--) {
                long double t = 0;
                for (uint32_t s = 0; s < trellis[i].size(); s++) {
                    long double q[2];
                    for (uint32_t j = 0; j < 2; j++) {
                        uint32_t sn = ((s << 1) & ((1 << m)-1)) | j; // Next states with edges from current states
                        vector<T> rl(&r[i*inv_rate], &r[i*inv_rate+inv_rate]);
                        vector<uint8_t> edge_symbols = get_output(s, sn & 1);
                        q[j] = trellis[i+1][sn].q * edge_p(rl, edge_symbols, channel_quality);
                    }
                    trellis[i][s].q = q[0] + q[1];
                }
                for (uint32_t s = 0; s < trellis[i+1].size(); s++) {
                    trellis[i+1][s].pi = (trellis[i+1][s].p * trellis[i+1][s].q) / t;
                }
            }
            /* --- Normalization -- */
            for (int i = 0; i <= n; i++) {
                long double t = 0;
                for (uint32_t s = 0; s < trellis[i].size(); s++) {
                    t += trellis[i][s].p * trellis[i][s].q;
                }
                for (uint32_t s = 0; s < trellis[i].size(); s++) {
                    trellis[i][s].pi = (trellis[i][s].p * trellis[i][s].q) / t;
                }
            }
            long double H = 0;
            for (int i = 1; i <= n; i++) {
                if (trellis[i-1].size() >= trellis[i].size()) {
                    for (uint32_t s = 0; s < trellis[i].size(); s++) {
                        H += trellis[i][s].h * trellis[i][s].pi;
                    }
                }
            }
            return H;
        }


        /**
         * Reliability Output Viterbi Algorithm (ROVA)
         *
         * @param[in]   r               Vector from channel
         * @param[in]   channel_quality Channel quality (e.g. SNR or error probability, depending on the channel)
         *
         * @return      Max over all x of p(x|y)
         */
        long double rova(vector<T> r, float channel_quality) {
            /* Construct trellis */
            int n = r.size() / inv_rate;
            vector<RNode> trellis[n+1];
            for (int i = 0; i <= n; i++) {
                for (int s = 0; s < 1 << min(i, (int)m); s++) {
                    RNode n;
                    n.state = s;
                    if (i == 0) {
                        n.P = 1;
                        n.PN = 0;
                        n.m = 0;
                    }
                    trellis[i].push_back(n);
                }
            }
            /* --- Viterbi decoding algorithm --- */
            for (int i = 0; i < n; i++) {
                for (uint32_t s = 0; s < trellis[i+1].size(); s++) {
                    float ms = INFINITY;
                    for (uint32_t j = 0; j < 2; j++) {
                        uint32_t sp = (s >> 1) | (j << (m-1)); // Previous states with edges to current state
                        if (sp >= trellis[i].size()) continue;
                        vector<T> ch_symbols(&r[i*inv_rate], &r[i*inv_rate+inv_rate]);
                        vector<uint8_t> edge_symbols = get_output(sp, s & 1);
                        float m = trellis[i][sp].m + metric(ch_symbols, edge_symbols, channel_quality);
                        if (m <= ms) {
                            ms = m;
                            trellis[i+1][s].bp = sp;
                        }
                    }
                    trellis[i+1][s].m = ms;
                }
            }
            /* --- ROVA from Raghavan and Baum (1998)--- */
            for (int i = 1; i <= n; i++) {
                long double delta = 0;
                for (uint32_t s = 0; s < trellis[i].size(); s++) {
                    if (i > n-m && (s & ((1 << (m-(n-i))) - 1)) != 0) continue; // Bypass unused states at the tail of the trellis.
                    for (int j = 0; j < 2; j++) {
                        uint32_t sp = (s >> 1) | (j << (m-1)); // Previous states with edges to current state
                        if (sp >= trellis[i-1].size()) continue;
                        vector<T> rl(&r[(i-1)*inv_rate], &r[i*inv_rate]);
                        vector<uint8_t> edge_symbols = get_output(sp, s & 1);
                        delta += edge_p(rl, edge_symbols, channel_quality)*(trellis[i-1][sp].P + trellis[i-1][sp].PN);
                    }

                }
                for (uint32_t s = 0; s < trellis[i].size(); s++) {
                    if (i > n-m && (s & ((1 << (m-(n-i))) - 1)) != 0) continue; // Bypass unused states at the tail of the trellis.
                    vector<T> rl(&r[(i-1)*inv_rate], &r[i*inv_rate]);
                    vector<uint8_t> edge_symbols = get_output(trellis[i][s].bp, s & 1); // Edge from previous survivor to current state.
                    long double ep = edge_p(rl, edge_symbols, channel_quality);
                    trellis[i][s].P = (1.0 / delta) * ep * trellis[i-1][trellis[i][s].bp].P;
                    long double sum = 0;
                    for (int j = 0; j < 2; j++) {
                        uint32_t sp = (s >> 1) | (j << (m-1)); // Previous states with edges to current state
                        if (sp >= trellis[i-1].size()) continue;
                        edge_symbols = get_output(sp, s & 1);
                        sum += edge_p(rl, edge_symbols, channel_quality)*(trellis[i-1][sp].P + trellis[i-1][sp].PN);
                    }
                    trellis[i][s].PN = (1.0/delta) * sum - trellis[i][s].P;
                }
            }
            return trellis[n][0].P;
        }
};

class BSCConvolutionalCode : public ConvolutionalCode<uint8_t> {
    private:
        /* Probability function for the BSC */
        long double edge_p(vector<uint8_t> ch_symbols, vector<uint8_t> edge_symbols, float p) override {
            long double edge_p = 1;
            for (uint32_t i = 0; i < ch_symbols.size(); i++) {
                edge_p *= ch_symbols[i] == edge_symbols[i] ? (1-p) : p;
            }
            return edge_p;
        }
    public:
        BSCConvolutionalCode(uint32_t inv_rate, uint32_t m, vector<uint32_t> generators) 
            : ConvolutionalCode<uint8_t>(inv_rate, m, generators) {}
};

class AWGNConvolutionalCode : public ConvolutionalCode<float> {
    private:
        /**
         * AWGN edge probability
         *
         * p(rl | vl) from Lin & Costello eq. 12.13
         */
        long double edge_p(vector<float> ch_symbols, vector<uint8_t> edge_symbols, float snr_dB) override {
            long double snr = pow(10, snr_dB/10.0)/inv_rate;
            long double Es = 1.0;
            long double N0 = 0.5 / snr;
            long double edge_p = 1;

            for (int i = 0; i < inv_rate; i++) {
                float x = edge_symbols[i] == 0 ? -1.0 : 1.0;
                edge_p *= sqrt(Es / (M_PI * N0)) * pow(M_E, -(Es/N0)*pow(ch_symbols[i] - x, 2));
            }
            
            return edge_p;
        }
    public:
        AWGNConvolutionalCode(uint32_t inv_rate, uint32_t m, vector<uint32_t> generators) 
            : ConvolutionalCode<float>(inv_rate, m, generators) {}
};

#endif
