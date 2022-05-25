#ifndef SYND_CODE_H
#define SYND_CODE_H

#include "util.hpp"

class SyndromeCode {
    private:
        uint8_t **G;
        uint8_t **H;
        uint8_t **HT;
        vector<uint8_t> zero_noise_v;

        /**
         * Gauss-Jordan elimination of matrix with elements in GF(2). In-place.
         * 
         * Beware: Does _not_ handle other fields!
         */
        void gauss_jordan(uint8_t **mtx, size_t m, size_t n) {
            int lead = 0;
            for (int r = 0; r < m; r++) {
                if (n <= lead) break;
                int i = r;
                while (mtx[i][lead] == 0) {
                    if (m == ++i) {
                        i = r;
                        if (n == ++lead) return;
                    }
                }
                // Swap i and r
                if (i != r) {
                    uint8_t *ri = mtx[i];
                    mtx[i] = mtx[r];
                    mtx[r] = ri;
                }

                for (i = 0; i < m; i++) {
                    if (i != r && mtx[i][lead] == 1) {
                        for (int j = 0; j < n; j++) {
                            mtx[i][j] ^= mtx[r][j];
                        }
                    }
                }
                lead++;

            }
        }

        /**
         * Get a noise vector corresponding to a given syndrome by
         * choosing a valid solution to the equation Hx = s.
         *
         * @param[in]   mtx     Reduced echelon form of H | s
         */
        vector<uint8_t> get_noise_v(uint8_t **mtx) {
            vector<uint8_t> noise_v(n, 0);
            int pivot = -1;
            for (int i = 0; i < n-k; i++) {
                for (int j = 0; j < n; j++) {
                    if (mtx[i][j] == 1) {
                        if (pivot == -1) {
                            pivot = j;
                        }
                        else if (pivot != j-1) {
                            for (int k = 1; k < j-pivot; k++) {
                                /* Set free variable */
                                noise_v[j-k] = MWC & 1;
                            }
                        }
                        pivot = j;
                        break;
                    }
                }
            }
            pivot = -1;
            for (int i = 0; i < n-k; i++) {
                for (int j = 0; j < n; j++) {
                    if (mtx[i][j] == 1) {
                        pivot = j;
                        uint8_t xj = 0;
                        for (int k = j+1; k < n; k++) {
                            xj ^= mtx[i][k] & noise_v[k];
                        }
                        noise_v[j] = xj ^ mtx[i][n];
                        break;
                    }
                }
            }
            return noise_v;
        }

        /**
         * Get a random codeword
         */
        vector<uint8_t> random_codeword() {
            vector<uint8_t> inf;
            vector<uint8_t> cw;
            for (int i = 0; i < k; i++) inf.push_back(MWC & 1);
            for (int i = 0; i < n; i++) {
                uint8_t sum = 0;
                for (int j = 0; j < k; j++) {
                    sum ^= inf[j] & G[j][i];
                }
                cw.push_back(sum);
            }
            return cw;
        }


    public:
        const int n, k;

        /**
         * Constructor
         */
        SyndromeCode(int n, int k, uint8_t **G_mtx, uint8_t **H_mtx) :
                G(new uint8_t*[k]), 
                H(new uint8_t*[n-k]), 
                HT(new uint8_t*[n]),
                n(n), k(k)
        {
            init_MWC();
            for (int i = 0; i < k; i++) G[i] = new uint8_t[n];
            for (int i = 0; i < n-k; i++) H[i] = new uint8_t[n];
            for (int i = 0; i < n; i++) HT[i] = new uint8_t[n-k];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n-k; j++) {
                    H[j][i] = H_mtx[j][i];
                    HT[i][j] = H_mtx[j][i];
                }
                for (int j = 0; j < k; j++) {
                    G[j][i] = G_mtx[j][i];
                }
            }

            /* Init zero vector */
            vector<uint8_t> inf_vector(k);
            for (int i = 0; i < k; i++) inf_vector[i] = 0;

            /* Allocate matrix H | s */
            uint8_t **mtx = new uint8_t*[n-k];
            for (int i = 0; i < n-k; i++)
                mtx[i] = new uint8_t[n+1];

            /* Fill matrix and do gauss-jordan elimination */
            for (int i = 0; i < n-k; i++)
                for (int j = 0; j < n; j++)
                    mtx[i][j] = H[i][j];
            for (int i = 0; i < n-k; i++)
                mtx[i][n] = inf_vector[i];
            gauss_jordan(mtx, n-k, n+1);

            /* Get noise vector and random codeword */
            zero_noise_v = get_noise_v(mtx);

            /* Delete allocated matrix */
            for (int i = 0; i < n-k; i++) {
                delete[] mtx[i];
            }
            delete[] mtx;
        }

        SyndromeCode(const SyndromeCode &obj) :
                G(new uint8_t*[obj.k]), 
                H(new uint8_t*[obj.n-obj.k]), 
                HT(new uint8_t*[obj.n]),
                zero_noise_v(obj.zero_noise_v),
                n(obj.n), k(obj.k)
        {
            for (int i = 0; i < k; i++) G[i] = new uint8_t[n];
            for (int i = 0; i < n-k; i++) H[i] = new uint8_t[n];
            for (int i = 0; i < n; i++) HT[i] = new uint8_t[n-k];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n-k; j++) {
                    H[j][i] = obj.H[j][i];
                    HT[i][j] = obj.HT[i][j];
                }
                for (int j = 0; j < k; j++) {
                    G[j][i] = obj.G[j][i];
                }
            }
        }

        /**
         * Destructor
         */
        ~SyndromeCode() {
            for (int i = 0; i < k; i++) {
                delete[] G[i];
            }
            for (int i = 0; i < n-k; i++) {
                delete[] H[i];
            }
            for (int i = 0; i < n; i++) {
                delete[] HT[i];
            }
            delete[] G;
            delete[] H;
            delete[] HT;
        }

        /**
         * Decode channel vector by calculating its syndrome
         */
        vector<uint8_t> decode(vector<float> ch_vector) {
            if (ch_vector.size() != n) throw;
            vector<uint8_t> ch_bits;
            for (int i = 0; i < n; i++)
                ch_bits.push_back(ch_vector[i] < 0 ? 0 : 1);

            vector<uint8_t> decoded;
            for (int i = 0; i < n-k; i++) {
                uint8_t res = 0;
                for (int j = 0; j < n; j++) {
                    res ^= ch_bits[j] & HT[j][i];
                }
                decoded.push_back(res);
            }
            return decoded;
        }

        /**
         * Encode zero information vector
         */
        vector<uint8_t> encode_zero() {
            vector<uint8_t> codeword = random_codeword();

            /* Add noise to codeword */
            for (int i = 0; i < n; i++)
                codeword[i] = codeword[i] ^ zero_noise_v[i];
            return codeword;
        }


        /**
         * Encode information vector as channel vector whose syndrome is the information
         */
        vector<uint8_t> encode(vector<uint8_t> inf_vector) {
            if (inf_vector.size() != n - k) throw;

            /* Allocate matrix H | s */
            uint8_t **mtx = new uint8_t*[n-k];
            for (int i = 0; i < n-k; i++)
                mtx[i] = new uint8_t[n+1];

            /* Fill matrix and do gauss-jordan elimination */
            for (int i = 0; i < n-k; i++)
                for (int j = 0; j < n; j++)
                    mtx[i][j] = H[i][j];
            for (int i = 0; i < n-k; i++)
                mtx[i][n] = inf_vector[i];
            gauss_jordan(mtx, n-k, n+1);

            /* Get noise vector and random codeword */
            vector<uint8_t> noise_v = get_noise_v(mtx);
            vector<uint8_t> codeword = random_codeword();

            /* Add noise to codeword */
            for (int i = 0; i < n; i++)
                codeword[i] = codeword[i] ^ noise_v[i];

            /* Delete allocated matrix */
            for (int i = 0; i < n-k; i++) {
                delete[] mtx[i];
            }
            delete[] mtx;

            return codeword;
        }
};

#endif
