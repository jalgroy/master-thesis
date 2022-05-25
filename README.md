## Master Thesis

This repository contains code used in my master thesis _Secure, reliable, and efficient communication over the wiretap channel_.

The code is used for simulation and analysis, and is not intended to be deployed in actual communication systems. The focus has been on correctness, not efficiency. For instance, techniques such as bit quantization, lookup tables, simpler Viterbi metrics, +++ should probably be used for a real-world deployment.

For Turbo Codes, custom scenarios for the "Coded Modulation Library" (http://iterativesolutions.com/Matlab.htm) are used.

User friendliness has not been a focus either. Some scripts may have hard-coded SNR values that have been read from simulation results. Many programs will have function calls that are commented out in main(). Simulation iterations, channel ranges, threads, etc. are also hard coded. 

Either way, you are free to use this any way you want.
