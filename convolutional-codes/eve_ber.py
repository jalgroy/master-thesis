#!/usr/bin/env python3

import csv

'''Calculate BER and FER as a function of the difference in SNR to someone with a given FER.
'''

codes = [
    "2_1_2_l128", 
    "2_1_4_l128", 
    "2_1_6_l128", 
]

bob_snr_10 = [
    3.5,
    2.65,
    2.05
] # Tuned for ~0.1 frame error rate

bob_snr_100 = [
    4.9,
    3.9,
    3.1
] # Tuned for ~0.01 frame error rate


bob_snr_1000 = [
    5.9,
    4.9,
    3.95
] # Tuned for ~0.001 frame error rate


types = ["ber","fer"]

for t in types:
    for i in range(len(codes)):
        with open(f'results/{t}_{codes[i]}.csv', 'r', newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            result = [['snr',t]]
            for row in reader:
                if row[0] == 'snr':
                    continue
                result.append([str(bob_snr_10[i]-float(row[0])), row[1]])
            with open(f'results/tuned_0.1_{t}_{codes[i]}.csv', 'w') as out:
                writer = csv.writer(out, delimiter=',')
                writer.writerows(result)

