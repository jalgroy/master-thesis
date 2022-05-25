#!/usr/bin/env python3

import numpy as np
import galois

def gaussian_elimination(m):
    lead = 0
    rc = len(m)
    cc = len(m[0])
    for r in range(rc):
        if cc <= lead:
            break
        i = r
        while m[i][lead] == 0:
            i += 1
            if rc == i:
                i = r
                lead += 1
                if cc == lead:
                    return
        # Swap row i and r
        if i != r:
            ri = m[i]
            m[i] = m[r]
            m[r] = ri

        for i in range(rc):
            if i != r and m[i][lead] == 1:
                for j in range(cc):
                    m[i][j] ^= m[r][j]
        lead += 1

def pad(l, n):
    while len(l) < n:
        l.append(0)

GF2 = galois.GF(2)

#n = 63
#k = 36
#generator = "1033500423"
#n = 31
#k = 16
#generator = "107657"
#n = 127
#k = 64
#generator = "1206534025570773100045"
n = 255
k = 131
generator = "215713331471510151261250277442142024165471"

bits = ''.join([bin(int(i))[2:].zfill(3) for i in generator])

bits = [int(i) for i in bits]
while bits[0] == 0:
    bits.pop(0)

g = galois.Poly(bits, field=GF2)


f = galois.Poly([1] + [0]*(n-1) + [1], field=GF2)

h = f / g

g_bits = g.coeffs.tolist()
g_bits.reverse()
G = []
for i in range(k):
    row = [0]*i
    row += g_bits
    row += [0]*(n-i-len(g_bits))
    G.append(row)

h_bits = h.coeffs.tolist()
h_bits.reverse()
H = []
for i in range(n-k):
    row = [0]*(n-i-len(h_bits))
    row += h_bits
    row += [0]*i
    H.append(row)

# Extend to create 64,36 code
for i in range(k):
    G[i].append(sum(G[i]) % 2)

for i in range(n-k):
    H[i].append(0)
H.append([1]*(n+1))

# Delete 4 rows to get 64,32 subcode
#G = G[:-4]

# Delete 3 rows to get 256,128 subcode
G = G[:-3]

n = len(G[0])
k = len(G)


gaussian_elimination(G)

print(f'G: {k} x {n}')
for row in G:
    print("{" + ",".join([str(r) for r in row]) + "},")

H = [[0]*n for _ in range(n-k)]

for i in range(n-k):
    for j in range(k):
        H[i][j] = G[j][k+i]

for i in range(n-k):
    H[i][k+i] = 1

print(f'H: {n-k} x {n}')
for row in H:
    print("{" + ",".join([str(r) for r in row]) + "},")

