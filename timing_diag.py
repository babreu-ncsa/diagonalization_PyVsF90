import numpy as np
from scipy.linalg import eigh
import time

its = 5
N = 200
samples = 100
t = its*[0]

for i in range(its):
    for j in range(samples):
        x = np.random.rand(N,N)
        s = time.time()
        w, v = eigh(x, lower=False)
        e = time.time()
        t[i] = t[i] + (e - s)
    N = N*2

print(t, sum(t))
