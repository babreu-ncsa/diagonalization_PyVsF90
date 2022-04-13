import numpy as np
import pandas as pd
from scipy.linalg import eigh
import time

its = 5
N = 200
samples = 100
t = np.zeros(samples)
df = pd.DataFrame()

for i in range(its):
    for j in range(samples):
        x = np.random.rand(N,N)
        s = time.perf_counter()
        w, v = eigh(x, lower=False)
        e = time.perf_counter()
        t[j] = (e - s)
    
    df[str(N)] = t
    N = N*2

df.to_csv("./py-diag-times.csv")
