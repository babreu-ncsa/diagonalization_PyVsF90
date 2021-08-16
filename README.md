# diagonalziation_PyVsF90
Performance comparison between SciPy's linalg.eigh and LAPACK's DSYEV

### Notes
I ran these codes on AKD EPYC 7742 nodes. The Fortran version is (5.4 +- 2.8)x faster than the Python version for 100 NxN random matrices, with N from 200 to 3200.

- Fortran: compiled with Intel's ifort 2019.144.1, using MKL.
- Python: v3.8.5, SciPy v1.5.2

LAPACK function: http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen_ga442c43fca5493590f8f26cf42fed4044.html

SciPy's funtion: https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.eigh.html#scipy.linalg.eigh


### Additional work
Python is smart enough to know which diagonalization function to call. I have tested it against the LAPACK's DOUBLE PRECISION version of the diagonalization routine. If numpy.random.rand() populates the matrices with single precision numbers, than a fair comparsion would be with SSYEV, not DSYEV, and you can add another 2x factor of speedup, making Fortran almost 10x faster.

