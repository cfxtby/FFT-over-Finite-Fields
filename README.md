# FFT-over-Finite-Fields
Research project

In this project, we have implemented Shuhong Gao and Todd Mateer's algorithm for computing the FFT of a polynomial over a finite field of size 2^m,  where m=2^k. The algorithm involves traversing a n x n matrix column-wise at first and then row-wise, where n=2^m. We show that, as n grows, it is more efficient to transpose the matrix before the traversal operation, using a cache-oblivious matrix transposition algorithm discussed by Frigo et al.

We have also designed and implemented a parallel algorithm that computes the FFT and investigated the speed-up achieved. We showed that our parallel implementation performs well for large n and described its performance on two distributed systems, a High Performance Computing (HPC) cluster and a personal network of computers.
