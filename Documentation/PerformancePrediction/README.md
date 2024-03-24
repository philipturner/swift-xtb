# Performance Prediction

This experiment will predict the performance of matrix diagonalization with different algorithms and hardware. There is a Swift script, which takes different parameters like block size, computer power, parallelism, etc. If written correctly, it should reproduce the performance behavior of LAPACK with one-stage tridiagonalization.

> The divide-and-conquer part will be omitted from the performance prediction. It is out of scope, and not the bottleneck being optimized.
