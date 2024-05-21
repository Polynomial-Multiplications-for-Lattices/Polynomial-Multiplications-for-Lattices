
# Polynomial-Multiplications-for-Lattices

This repository collects several tricks for polynomial multiplications that have found practical importance in lattice-based cryptosystems.
In the literature, researchers argued the efficiency of specific approaches multiplying polynomials on certain platforms.
Since hardware improves overtime, we eventually need to transfer the knowledge to future platforms.
There are two parts in this repository:
- A collection of C implementations demonstrating each of the tricks.
- A collection of implementations with platform-specific optimizations explaining how researchers optimized for the target platforms.

# Survey

[Ber01]
Daniel J. Bernstein. Multidigit multiplication for mathematicians. 2001. https://cr.yp.to/papers.html#m3.

[Ber08]
Daniel J. Bernstein. Fast multiplication and its applications. Algorithmic number theory, 44:325–384, 2008. https://cr.yp.to/papers.html#multapps.

[DV90]
Pierre Duhamel and Martin Vetterli. Fast Fourier transforms: a tutorial review and a state of the art. Signal processing, 19(4):259–299, 1990. https://www.sciencedirect.com/science/article/pii/016516849090158U.

[LZ22]
Zhichuang Liang and Yunlei Zhao. Number Theoretic Transform and Its Applications in Lattice-based Cryptosystems: A Survey. arXiv preprint arXiv:2211.13546, 2022. https://arxiv.org/abs/2211.13546.

[Nus82]
Henri J. Nussbaumer. Fast Fourier Transform and Convolution Algorithms. Springer Berlin, Heidelberg, 2nd edition, 1982. https://doi.org/10.1007/978-3-642-81897-4.

[Win80]
Shmuel Winograd. Arithmetic Complexity of Computations, volume 33. Society for Industrial and Applied Mathematics, 1980. https://epubs.siam.org/doi/10.1137/1.9781611970364.

# Call for Contributors

TBA

# TODOs

- See if we want to include Kronecker substitution. There are two incentives.
    - Paper [Post-Quantum Cryptography with Contemporary Co-Processors: Beyond Kronecker, Schönhage-Strassen & Nussbaumer](https://www.usenix.org/conference/usenixsecurity22/presentation/bos) reduces polynomial multiplications to large integer multiplications.
    - The program [here](https://github.com/mupq/pqm4/blob/Round3/crypto_kem/sntrup761/m4f/jump16divsteps_mod3_asm.S) shows how to multiply polynomials over Z_3 with the long multiplication instructions `umull` and `umlal`.




