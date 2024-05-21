
# Examples of Modular Multiplications

- `Montgomery_acc`: Accumulative variant of Montgomery multiplication.
    - Assumed knowledge: Integer arithmetic.
    - References: [Mon85], [ACC+21], [GKS21]
    - Additional references:
- `Montgomery_sub`: Subtractive variant of Montgomery multiplication.
    - Assumed knowledge: Integer arithmetic.
    - References: [Seil18]
    - Additional references:
- `Barrett`: Barrett multiplication.
    - Assumed knowledge: Integer arithmetic.
    - References: [BHK+22]
    - Additional references:
- `Barrett_approx`: Barrett multiplication with generalized approximation.
    - Assumed knowledge: Integer arithmetic.
    - References: [HKS23]
    - Additional references:
- `Barrett_Montgomery_cmp`: Barrett multiplication.
    - Assumed knowledge: Integer arithmetic.
    - References: [BHK+22]
    - Additional references:

# References

[ACC+21]
Erdem Alkim, Dean Yun-Li Cheng, Chi-Ming Marvin Chung, Hülya Evkan, Leo Wei-Lun Huang, Vincent Hwang, Ching-Lin Trista Li, Ruben Niederhagen, Cheng-Jhih Shih, Julian Wälde, and Bo-Yin Yang. Polynomial Multiplication in NTRU Prime Comparison of Optimization Strategies on Cortex-M4. IACR Transactions on Cryptographic Hardware and Embedded Systems, 2021(1):217–238, 2021. https://tches.iacr.org/index.php/TCHES/article/view/8733.

[BHK+22]
Hanno Becker, Vincent Hwang, Matthias J. Kannwischer, Bo-Yin Yang, and Shang-Yi Yang. Neon NTT: Faster Dilithium, Kyber, and Saber on Cortex-A72 and Apple M1. IACR Transactions on Cryptographic Hardware and Embedded Systems, 2022(1):221–244, 2022. https://tches.iacr.org/index.php/TCHES/article/view/9295.

[GKS21]
Denisa O. C. Greconici, Matthias J. Kannwischer, and Daan Sprenkels. Compact Dilithium Implementations on Cortex-M3 and Cortex-M4. IACR Transactions on Cryptographic Hardware and Embedded Systems, 2021(1):1–24, 2021. https://tches.iacr.org/index.php/TCHES/article/view/8725.

[HKS23]
Vincent Hwang, YoungBeom Kim, and Seog Chung Seo. Barrett Multiplication for Dilithium on Embedded Devices. 2023. https://eprint.iacr.org/2023/1955.

[Mon85]
Peter L. Montgomery. Modular Multiplication Without Trial Division. Mathematics of computation, 44(170):519–521, 1985. https://www.ams.org/journals/mcom/1985-44-170/S0025-5718-1985-0777282-X/?active=current.

[Seil18]
Gregor Seiler. Faster AVX2 optimized NTT multiplication for Ring-LWE lattice cryptography. 2018. https://eprint.iacr.org/2018/039.





