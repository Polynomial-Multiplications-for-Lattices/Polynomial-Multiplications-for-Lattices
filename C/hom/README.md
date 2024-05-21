
# Examples of Homomorphisms

- `DWT.c`: This file demonstrates Cooley--Tukey FFT for discrete weighted transform.
    - Assumed knowledge: Chinese remainder theorem for polynomial rings.
    - References: [CT65], [CF94].
    - Additional references:
- `FNT.c`: This file demonstrates Fermat number transform.
    - Assumed knowledge: Chinese remainder theorem for polynomial rings.
    - References: [AB74].
    - Additional references: [SS71], [CT65].
- `GT.c`: This file demonstrates Good--Thomas FFT.
    - Assumed knowledge: Multi-variate polynomial rings (minimum); tensor product of associate algebras (recommended).
    - References: [Goo58].
    - Additional references:
- `Karatsuba.c`: This file demonstrates Karatsuba.
    - Assumed knowledge: Chinese remainder theorem for polynomial rings and evaluation at infinity; module homomorphism (recommended).
    - References: [KO62].
    - Additional references: [Too63].
- `Karatsuba-striding.c`: This file demonstrates striding followed by Karatsuba.
    - Assumed knowledge: Chinese remainder theorem for polynomial rings and evaluation at infinity; module homomorphism (recommended).
    - References: [KO62], [Section 3, Ber01].
    - Additional references: [Too63].
- `TC.c`: This file demonstrates Toom-4.
    - Assumed knowledge: Chinese remainder theorem for polynomial rings and evaluation at infinity; module homomorphism (recommended).
    - References: [Too63].
    - Additional references:
- `TC-striding.c`: This file demonstrates striding followed by Toom-4.
    - Assumed knowledge: Chinese remainder theorem for polynomial rings and evaluation at infinity; module homomorphism (recommended).
    - References: [Too63], [Section 3, Ber01].
    - Additional references:
- `Toeplitz-TC.c`: This file demonstrates Toeplitz matrix-vector product built upon Toom-4.
    - Assumed knowledge: Module-theoretic dual of algebra homomorphisms over commutative rings.
    - References: [Too63], [Fid73], [Win80].
    - Additional references:
- `TMVP-polymul`: TBA.
    - Assumed knowledge:
    - References: [Win80].
    - Additional references:
- `Nussbaumer.c`: This file demonstrates Nussbaumer FFT.
    - Assumed knowledge: Chinese remainder theorem for multi-variate polynomial rings.
    - References: [Nus80].
    - Additional references: [SS71], [Sch77], [CK91].
- `Schoenhage.c`: This file demonstrates Schoenhage FFT.
    - Assumed knowledge: Chinese remainder theorem for multi-variate polynomial rings.
    - References: [Sch77].
    - Additional references: [SS71].
- `Rader.c`: TBA.
    - Assumed knowledge: Galois theory.
    - References: [Rad68].
    - Additional references: [Ber22].
- `Bruun.c`: TBA.
    - Assumed knowledge: Chinese remainder theorem for polynomial rings (minimum). Galois theory (recommended).
    - References: [BGM93].
    - Additional references: [Bru78].

# References

[AB74]
Ramesh C. Agarwal and Charles S. Burrus. Fast convolution using Fermat number transforms with applications to digital filtering. IEEE Transactions on Acoustics, Speech, and Signal Processing, 22(2):87–97, 1974. https://ieeexplore.ieee.org/abstract/document/1162555.

[Ber01]
Daniel J. Bernstein. Multidigit multiplication for mathematicians. 2001. https://cr.yp.to/papers.html#m3.

[Ber22]
Daniel J. Bernstein. Fast norm computation in smooth-degree abelian number fields. 2022. https://eprint.iacr.org/2022/980.

[BGM93]
Ian F. Blake, Shuhong Gao, and Ronald C. Mullin. Explicit Factor- ization of x2k + 1 over Fp with Prime p ≡ 3 mod 4. Applicable Alge- bra in Engineering, Communication and Computing, 4(2):89–94, 1993. https://link.springer.com/article/10.1007/BF01386832.

[Bru78]
Georg Bruun. z-transform DFT Filters and FFT’s. IEEE Transactions on Acoustics, Speech, and Signal Processing, 26(1):56–63, 1978. https://ieeexplore.ieee.org/document/1163036.

[CF94]
Richard Crandall and Barry Fagin. Discrete Weighted Trans- forms and Large-integer Arithmetic. Mathematics of computa- tion, 62(205):305–324, 1994. https://www.ams.org/journals/mcom/1994-62-205/S0025-5718-1994-1185244-1/?active=current.

[CK91]
David G. Cantor and Erich Kaltofen. On Fast Multiplication of Polynomials over Arbitrary Algebras. Acta Informatica, 28(7):693–701, 1991. https://link.springer.com/article/10.1007/BF01178683.

[CT65]
James W. Cooley and John W. Tukey. An Algorithm for the Ma- chine Calculation of Complex Fourier Series. Mathematics of Com- putation, 19(90):297–301, 1965. https://www.ams.org/journals/mcom/1965-19-090/S0025-5718-1965-0178586-1/.

[Fid73]
Charles M. Fiduccia. On the Algebraic Complexity of Matrix Multiplication.
1973. https://cr.yp.to/bib/entries.html#1973/fiduccia-matrix.

[Goo58]
I. J. Good. The Interaction Algorithm and Practical Fourier Analysis. Journal of the Royal Statistical Society: Series B (Methodological), 20(2):361– 372, 1958. https://www.jstor.org/stable/2983896.

[GS66]
W. M. Gentleman and G. Sande. Fast Fourier Transforms: For Fun and Profit. In Proceedings of the November 7-10, 1966, Fall Joint Computer Conference, AFIPS ’66 (Fall), pages 563–578. Association for Computing
Machinery, 1966. https://doi.org/10.1145/1464291.1464352.

[KO62]
A. Karatsuba and Yu. Ofman. Multiplication of many-digital numbers by automatic computers. In Doklady Akademii Nauk, volume 145(2), pages 293–294, 1962. http://cr.yp.to/bib/1963/karatsuba.html.

[Nus80]
Henri J. Nussbaumer. Fast Polynomial Transform Algorithms for Digital Convolution. IEEE Transactions on Acoustics, Speech, and Signal Pro- cessing, 28(2):205–215, 1980. https://ieeexplore.ieee.org/document/1163372.

[Rad68]
Charles M. Rader. Discrete Fourier Transforms When the Number of Data Samples Is Prime. Proceedings of the IEEE, 56(6):1107–1108, 1968. https://ieeexplore.ieee.org/document/1448407.

[Sch77]
Arnold Schönhage. Schnelle Multiplikation von Polynomen über Körpern der Charakteristik 2. Acta Informatica, 7(4):395–398, 1977. https://link.springer.com/article/10.1007/bf00289470.

[SS71]
Arnold Schoenhage and Volker Strassen. Schnelle Multiplikation großer Zahlen. Computing, 7(3-4):281–292, 1971. https://link.springer.com/article/10.1007/BF02242355.

[Tho63]
Llewellyn Hilleth Thomas. Using a computer to solve problems in physics. Applications of digital computers, pages 44–45, 1963.

[Too63]
Andrei L. Toom. The Complexity of a Scheme of Functional Elements Realizing the Multiplication of Integers. Soviet Mathematics Doklady, 3:714–716, 1963. http://toomandre.com/my-articles/engmat/MULT-E.PDF.

[Win80]
Shmuel Winograd. Arithmetic Complexity of Computations, volume 33. Society for Industrial and Applied Mathematics, 1980. https://epubs.siam.org/doi/10.1137/1.9781611970364.




