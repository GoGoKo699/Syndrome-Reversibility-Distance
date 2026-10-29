# Exact sparse synthesis for two-error gates

The [growing-distance construction](growing_two_error_gates.md) proves that the depth-$m$ concatenated five-qubit code admits a two-error-transparent logical Pauli gate with

```math
n=5^m,\qquad d=3^m,\qquad
\ell=d+4,\qquad
\frac{J}{\omega}=\frac{4d^2+21d+8}{8}.
```

This note separates three further statements: positive sparse certificates for each of the four twenty-five-qubit classes; an exact pairwise coupling method sufficient for recursive synthesis; and a finite 125-qubit gate certificate with 69,125 native Pauli terms per axis. The first two are reusable ingredients. The finite certificate has independent exact signed-membership, coefficient and complete error-response checks.

The model is the [fixed-code full error-image transparency model](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/model.md). Interaction order and native absolute coefficient sum are different resources from the number of Pauli terms, coefficient bit length, and the runtime needed to find those terms. All statements here protect physical Pauli errors of weight at most two.

## 1. Positive class certificates at twenty-five qubits

Write $A_0,A_1,A_2,F$ for the four normalized classes in the [all-depth construction](growing_two_error_gates.md#2-four-classes-and-the-frontier-padding-rule). Every atom is an actual Hermitian signed Pauli that acts as the prescribed positive logical axis on the code. Within each class the masses are strictly positive and sum to one. The negative signs in the final Hamiltonian enter only when these normalized classes are combined.

The [stored positive certificates](../verification/data/two_error_positive_classes.json) give the logical-$Z$ classes as orbit averages under the coordinate translations $C_5\times C_5$. Global cyclic relabeling $X\mapsto Y\mapsto Z\mapsto X$ supplies the other axes. The [loader and independent checker](../verification/positive_two_error.py) retain the actual Hermitian coset signs, expand the orbits, and check each signed logical membership.

| Class | Stored orbit seeds | Expanded native terms | Common denominator bits of orbit masses |
|---|---:|---:|---:|
| $A_0$ | 78 | 1,950 | 48 |
| $A_1$ | 96 | 2,400 | 79 |
| $A_2$ | 88 | 2,200 | 63 |
| $F$ | 93 | 2,325 | 126 |

These are certificates for the individual class response vectors. They are not obtained merely by separating an arbitrary sparse final gate into four pieces: the final gate equations alone do not force each piece to have the dense class response vector.

For a Pauli word $Q$ and physical Pauli error $E$, let $\chi(E,Q)\in\{+1,-1\}$ be its conjugation character. The positive certificates preserve

```math
f_T(E)=\mathbb E_{Q\sim T}[\chi(E,Q)]
\qquad (\operatorname{wt}(E)\le2)
```

for every class $T$. A translation-invariant twenty-five-qubit distribution has only

```math
1+3+\frac{25-1}{2}\,3^2=112
```

distinct error-orbit coordinates: identity, three single-site axes, and nine ordered axis pairs for each nonzero displacement modulo reversal. Thus preserving all 112 coordinates preserves all 2,776 physical error responses, including identity. The checked profiles are

| Error pattern | $A_0$ | $A_1$ | $A_2$ | $F$ |
|---|---:|---:|---:|---:|
| Identity | $1$ | $1$ | $1$ | $1$ |
| One nonidentity Pauli | $13/25$ | $31/75$ | $23/75$ | $23/75$ |
| Two sites in one bottom block | $9/25$ | $89/225$ | $97/225$ | $11/75$ |
| Sites in distinct bottom blocks | $29/125$ | $13/125$ | $1/225$ | $23/375$ |

The checker verifies every expanded native word and every physical error response on each of the three logical axes. Translation symmetry compresses the calculation; it does not replace the individual-error guarantee by an average over errors.

## 2. Pairwise coupling is sufficient

Fix one outer signed minimum logical representative and one choice of child classes in a parent mixture. Let $Q_1,\ldots,Q_s$ be the nonidentity child factors, with $s\le4$. The first three factors occupy active children; the fourth can be a localized stabilizer-padding factor in an inactive child. For a physical error $E$ restricted to these children,

```math
\chi(E,Q_1\otimes\cdots\otimes Q_s)
=\prod_{j=1}^s\chi(E_j,Q_j).
```

If $\operatorname{wt}(E)\le2$, at most two child restrictions are nonidentity. Two errors in one child require its marginal response vector through weight two. Errors in two different children require the corresponding product of single-child responses. Consequently, preserving each child marginal and every pairwise product marginal preserves the complete parent response vector through weight two. Full joint independence of the child choices is unnecessary.

Apply the coupling separately after fixing the root representative and child-class mixture. Every output atom is still a tensor product of allowed signed child words. Its logical coset, physical class pattern and interaction order are preserved. The argument uses only response vectors and signed logical membership; it does not require sparse child words to carry an auxiliary sample of the original active support tree.

Preserving only the one-child marginals is insufficient. For example, let the character of a fixed error in each of two children be uniformly $+1$ or $-1$. Independent child choices give a two-child response of zero. Choosing the two character signs to be equal preserves both marginals but gives a two-child response of one. This is the obstruction to reusing a coupling method justified only by single-error profiles.

### An exact rational construction

The standard-library helper [pairwise_coupling.py](../verification/pairwise_coupling.py) takes one to four finite positive rational distributions and returns a finite joint distribution with all the required product pair marginals. For three or four inputs, let $U,V$ be uniform on the unit torus and use the quantile functions of the inputs at

```math
U,\qquad V,\qquad U+V\pmod1,\qquad U+2V\pmod1.
```

Every pair of these linear forms is jointly uniform. Its integer coefficient matrix has nonzero determinant, so the induced torus homomorphism is onto and sends Haar measure to Haar measure. Determinant two for the pair $(U,U+2V)$ is allowed; invertibility as a one-to-one map is not required. Applying the separate quantile maps gives the exact product marginal for every pair of input distributions, including unequal rational masses.

Quantile breakpoints cut the unit square along rational lines of the forms

```math
U=a,\qquad V=b,\qquad U+V=c+z,\qquad U+2V=e+z,
```

where $z$ ranges over the finitely many relevant integer periods. Each resulting polygon has a constant tuple of quantile values. Its exact rational area is the mass assigned to that tuple. The helper clips these polygons with rational arithmetic and merges pieces that give the same tuple; it does not sample a grid whose size is a coefficient denominator.

If $K$ is the sum of the input support sizes, there are $O(K)$ lines and $O(K^2)$ cells because the number of factors is at most four. Hence the output has $O(K^2)$ atoms. Rational line intersections and polygon areas have bit length polynomial in the input description. One input is returned directly; for two inputs the helper returns their exact product distribution.

## 3. Recursive synthesis and its computational scope

Let a child code have $k'=3^{m-2}$ active bottom blocks and let the parent have $k=3k'$. The exact class recurrences are

| Parent class | Child-class mixture |
|---|---|
| $A_0$ | Three $A_0$ children. |
| $A_1$ | Choose one active child uniformly for $A_1$; the other two use $A_0$. |
| $A_2$ | With probability $(k'-1)/(3k'-1)$, choose one active child uniformly for $A_2$. With probability $2k'/(3k'-1)$, choose two active children uniformly for $A_1$. All remaining active children use $A_0$. |
| $F$ | With probability $3(k'-1)/(k-1)$, choose one active child uniformly for $F$. With probability $2/(k-1)$, choose one inactive root child uniformly and place $K_4$ in a uniformly chosen bottom block inside it. All remaining active children use $A_0$. |

Each row is also averaged over the ten signed outer minimum representatives of the desired logical axis. The two $F$ alternatives count the descendant inactive frontiers and the two frontiers created at the root, respectively. A localized $K_4$ factor has an explicit positive distribution on $15\cdot5^{m-2}=O(n)$ words. The [frontier identity](growing_two_error_gates.md#3-the-pair-correlation-identity) proves that this recurrence has the required dense class profile.

At each depth, maintain a separate positive normalized sparse distribution for every class and logical axis. Couple its child factors by the preceding procedure. The reference compiler [growing_two_error_sparse.py](../verification/growing_two_error_sparse.py), through `synthesize_classes(level)`, averages each candidate over coordinate translations and streams the resulting integer orbit-response column into an exact conic basis reduction. The child orbit distributions are expanded to native words before coupling; coupling orbit seeds alone would omit the required independent or pairwise coordinate translations.

The full response vector is

```math
\left(\mathbb E[\chi(E,Q)]\right)_{\operatorname{wt}(E)\le2}
\in\mathbb Q^{D_2(n)},\qquad
D_2(n)=1+3n+9\binom n2.
```

The reference compiler preserves its translation-orbit compression, of dimension

```math
C_2(n)=4+\frac{9(n-1)}2=O(n).
```

The identity coordinate preserves normalization. Exact elimination of dependent orbit columns leaves at most $C_2(n)$ positive orbit seeds per class. Expansion therefore gives at most

```math
nC_2(n)=D_2(n)+n-1=O(n^2)
```

native words per class. This is the executable compiler's bound. Direct reduction of unrestricted native response columns gives the slightly smaller theoretical bound $D_2(n)$; the compiler does not implement that additional final reduction. The checked twenty-five-qubit classes provide the recursion base.

If the child lists have $O(n^2)$ support, the pairwise coupling has $O(n^4)$ candidate atoms before support reduction. The finite number of root and class mixtures changes only the constant. This improves the loose $O(n^8)$ candidate bound obtained by expanding four fully independent child lists. Both bounds are polynomial, and neither is a practical runtime estimate for a large concatenation depth.

The class target profiles have polynomially bounded rational bit length, as shown in the [all-depth description-size argument](growing_two_error_gates.md#6-description-size-and-scope). Orbit-response columns have integer entries bounded by $n$. Independent-support coefficients have polynomial bit length by Cramer's rule and the Hadamard bound. The same tools give the conservative bound $O(D_2(n)\log n)$ for unrestricted native reduction. Together with the polynomial candidate bound, this gives an exact algorithm with polynomial bit complexity and a high polynomial degree. The streaming implementation uses rational arithmetic throughout.

Finally combine the four positive classes with the forced masses

```math
\left(\frac{(d+2)(d+4)}8,
-\frac{d(d+4)}4,
\frac{d(2d+9)}{16},
-\frac{5d}{16}\right).
```

The implemented `synthesize_hamiltonian(level, axis="Z")` mixes these class distributions and reduces their signed orbit-response columns. It returns at most $C_2(n)$ signed orbits, hence at most $D_2(n)+n-1$ native terms in the final gate. It checks the four forced class masses exactly. The block-moment equations ensure that the final reduction preserves the optimum $J/\omega=(4d^2+21d+8)/8$.

A further unrestricted reduction of native response columns can leave at most $D_2(n)$ native terms while preserving those masses and the same optimum. This slightly stronger support bound remains theoretical; that final native reduction is not implemented.

The positive base loader, rational pairwise coupling and general reference compiler are executable. An end-to-end depth-two reference compilation returns 103 signed orbits, or 2,575 native terms. The full depth-three reference compilation has not been used as a practical synthesis procedure; the finite certificate below uses a numerical support search over recursive candidates instead. Polynomial complexity does not establish fast large-depth synthesis, small coefficient denominators, minimal term count, or an efficient physical implementation of the resulting many-body terms.

## 4. The 125-qubit finite certificate

At depth three the native targets are

```math
n=125,\qquad d=27,\qquad \ell=31,\qquad
\frac J\omega=\frac{3491}{8}.
```

The four signed masses are $899/8,-837/4,1701/16,-135/16$ in the order $A_0,A_1,A_2,F$. The [all-depth moment proof](growing_two_error_gates.md#5-optimality-for-arbitrary-hamiltonians-on-this-family) establishes the interaction-order and coupling lower bounds independently of any numerical search.

Index physical sites by three base-five coordinates. The coordinatewise translation group $C_5^3$ has 125 elements and preserves the signed code. An invariant gate has

```math
1+3+\frac{125-1}{2}\,3^2=562
```

error-orbit response coordinates: identity, three single-error axes, and nine ordered axis pairs for each of 62 nonzero displacement pairs $\{\delta,-\delta\}$. Every physical Pauli of weight at most two belongs to one of these orbits. Their total number is

```math
1+3(125)+9\binom{125}{2}=70,126.
```

A numerical linear program selected support from recursive signed logical-coset candidates using the four fixed class signs and these 562 response equations. Its floating-point coefficients are not used in the certificate. Exact rational reconstruction by fraction-free elimination produced the [stored orbit masses](../verification/data/two_error_125_orbits.json). The [independent standard-library verifier](../verification/verify_two_error_synthesis.py) checks positivity, signed logical membership, interaction order and every response equation with exact arithmetic. No numerical optimizer is required to load or verify this gate.

All terms have weight 27, 29 or 31. Every nonidentity coordinate translation has order five; an invariant Pauli label would have weight divisible by five. These weights therefore force each native term orbit to have 125 distinct words. After checking that the stored orbit seeds represent disjoint orbits, $N$ seeds give exactly $125N$ native terms. An independent-support orbit certificate uses at most 562 seeds, hence at most 70,250 native terms. This orbit-symmetric bound is separate from the unrestricted $D_2(125)=70,126$ support bound.

The exact finite certificate has the following verified resources:

| Certificate item | Verified value |
|---|---|
| Stored data and verifier | [Orbit certificate](../verification/data/two_error_125_orbits.json), [exact verifier](../verification/verify_two_error_synthesis.py) |
| Logical axes checked | $X,Y,Z$ |
| Disjoint orbit seeds | $553$ |
| Expanded native terms per axis | $69,125$ |
| Common denominator bits of native coefficients | $1,458$ |
| Protected error responses | $562$ orbit coordinates covering all $70,126$ physical errors |
| Exact interaction order and coupling | $31$ and $3491/8$ |

For each axis, the verifier expands all 69,125 native words and checks their actual signed logical-coset membership against the 124-generator code. It independently classifies all 70,126 protected physical errors into the 562 coordinate orbits and evaluates their exact integer character sums. Literal site-by-site translated characters independently check a representative column from each physical class. Thus the scalar response calculation certifies the full error-image identity $H_\alpha EV=EV\alpha$, including error-image preservation.

The orbit masses have a 1,451-bit common denominator; division by the orbit size produces a 1,458-bit common denominator for the native coefficients. This is an exact sparse witness at the optimal interaction order and coupling budget, with substantial coefficient complexity. Neither its term count nor its denominator is claimed minimal. The [loader](../verification/growing_two_error_sparse.py), through `load_hamiltonian(3, axis)`, returns each collected native label with its coset sign and relative coefficient; their product is the bare Hermitian Pauli coefficient.

The same verifier checks the positive depth-two class profiles, pairwise coupling on unequal distributions including a 127-bit denominator, tensor responses with two errors in one factor, recursive class identities through depth six, and the end-to-end depth-two reference compiler. These finite checks support the supplied programs; the all-depth conclusion follows from the recursive proof.

The code retains $r_S=\infty$ for each specified known product unitary followed by full ideal original-syndrome measurement, as established for the [underlying family](minimum_readout_safe_gates.md). The gate guarantee remains fixed at two physical Pauli errors. The finite 125-qubit certificate and the general polynomial reference algorithm do not establish fast synthesis at substantially greater depths or a joint noisy gate/readout threshold.
