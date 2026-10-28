# Sparse exact protected gates

The [concatenated construction](minimum_readout_safe_gates.md) can be made sparse without changing its exact single-error transparency, interaction order, or native coupling budget. For each logical Pauli axis on the code with $n=5^m$ and $d=3^m$, there are two complementary conclusions:

```math
\boxed{\ell=d+2,\qquad J=(d+1)\omega,\qquad r_S=\infty,}
```

- A general support-reduction argument gives a Hamiltonian with at most $3n+1$ nonzero Pauli terms per axis.
- An explicit recursive algorithm uses at most four translation orbits, hence at most $4n$ terms per axis. It operates on a constant number of Pauli words at each recursion level, uses rational arithmetic of polynomial bit complexity, and needs no exponentially large candidate list.

The second algorithm also bounds the bit lengths of its rational coefficients by $O(\log n)$. A final reduction of its explicit native terms makes the $3n+1$ bound polynomial-time constructive too, with a larger general coefficient-size bound. A separate explicit twenty-five-qubit example uses 75 terms per axis with bare coefficients $2/25$, $7/50$, and $9/50$, improving on the earlier 28,000-term expansion. These are distinct support and coefficient guarantees.

The linear support-reduction principle is standard finite-dimensional linear algebra and linear programming. The result here applies it to signed logical Pauli representatives and then uses the concatenated code's translation symmetries to make sparse synthesis constructive with controlled coefficient size.

## 1. Exact transparency as linear incidence equations

Fix a logical Pauli axis $\alpha\in\{X,Y,Z\}$ and the concatenated encoding $V$. For every physical Pauli label $P_Q$ used below, retain its sign $\eta_Q\in\{+1,-1\}$ so that the signed Hermitian representative

```math
Q=\eta_QP_Q\qquad\text{satisfies}\qquad QV=V\alpha.
```

Every representative lies in the same nontrivial logical stabilizer coset. Define its local-axis incidence vector $v_Q\in\{0,1\}^{3n}$: coordinate $(i,a)$ is one precisely when its local label at site $i$ is $a\in\{X,Y,Z\}$. Global Pauli signs do not enter this vector.

Write a real Hamiltonian in these signed-coset coordinates as $H=\sum_Q c_Q Q$. The conditions

```math
\sum_Qc_Q=1,\qquad
\sum_Qc_Qv_Q=0
\qquad\text{(S1)}
```

are equivalent to the full identities $HV=V\alpha$ and $HEV=EV\alpha$ for all single-site Paulis $E$. Indeed, each individual term satisfies $QEV=\chi(E,Q)EV\alpha$. Conjugation by $X_i$, for example, changes the sign of precisely the $Y$ and $Z$ labels at site $i$. Equation (S1) therefore makes every conjugated coefficient sum one. Conversely, equality of the three conjugated sums forces all three local signed incidence masses to vanish. These equations establish the full physical-space intertwining, including preservation of the protected error images; they are not only compressed logical-energy tests.

The construction uses only weights $d$ and $d+2$, with $c_Q\ge0$ in the lower class and $c_Q\le0$ in the upper class. This sign refers to $c_Q$, not necessarily to the coefficient of the literal bare word $P_Q$. Put $\varepsilon_Q=+1$ or $-1$ according to these two classes and write $c_Q=\varepsilon_Qx_Q$, with $x_Q\ge0$.

## 2. At most $3n+1$ native Pauli terms

Form the real matrix with columns

```math
a_Q=\varepsilon_Q\begin{pmatrix}1\\v_Q\end{pmatrix},
\qquad
Ax=\begin{pmatrix}1\\0\end{pmatrix},
\qquad x\ge0.
\qquad\text{(S2)}
```

The dense two-class Hamiltonian proves feasibility. If a feasible solution has linearly dependent positive-support columns, choose a nonzero dependence $Ah=0$ on that support, orient it to have a positive entry, and subtract

```math
\lambda h,\qquad
\lambda=\min_{h_Q>0}\frac{x_Q}{h_Q}.
```

No coefficient becomes negative, the equations remain true, and at least one positive coefficient vanishes. Repetition leaves independent support columns. Since $A$ has $3n+1$ rows, the final Hamiltonian uses at most $3n+1$ terms. Integer matrix entries and an integer right-hand side give rational coefficients on every such independent support.

The support reduction preserves the native coupling budget automatically. Let $M_d$ and $M_{d+2}$ be the sums of the nonnegative variables in the two classes. The normalization row and the sum of all incidence rows give

```math
M_d-M_{d+2}=1,\qquad
dM_d-(d+2)M_{d+2}=0.
```

Thus every feasible solution has

```math
M_d=\frac{d+2}{2},\qquad
M_{d+2}=\frac d2,\qquad
J=M_d+M_{d+2}=d+1.
\qquad\text{(S3)}
```

Both masses are positive, so both weights remain present and the interaction order is exactly $d+2$. The code and its $r_S=\infty$ property are unchanged. The [locality and coupling optimality proofs](minimum_readout_safe_gates.md#4-full-error-transparency-and-optimal-resources) therefore continue to apply. This argument by itself does not give an efficient way to list the exponentially many columns of the original construction; the next construction removes that requirement.

## 3. Translation orbits reduce the constraint dimension

Label the physical qubits by vectors in $\mathbb F_5^m$. Let

```math
G_m=(C_5)^m
```

act by coordinatewise translations. This group has $n=5^m$ elements and acts regularly: for every pair of sites, exactly one group element takes the first to the second.

These permutations preserve the signed concatenated stabilizer and every logical Pauli coset. For the five-qubit seed, a cyclic shift preserves the positive stabilizer group and fixes $XXXXX$ and $ZZZZZ$. Inductively, a translation of the first coordinate cyclically permutes the five child blocks, while translations of the remaining coordinates act identically within those blocks. Inner signed stabilizers are preserved, and the lifted outer generators transform by the seed's signed stabilizer automorphism. Canonical logical representatives $X^{\otimes n}$, $Y^{\otimes n}$, and $Z^{\otimes n}$ are fixed by permutations. Hence a signed word implementing $\alpha$ still implements the same $\alpha$ after translation.

For a signed seed word $Q$, write its global axis-count vector as

```math
c(Q)=(c_X(Q),c_Y(Q),c_Z(Q)),
\qquad c_X+c_Y+c_Z=\operatorname{wt}(Q),
```

and define its orbit average

```math
\mathcal O(Q)=\frac1n\sum_{g\in G_m}g(Q).
\qquad\text{(S4)}
```

Every term in the average acts as the same logical axis. At any specified physical site, the coefficient mass of local label $a$ is exactly $c_a(Q)/n$: regularity presents each original site exactly once. Consequently a Hamiltonian $\sum_Q a_Q\mathcal O(Q)$ is transparent precisely when

```math
\sum_Qa_Q=1,\qquad
\sum_Qa_Qc(Q)=0.
\qquad\text{(S5)}
```

The $3n+1$ incidence conditions have reduced to four scalar equations. Orbit averages are algebraic sums of concurrent native terms, not physical randomization or a prescription to apply sequential pulses.

## 4. Constructing the orbits without a large candidate list

For each logical axis and each of the two weights $w=d,d+2$, maintain a probability distribution on at most three signed seed words such that

```math
\mathbb E[c(Q)]=\left(\frac w3,\frac w3,\frac w3\right).
\qquad\text{(S6)}
```

Their orbit averages then have exactly the isotropic local-axis counts used in the original two-class proof. The three-seed bound follows because all count vectors lie in the affine plane $c_X+c_Y+c_Z=w$: normalization together with the three count equations has rank at most three. Apply the same dependence-elimination argument as above to reduce any feasible positive mixture to independent support of size at most three. At the seed level the finite ten-word and six-word classes supply the required distributions.

To pass from level $m-1$ to $m$, fix one signed weight-three outer representative of the desired logical axis. It has three active child blocks with local labels $\beta_1,\beta_2,\beta_3$. In the lower class, use the three corresponding child lower-class seed distributions. In the upper class, fix one active child to use its upper-class distribution and let the other two use their lower classes. There is no need to average all ten outer representatives or all three promotion positions: the expected global count of each child is independent of its logical label.

Couple the three child distributions by putting their cumulative probabilities on the same interval $[0,1)$. Refine the interval at the union of their breakpoints, and on each subinterval choose the three corresponding child words. If the support sizes are $k_1,k_2,k_3$, the resulting joint distribution has at most $k_1+k_2+k_3-2\le7$ outcomes. It preserves all three marginals exactly and uses rational interval lengths. Independence of the children is unnecessary.

Tensor the three selected signed child words into their prescribed blocks, with the sign of the outer representative. Every resulting seed implements the desired logical axis. Its weight is $3d_{m-1}=d_m$ in the lower class or $3d_{m-1}+2=d_m+2$ in the upper class. The expected count vector is the sum of the child expectations, giving (S6) in either case. Reduce these at most seven candidate seeds back to at most three.

Maintain these separate class distributions at every level. For the final Hamiltonian, multiply the lower mixture by $(d+2)/2$ and the upper mixture by $-d/2$. The at most six resulting orbit columns

```math
\varepsilon_Q\begin{pmatrix}1\\c(Q)\end{pmatrix}
```

satisfy (S5). A final nonnegative support reduction leaves at most four independent orbit columns. Equations (S3) still fix the two total masses and the native budget. Since every orbit has at most $n$ distinct words, this constructs a Hamiltonian with at most $4n$ terms per axis.

### Coefficient size and computation

The reduction is exact and deterministic once seed order and pivot choices are fixed. The maintained class mixtures solve integer-column systems of rank at most three with right-hand side denominator at most three. Independent-support coefficients are therefore rational numbers whose numerator and denominator bit lengths are $O(\log n)$ by Cramer's rule and the determinant bound for fixed-size integer matrices. Coupling three such distributions and reducing at most seven candidates uses a constant number of rational operations of polynomial bit complexity.

For the final orbit coefficients, select an independent row basis including the normalization row. Its entries in that row are $\pm1$ and its remaining entries have magnitude at most $d+2$. For rank $r\le4$, the Hadamard bound gives a common coefficient denominator $D$ satisfying the safe bound

```math
D\le16(d+2)^3.
```

Expanding each orbit divides by $n$, so every collected bare Pauli coefficient has denominator dividing $nD$, in particular at most $16n(d+2)^3\le16n^4$. Repeated words within an orbit or overlapping equal-weight orbits have consistent signed-coset signs; collecting them preserves this common-denominator bound and the absolute coefficient sum. The normalized coefficients have magnitude at most $J=d+1$, so their numerators also have $O(\log n)$ bits.

There are only six maintained class distributions, each with at most three seed words, at every level. The recursion has depth $\log_5 n$. Constructing these words, performing the fixed-dimensional reductions, and explicitly listing the final at most $4n$ words of length $n$ all have polynomial bit complexity in $n$. This is a constructive statement about the native Pauli description; it does not bound geometric range or experimental control overhead.

If the smaller support bound is desired, apply the full incidence reduction (S2) to this explicit list of at most $4n$ terms. It has only $O(n)$ columns, so the reduction reaches at most $3n+1$ terms in polynomial bit complexity without an exponential listing step. An independent-support solution of the resulting integer system has numerator and denominator bit lengths $O(n\log n)$ by Cramer's rule and the Hadamard bound. Thus the four-orbit output gives the stronger $O(\log n)$ coefficient-size guarantee, while the final native-term reduction gives the smaller general support bound. Neither output is asserted to minimize the number of terms or rational denominators.

## 5. A twenty-five-qubit example with 75 terms

The following three seeds give a simpler explicit certificate at $n=25,d=9$. Spaces separate the five consecutive child blocks. The sign column specifies the signed representative $Q=\eta_QP_Q$ acting as logical $Z$.

| Seed | Literal bare Pauli word $P_Q$ | $\eta_Q$ | Weight | $(c_X,c_Y,c_Z)$ | Coefficient of $\mathcal O(Q)$ |
|---|---|---:|---:|---|---:|
| $A$ | `ZIYIZ IIIII IIIII ZIYIZ ZIXXI` | $+1$ | 9 | $(2,2,5)$ | $2$ |
| $B$ | `YXIIX IIIII IIIII YXIIX YIIYZ` | $+1$ | 9 | $(4,4,1)$ | $7/2$ |
| $C$ | `ZXXZY IIIII IIIII YXIIX YIIYZ` | $-1$ | 11 | $(4,4,3)$ | $-9/2$ |

Each translation orbit has size 25, and the three orbits are disjoint. Set

```math
H_Z=2\mathcal O(A)+\frac72\mathcal O(B)-\frac92\mathcal O(C).
```

Normalization and count balance are immediate:

```math
2+\frac72-\frac92=1,\qquad
2(2,2,5)+\frac72(4,4,1)-\frac92(4,4,3)=(0,0,0).
```

Thus all 76 no-error/single-Pauli intertwining identities hold exactly. Because $C$ is the negative of its bare word, every literal bare coefficient in this particular Hamiltonian is positive: the three orbit coefficients per term are $2/25$, $7/50$, and $9/50$. Their absolute sum is $2+7/2+9/2=10$. This gives 75 terms, interaction order 11, and $J/\omega=10$.

The two weight classes in this example are not individually isotropic. Their signed global count vectors cancel in the displayed combination, which is sufficient by (S5). The recursive class invariant (S6) is a device for producing feasible orbit candidates; it need not survive the final joint reduction.

Uniformly cycling physical labels $X\mapsto Y\mapsto Z\mapsto X$ preserves the signed five-qubit stabilizer and its concatenation. It maps this logical-$Z$ example into logical-$X$ and logical-$Y$ examples with the same resources. The finite verification checks these signs and identities for each axis independently. The recursive four-orbit algorithm need not choose these particular three orbits; its worst-case bound is $4n$, or 100 terms at $n=25$.

## 6. Scope and verification

The construction is exact for one arbitrary single-qubit fault at an unknown time under the [fixed-code gate model](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/model.md). The unchanged code retains exact post-syndrome recoverability for every specified product unitary when recovery knows that unitary. These are separate operational guarantees. The [higher-error obstruction](minimum_readout_safe_gates.md#6-scope-and-verification) still requires $\ell\ge d+2t$ on this family, so the present $\ell=d+2$ gates do not establish protection through weight two when the code can correct that many errors.

The locality and coupling optima apply separately to the three Pauli axes. Real linear combinations preserve transparency and allow arbitrary continuous control of one logical qubit, with budget $(d+1)\sum_\alpha|\Omega_\alpha|$. This does not assert coupling optimality for every arbitrary logical direction, entangling logical gates, or multiqubit computational universality.

Sparse term count does not remove the need for simultaneous cancellations throughout a pulse. A direct calibration bound follows in the coefficient-error model $\Delta H=\sum_Q\delta c_QQ$, using signed positive logical representatives. For a single-site error $E$, the difference between its perturbed logical rate and the perturbed no-fault rate is

```math
\Delta\omega_E-\Delta\omega_I
=-2\sum_{Q:\,QE=-EQ}\delta c_Q,
\qquad
|\Delta\omega_E-\Delta\omega_I|
\le2\sum_Q|\delta c_Q|.
```

If the no-fault logical rate is calibrated to unit strength and the remaining coefficient errors satisfy $|\delta c_Q|\le\eta|c_Q|$, this is at most $2\eta J=2\eta(d+1)$. This algebraic bound applies to both dense and sparse descriptions. The rational bit-length bounds describe exact coefficient specification; they do not supply a hardware calibration guarantee. Geometric locality, synthesis from other interactions, recovery-device faults, and noisy syndrome extraction require separate analysis.

The [verification guide](../verification/README.md) records the exact sparse certificates. They check signed stabilizer membership, orbit sizes, count balance, native coefficients, and the full single-error identities. Recursive instances support the implemented algorithm; its all-depth guarantee and coefficient-size bound follow from the proofs above.
