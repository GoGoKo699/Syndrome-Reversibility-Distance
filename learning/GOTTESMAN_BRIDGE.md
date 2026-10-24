# From Gottesman's stabilizer tutorial to syndrome reversibility

An error syndrome is useful because it can reveal what happened to a physical system while preserving the quantum information encoded in it. This project asks exactly when that separation survives coherent rotations on every physical qubit. The answer depends on how stabilizer supports and logical-operator supports fit together.

This bridge follows the [30-milestone map](MILESTONES.md). It develops the connections needed to read the complete proofs under [science/](../science/model.md), with a small measurement calculation before the general theorems. The scientific pages remain the authoritative statements of assumptions, quantifiers, and proofs.

## 1. The single educational anchor

**Daniel Gottesman, *An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation*.** We use [arXiv:0904.2557v1](https://arxiv.org/abs/0904.2557v1), dated 16 April 2009, published as a chapter in *Quantum Information Science and Its Contributions to Mathematics*, Proceedings of Symposia in Applied Mathematics **68**, 13–58 (2010). The [author's text](https://arxiv.org/pdf/0904.2557) supplies the stabilizer foundation.

The source route is selective:

| Source section | Relevant foundation | Where it enters this project |
|---|---|---|
| §2.1 | Pauli expansions and coherent errors | Physical rotations and their error paths |
| §2.2 | Error-correction conditions | Preserving arbitrary encoded states |
| §2.3 | Stabilizers, syndrome measurement, distance | The fixed code and readout |
| §3.1, especially Table 4 | Five-qubit code and binary representation | Explicit families and syndrome columns |
| §3.4 | Logical Pauli quotient | Distinguishing detectable, trivial, and logical operators |
| §3.5 | Clifford conjugation | Local axis changes and the physical CNOT |
| §5.2, concatenation discussion | Repeated encoding | The growing-distance construction |

The general fault-tolerance proof in §5 is not needed for the present result. Our measurement is ideal, and our small-angle bound is for a fixed code.

Three additions are supplied here rather than assigned to the anchor: the precise optimization over recovery after a recorded measurement; the fidelity and diamond-distance estimates needed for that optimization; and the support/circuit theorems proved in this repository. The first addition uses established results of Gregoratti and Werner. The even-repetition mechanism used later is due to earlier coherent-error work, including Huang, Doherty and Flammia. These primary research sources are attributed in [Related work](../references/RELATED_WORK.md); they do not constitute additional educational anchors.

## 2. Route through the milestones

The map has 12 background/preliminaries milestones and 18 narrative/innovation milestones. The first phase combines individual physical and mathematical objects into the complete operational question. The second combines the proof mechanisms into the code-family separation and its interpretation.

| Milestones | Connection developed here |
|---|---|
| M01–M06 | Encoding, Pauli labels, checks, logical operators, and coherent rotations: §§3–5 |
| M07–M12 | Recorded measurements, recovery, and performance measures: §§4–6 |
| M13–M17 | From measured checks to the exact support theorem: §7 |
| M18–M20 | From a support collision to an even logical circuit: §8 |
| M21–M22 | A uniform bound that includes rare outcomes: §9 |
| M23–M25 | The exact repetition witness and matching exponent: §10 |
| M26–M29 | Five-qubit blocks, concatenation, and one CNOT: §11 |
| M30 | The complete claim, its attribution, and its limits: §12 |

## 3. Translate the notation before comparing statements

A code stores $k$ logical qubits inside $n$ physical qubits. Its logical Hilbert-space dimension is $D=2^k$. An encoding isometry $V$ maps logical states into the physical codespace; $V^\dagger V=I_D$ and $P=VV^\dagger$ is the projector onto that subspace.

| Concept | Convention in this repository |
|---|---|
| Stabilizer | The actual physical group consists of signed Hermitian Paulis fixing the code. For support counting, $S$ denotes its labels modulo global phases. |
| Centralizer | $N$ is the phase-free set of Paulis commuting with every stabilizer, corresponding to the anchor's $N(S)$. |
| Logical Pauli | A label in $N\setminus S$ acts nontrivially on the code. A bar denotes its induced logical operator. |
| Support and weight | Support is the set of nonidentity sites; weight is its size. Equal weight does not imply equal support. |
| Syndrome | $s$ is a binary string recording the eigenvalues of the independent checks. The actual projectors retain the physical signs. |
| Rotation scale | A local rotation is $\exp(-i\theta h_j/2)$. For a Pauli axis, its error amplitude is $\sin(\theta/2)$. |
| Recovery information | Each specified physical unitary may have its own recovery, applied after the prescribed measurement. |
| Error convention | Diamond error means one half of the diamond norm of the difference from the identity channel. |

Discarding phases is legitimate for support and binary commutation calculations because those calculations retain the symplectic commutation relation separately. It is not legitimate when forming the physical projectors. For example, $+XX$ and $-XX$ specify different positive-eigenvalue codespaces even though their phase-free support labels agree. The [model](../science/model.md#1-code-phases-and-syndrome) keeps these two roles separate.

The measurement order is fixed throughout: encode, apply $U=\bigotimes_jU_j$, measure all original checks, then recover using the complete classical record and surviving quantum state. In symbols, the recovery quantifier is

```math
\text{for every specified }U,\quad
\text{there exists a recovery }\mathcal R_U.
```

The quantifier does not supply one universal decoder for an unknown rotation. It also places recovery after readout. Applying a known $U^\dagger$ before readout would answer a different, immediate reversal problem.

## 4. The syndrome induces a measurement on the logical state

For each syndrome $s$, let $\Pi_s$ project onto its physical eigenspace. Choose an isometry $W_s$ identifying that eigenspace with a copy of the logical space. A Pauli return from the syndrome sector provides one such identification. The conditional map and its positive effect are

```math
K_s=W_s^\dagger\Pi_sUV,
\qquad
M_s=K_s^\dagger K_s
=V^\dagger U^\dagger\Pi_sUV.
```

On logical density matrix $\rho$, the outcome has probability

```math
\Pr(s\mid\rho)=\mathrm{Tr}(M_s\rho).
```

This turns a physical syndrome question into a logical measurement question. If $M_s=p_sI_D$, the outcome probability is independent of the encoded state. The corresponding $K_s$ has equal singular values, so it is a scalar times a unitary identification and can be reversed conditionally. If $M_s$ is nonscalar, different logical states produce different outcome probabilities, and the branch changes their relative amplitudes.

The complete readout has one Kraus operator per recorded outcome. This special structure is why input-independent outcome probabilities suffice for exact recovery here. An instrument with additional unobserved errors within an outcome can lose information even when that recorded probability is input independent.

Retaining the syndrome as classical data does not retain coherence between its different values. Conditional recovery can use the recorded value; it cannot recreate arbitrary phases already removed by the measurement. The next example makes that loss visible in two qubits.

## 5. A complete two-qubit example

### 5.1 Code and rotation

Take the positive eigenspace of the stabilizer $XX$. It encodes one qubit with basis

```math
\lvert0_L\rangle=\lvert++\rangle,
\qquad
\lvert1_L\rangle=\lvert--\rangle,
\qquad
\lvert\pm\rangle=
\frac{\lvert0\rangle\pm\lvert1\rangle}{\sqrt{2}}.
```

The physical operator $E=ZZ$ swaps these basis states, so its induced logical action is $\overline X$. A physical $XI$ acts as a nontrivial logical phase operator, so this toy code has ordinary distance $d=1$. It illustrates the measurement mechanism rather than supplying protection against every single-qubit error.

Rotate both physical qubits about $Z$ by the same angle. With $c=\cos(\theta/2)$ and $s=\sin(\theta/2)$,

```math
U_\theta=(cI-isZ)\otimes(cI-isZ)
=c^2II-ics(ZI+IZ)-s^2ZZ.
```

Measure $XX$. The terms $II$ and $ZZ$ commute with the check and give the positive syndrome. The terms $ZI$ and $IZ$ anticommute with it and give the negative syndrome. Thus the readout distinguishes two pairs of coherent paths, without distinguishing the two members of either pair.

### 5.2 The two logical branches

Use the identity return for the positive outcome and $ZI$ for the negative outcome. The logical branch operators are

```math
K_+=c^2I-s^2\overline X,
\qquad
K_-=-ics(I+\overline X).
```

Their effects are

```math
\begin{aligned}
M_+&=(c^4+s^4)I-2c^2s^2\overline X,\\
M_-&=2c^2s^2(I+\overline X),
\end{aligned}
\qquad M_++M_-=I.
```

Let $x=\mathrm{Tr}(\overline X\rho)$ denote the encoded $X$ expectation and let $q=s^2$. Then

```math
\Pr(-\mid\rho)=2q(1-q)(1+x).
```

An encoded $+1$ eigenstate of $\overline X$ gives the negative outcome with probability $4q(1-q)=\sin^2\theta$. An encoded $-1$ eigenstate never gives it. The syndrome is therefore learning a logical observable whenever $cs\ne0$. At $\theta=\pi/2$, it distinguishes those two eigenspaces perfectly.

### 5.3 The residual error and the invariant

For $0\leq\theta\leq\pi/2$, the two-branch fidelity optimum described in §6 is $F_{\mathrm{opt}}=c^2=1-q$. The above Pauli returns already give the averaged channel

```math
\Lambda_\theta(\rho)
=(1-q)\rho+q\overline X\rho\overline X.
```

To see the cancellation, expand both $K_+\rho K_+^\dagger$ and $K_-\rho K_-^\dagger$. Their cross terms have equal size and opposite signs. The identity coefficient is $c^4+c^2s^2=c^2$, and the $\overline X$ coefficient is $s^4+c^2s^2=s^2$. This is dephasing in the eigenbasis of $\overline X$.

Its half-diamond error is $q$, and no recovery can do better because $1-F_{\mathrm{opt}}=q$ is a lower bound on that error. Thus

```math
\varepsilon_{\mathrm{opt}}(U_\theta)
=e_{\mathrm{opt}}(U_\theta)
=q=\frac{\theta^2}{4}+O(\theta^4).
```

Here $g=XX$ and $L=YY$ have equal support; their phase-free product is $ZZ$. Each proper nonempty subword, $ZI$ or $IZ$, has nonzero syndrome. Consequently this is an even logical circuit of weight two, giving $r=2$ although $d=1$. The small example already contains the distinction between ordinary distance and the first order of irreversible loss.

## 6. Recovery and error measures supplied by the project

### Polar decomposition separates rotation from distortion

Write each branch as $K_s=Q_sA_s$, where $A_s=\sqrt{M_s}$ is positive. The unitary factor $Q_s$ rotates the branch output; the positive factor changes amplitudes along its eigenvectors. Reversing $Q_s$ leaves a channel with positive Kraus operators:

```math
\Phi(\rho)=\sum_s A_s\rho A_s,
\qquad \sum_sA_s^2=I_D.
```

A singular branch causes no obstruction to defining this recovery: its polar partial isometry extends to a unitary on the equal-dimensional syndrome and logical spaces. The extension on the branch kernel has no effect on any state actually reaching that branch.

Entanglement fidelity tests preservation of the maximally entangled logical/reference state

```math
\lvert\Omega_D\rangle
=\frac{1}{\sqrt{D}}\sum_{j=1}^D
\lvert j\rangle_L\lvert j\rangle_R.
```

Preserving this state checks the channel's action on superpositions as well as basis states. The optimal fidelity over all outcome-dependent recovery channels is

```math
F_{\mathrm{opt}}(U)
=\frac{1}{D^2}\sum_s
\left(\mathrm{Tr}\sqrt{M_s}\right)^2.
```

This is the Gregoratti–Werner recovery theorem, with a self-contained derivation in [the model](../science/model.md#4-exact-optimal-entanglement-fidelity). Polar recovery attains it. Equality $F_{\mathrm{opt}}=1$ holds exactly when every $M_s$ is scalar, including the zero operator for an impossible outcome.

### Diamond error checks the worst input, including a reference

For a channel $\Lambda$ on the logical system, its half-diamond distance from the ideal channel can be written

```math
\frac{1}{2}\lVert\Lambda-\mathrm{id}\rVert_\diamond
=\sup_{\rho_{LR}}
\frac{1}{2}\left\lVert
(\Lambda\otimes\mathrm{id}_R)(\rho_{LR})-\rho_{LR}
\right\rVert_1.
```

The reference may have dimension $D$. The trace norm $\lVert\cdot\rVert_1$ is the sum of singular values; half of the trace norm of a state difference measures their optimal distinguishability. The diamond distance takes the largest such difference over all logical/reference input states. This is stronger than checking a few logical basis states.

Entanglement infidelity $e_{\mathrm{opt}}=1-F_{\mathrm{opt}}$ is no larger than the optimized half-diamond error $\varepsilon_{\mathrm{opt}}$: the maximally entangled input and the test onto its ideal output are among the possibilities allowed by the latter. The two optima need not agree for a general rotation. They agree for the repetition witnesses because a recovery realizes a Pauli dephasing channel whose half-diamond error saturates the fidelity lower bound.

The source's general error-correction discussion motivates preserving arbitrary encoded inputs. The exact record-assisted optimizer and the reference-sensitive estimates above are additional ingredients supplied by [model.md](../science/model.md), not statements attributed to the selected tutorial.

## 7. Why exact physical support classifies reversibility

The support theorem combines three observations. Each can be understood separately before reading their complete proof.

First, the full syndrome projectors and the full stabilizer group determine one another by binary Fourier inversion:

```math
\Pi_s=\frac{1}{|S|}\sum_{g\in S}\chi_s(g)g,
\qquad
g=\sum_s\chi_s(g)\Pi_s.
```

The signs $\chi_s(g)$ are the eigenvalues of $g$ in syndrome sector $s$. Consequently, every logical effect is scalar precisely when every compressed conjugated stabilizer $V^\dagger U^\dagger gUV$ is scalar. Looking only at the independent generator list would discard information needed for this equivalence.

Second, a one-qubit unitary maps a nonidentity Pauli into a linear combination of $X,Y,Z$ while fixing the identity. Product conjugation can therefore change axes on a stabilizer's active sites but cannot introduce or remove an active site. As the rotations vary, the span contains exactly all Pauli words with that same support.

Third, compression of a Pauli into the code has three possible outcomes:

| Physical Pauli label | Logical compression | Reason |
|---|---|---|
| Outside $N$ | Zero | It changes at least one check eigenvalue and leaves the code sector. |
| In $S$ | A scalar identity | The physical signed stabilizer fixes every encoded state. |
| In $N\setminus S$ | A nontrivial logical Pauli | It preserves the sector while changing encoded information. |

The only obstruction is thus a nontrivial logical Pauli lying in the support orbit of a stabilizer. If no exact supports match, all branches stay reversible for every product unitary. If supports match, local Clifford changes of axis can turn that stabilizer into the logical Pauli; syndrome readout then measures its eigenvalue. The [complete theorem](../science/support_criterion.md#1-the-support-criterion) proves necessity and sufficiency.

This also explains why total weight is too coarse. Two operators of equal weight can act on different physical qubits. Conversely, checking only global even/odd weight can miss a matching support structure in codes with multiple logical qubits.

## 8. From a support collision to a logical circuit

For a matching-support pair $g\in S$ and $L\in N\setminus S$, form the phase-free product $E=gL$. On a site where $g$ and $L$ agree, $E$ is identity. On a site where they differ, both are nonidentity and $E$ anticommutes with the local factor of $g$. Thus the weight of $E$ counts exactly how many local Pauli axes must change in this collision.

Let $\mathcal C_S$ be the set of all such matching-support pairs. The invariant minimizes that count:

```math
r(S)=\min_{(g,L)\in\mathcal C_S}
\mathrm{wt}(gL),
\qquad \min\varnothing=\infty.
```

Every eligible $E$ is a logical Pauli, so $r\ge d$. It anticommutes locally with $g$ at all its active sites but commutes with $g$ globally. The number of those local anticommutations must therefore be even.

To understand the circuit description, assign to each active local factor $E_j$ its binary syndrome column $v_j$. A product over a subset of factors has zero syndrome exactly when its columns sum to zero. A logical circuit is a nontrivial logical Pauli for which the full collection is the only nonempty dependent subset. Here “circuit” means a minimal binary dependence, not a sequence of quantum gates.

A minimum witness has this property. A proper zero-syndrome subword would either be a smaller logical witness, or a stabilizer that could be removed to leave one. Conversely, for an even circuit, the only nonzero vector in the syndrome matrix's kernel is the all-ones vector. Its even length makes that vector orthogonal to itself. Binary row-space/kernel duality then provides a product of stabilizer checks that anticommutes locally at every active site. This supplies the required witness.

The invariant is consequently the minimum even logical-circuit weight. A minimum-distance logical is always a circuit, so even ordinary distance gives $r=d$. For odd distance, the relevant even circuit can be heavier or absent. The exact argument and the meaning of a fixed-label subword are in [the circuit characterization](../science/support_criterion.md#4-equivalent-characterization-by-even-logical-circuits).

## 9. Why rare syndrome outcomes do not weaken the upper exponent

The weak-rotation family is

```math
U_{\boldsymbol h}(\theta)
=\bigotimes_{j=1}^n e^{-i\theta h_j/2},
\qquad h_j=h_j^\dagger,
\qquad \lVert h_j\rVert_\infty\le1.
```

The norm bound fixes the meaning of the shared angle scale. Changing a local Pauli axis in the conjugation expansion costs at least one power of $|\theta|$. Exposing a logical Pauli through a stabilizer therefore costs at least $r$ powers. The collision polynomial keeps track of all such possibilities:

```math
B_S(x)=\sum_{(g,L)\in\mathcal C_S}
x^{\mathrm{wt}(gL)}.
```

Write each effect as $M_s=p_sI_D+\Delta_s$, where $p_s=\mathrm{Tr}(M_s)/D$. The Fourier factor and the Pauli coefficient bound yield

```math
\sum_s\lVert\Delta_s\rVert_\infty
\le B_S(|\theta|).
```

The operator norm $\lVert\cdot\rVert_\infty$ is the largest singular value. This estimate controls the unnormalized branch effects. Dividing by $p_s$ would be dangerous for a rare outcome and is unnecessary.

The positive-Kraus recovery channel has the exact identity

```math
\Phi-\mathrm{id}
=-\frac{1}{2}\sum_s[A_s,[A_s,\,\cdot\,]],
\qquad A_s=\sqrt{M_s}.
```

Subtracting the scalar $\sqrt{p_s}I_D$ leaves both commutators unchanged. The double commutator makes the error quadratic in that difference. At the same time, scalar spectral calculus gives

```math
\left\lVert\sqrt{M_s}-\sqrt{p_s}I_D\right\rVert_\infty^2
\le\lVert M_s-p_sI_D\rVert_\infty.
```

The square root and the quadratic channel estimate balance. As a result,

```math
\varepsilon_{\mathrm{opt}}(U_{\boldsymbol h}(\theta))
\le\min\{1,B_S(|\theta|)\}
=O(|\theta|^r)
```

for a fixed finite-$r$ code and $|\theta|\le1$. A zero-probability effect is the zero matrix by positivity, so it is covered as well. This mechanism is developed in [the upper-bound proof](../science/error_exponent.md#2-a-finite-code-collision-polynomial). The estimate gives one explicit code-dependent constant; it does not bound that constant uniformly over larger codes.

## 10. Why the exponent is attained against every recovery

Choose a minimum even circuit $E$ of weight $r$ and rotate about its local factors. The unitary expands into Pauli subwords $E_A$, indexed by subsets of its support. Two subwords have the same syndrome exactly when their symmetric difference has zero syndrome. Circuit minimality leaves only identical subsets and complementary subsets.

Every reached syndrome therefore contains one pair $A,A^c$. With $c=|\cos(\theta/2)|$, $s=|\sin(\theta/2)|$, and $|A|=w$, a Pauli return gives a logical operator with two amplitudes of magnitudes

```math
a=c^{r-w}s^w,
\qquad b=c^ws^{r-w}.
```

Because $r$ is even, their relative phase is real. The branch has the form $aI\pm b\overline E$, up to an overall phase. Its two singular values are $|a+b|$ and $|a-b|$, each repeated $D/2$ times. The exact recovery formula makes its optimal fidelity contribution $\max(a^2,b^2)$. Thus even an arbitrary quantum recovery loses the smaller probability in each pair.

For $q=s^2\le1/2$, this gives the even-repetition tail

```math
p_r(q)=\sum_{w>r/2}\binom{r}{w}q^w(1-q)^{r-w}
+\frac{1}{2}\binom{r}{r/2}q^{r/2}(1-q)^{r/2}.
```

The factor one half counts equal-weight complementary pairs once. Returning the lower-weight subword realizes a Pauli dephasing channel with rate $p_r(q)$, so its half-diamond error equals the all-recovery fidelity lower bound. Both optimized errors are exactly $p_r(q)$ for this witness.

The tied pairs lead at small angle:

```math
p_r\!\left(\sin^2(\theta/2)\right)
=\frac{1}{2}\binom{r}{r/2}\left(\frac{\theta}{2}\right)^r
+O(\theta^{r+2}).
```

This positive coefficient and the preceding upper bound establish the sharp worst-case exponent for each fixed finite-$r$ code:

```math
\sup_{\boldsymbol h}e_{\mathrm{opt}}(U_{\boldsymbol h}(\theta))
=\Theta(|\theta|^r),
\qquad
\sup_{\boldsymbol h}\varepsilon_{\mathrm{opt}}(U_{\boldsymbol h}(\theta))
=\Theta(|\theta|^r).
```

For $r=\infty$, both errors vanish exactly for every product unitary. At $r=2$ the toy example gives $p_2(q)=q$; at $r=4$ it gives $p_4(q)=3q^2-2q^3$. The repetition mechanism is inherited; the code-wide statement that every finite-$r$ code contains this witness is the new connection. The [complete lower bound](../science/error_exponent.md#3-an-embedded-even-repetition-witness) distinguishes this exact witness equality from the two generally different recovery metrics.

## 11. From the five-qubit code to the one-CNOT separation

The source's Table 4 supplies the standard five-qubit stabilizer. The [family proof](../science/code_families.md#2-the-five-qubit-building-block) checks its distance and the required parity properties directly. Its stabilizers all have even weight, while every representative of a nontrivial logical Pauli has odd weight. These are properties of complete stabilizer and logical cosets, not only their displayed generators.

Concatenation replaces each outer physical qubit by an inner encoded qubit. A nontrivial outer logical must act on at least three inner blocks; each such action costs the inner distance. Minimum representatives attain that product, yielding $[[5^t,1,3^t]]$ for every integer $t\ge1$. The parity property also persists: an even number of odd-weight inner logical factors produces an even-weight outer stabilizer, and odd logical representatives remain odd. These familiar coding ingredients support the new comparison.

Now take two such blocks. A matching stabilizer/logical support would need to match separately inside each block. On at least one block the logical is nontrivial and has odd weight, whereas the stabilizer restriction has even weight. This is impossible, giving $r=\infty$. Total parity alone would be insufficient: a logical acting nontrivially on both blocks can have even total weight.

Conjugate the code by CNOT from physical qubit $0$ in the first block to physical qubit $5^t$ in the second. The construction chooses a minimum logical in the first block that avoids the control, so one weight-$3^t$ logical survives unchanged. A hypothetical shorter logical after CNOT would inverse-conjugate to a weight-$3^t$ logical of the product code touching both blocks. A minimum-weight product-code logical cannot do that, because its nontrivial logical block already uses the entire weight budget. The ordinary distance stays $3^t$.

Meanwhile, a minimum logical $E_B$ with a $Z$ factor at the target transforms into $Z_0E_B$. Its weight is $3^t+1$. A normalizing subword must select the new control factor exactly when it selects the target factor; otherwise inverse CNOT leaves a detectable single-qubit $Z_0$ in the first block. Minimality of $E_B$ then permits only the empty and full subwords. The transformed operator is an even logical circuit. Since $r$ is even and at least the odd distance $3^t$, this circuit forces $r=3^t+1$.

The complete result is

| Physical code | Common parameters | Syndrome reversibility distance |
|---|---|---|
| Two concatenated five-qubit blocks | $[[2\cdot5^t,2,3^t]]$ | $\infty$ |
| Their conjugate under one physical CNOT | $[[2\cdot5^t,2,3^t]]$ | $3^t+1$ |

For $t=1$, the explicit ten-qubit rotation gives the exact error $3q^2-2q^3$. The [verification suite](../verification/README.md) constructs its physical state and all 256 signed syndrome sectors, then checks the recovered channel. The proof supplies the statement for every concatenation level.

The CNOT changes the code and the observables measured afterward. It is not an additional error in the product-unitary noise model. An entangling Clifford can change support relative to physical qubits, so it can change this invariant even when the conventional parameters are preserved.

## 12. Assemble the complete statement

The conceptual chain is now explicit: the prescribed readout induces logical effects; their scalarity decides exact reversibility; product conjugation restricts the possible nonscalar terms to exact support orbits; the smallest even circuit identifies the first unavoidable weak-rotation loss; and the CNOT construction separates that property from ordinary distance.

Use the [map's final synthesis](MILESTONES.md) alongside [physical interpretation](../science/interpretation.md) to retain the operational boundaries. The code, physical tensor-product structure, and ideal full measurement are fixed. Recovery knows the specified unitary and acts afterward with unrestricted resources. The asymptotic theorem fixes the code before taking $\theta\to0$. The growing family encodes two logical qubits and supplies exact values of $r$, with no claim of a constant-rate family or uniform threshold.

The single anchor supplies the language of stabilizer quantum error correction. Established research supplies the outcome-assisted recovery optimum and the even-repetition mechanism. This repository proves the exact support classification, its circuit and exponent consequences, and the equal-parameter family. The [related-work account](../references/RELATED_WORK.md) records those boundaries, while the topic-based [scientific reading route](../README.md#read-the-science) leads directly to every complete proof.

