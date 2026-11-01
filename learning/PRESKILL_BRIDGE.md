# From Preskill's QEC foundations to syndrome reversibility

An error syndrome is useful because it can reveal what happened to a physical system while preserving the quantum information encoded in it. This project asks exactly when that separation survives coherent rotations on every physical qubit. The answer depends on how stabilizer supports and logical-operator supports fit together. The final sections connect that readout property to protected logical control on the same code.

This bridge follows the [30-milestone map](MILESTONES.md). It develops the connections needed to read the complete proofs under [science/](../science/model.md), with a small measurement calculation before the general theorems. The scientific pages remain the authoritative statements of assumptions, quantifiers, and proofs.

## 1. The single educational anchor

This project and [Error-Transparent Locality](https://github.com/GoGoKo699/Error-Transparent-Locality) use one common educational anchor: **John Preskill, *Quantum Information*, Chapter 7, “Quantum Error Correction,” updated March 2026**. The [author's course index](https://preskill.caltech.edu/ph229/index.html) links the [official chapter](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf). Printed page numbers below are six smaller than PDF viewer page numbers.

The common foundation starts with exact correction, information leakage and a reference system, then develops stabilizers, logical operators, the five-qubit code and concatenation. This bridge follows those ideas into the readout experiment; the [companion gate bridge](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/learning/PRESKILL_BRIDGE.md) follows them into protected dynamics.

| Verified chapter location | Printed pages | Foundation used here |
|---|---:|---|
| [§7.2, “Conditions for quantum error correction”](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf#page=10), including Eq. (7.17) | 4–7 | Correctable error spans, recovery and degeneracy |
| [§§7.3.1–7.3.3](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf#page=14) | 8–10 | Distinguishability, leakage of logical information and reference decoupling |
| [§7.4, “Some general properties of QECC's”](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf#page=18), especially §§7.4.1 and 7.4.3 | 12–14 | Pauli expansions, weight, distance and detection |
| [§§7.6–7.8](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf#page=24) | 18–24 | Classical linear codes, CSS checks and encoded superpositions |
| [§§7.10.1–7.10.2](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf#page=33), including Eqs. (7.116) and (7.129)–(7.132) | 27–31 | Stabilizers, the real-Pauli convention and binary symplectic notation |
| [§7.10.4, “Encoded qubits”](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf#page=39) | 33–34 | The logical Pauli quotient and encoded operators |
| [§7.11, “The 5-Qubit Code”](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf#page=40), including Eq. (7.150) | 34–37 | The standard five-qubit seed and syndrome measurement |
| [§7.17, “Concatenated codes”](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf#page=54) | 48–50 | Inner and outer encoding and multiplicative distance bounds |

The chapter's fault-tolerance discussions are background context. The present readout measurement is ideal, and the weak-angle limit fixes the code before taking the angle to zero. The local Clifford and CNOT rules needed below are derived in §3; channel metrics and the record-assisted recovery optimizer are supplied in §§4–6. The chapter is not assigned as a source for those additional theorems.

The following additions are supplied here rather than assigned to the anchor: the precise optimization over recovery after a recorded measurement; the fidelity and diamond-distance estimates needed for that optimization; the support/circuit theorems proved in this repository; and the Hamiltonian resource and tensor-increment arguments connecting readout to protected gates. The first addition uses established results of Gregoratti and Werner. The even-repetition mechanism used later is due to earlier coherent-error work, including Huang, Doherty and Flammia. These primary research sources are attributed in [Related work](../references/RELATED_WORK.md); they do not constitute additional educational anchors.

## 2. Route through the milestones

The map has 12 background/preliminaries milestones and 18 narrative/innovation milestones. The first phase combines individual physical and mathematical objects into the complete operational question. Most of the second phase develops the support classification, sharp exponent and same-parameter family separation. Its final milestones connect those results to gate-resource constraints and fixed-code coexistence.

| Milestones | Connection developed here |
|---|---|
| M01–M06 | Encoding, Pauli labels, checks, logical operators and coherent rotations: §§3–5 |
| M07–M12 | Recorded measurements, recovery and performance measures: §§4–6 |
| M13–M14 | From the full syndrome to the exact support theorem: §7 |
| M15–M16 | From support collisions to even logical circuits: §8 |
| M17 | A uniform bound that includes rare outcomes: §9 |
| M18–M19 | Exact repetition witness and matching exponent: §10 |
| M20–M24 | Five-qubit blocks, concatenation, one physical CNOT and its finite channel: §11 |
| M25–M26 | Separate gate/readout models and their joint locality theorem: §§12–13 |
| M27–M29 | Optimal protected axes, fixed-code entanglers and resource scope: §14 |
| M30 | The complete scientific account: §15 |

## 3. Translate the notation before comparing statements

A code stores $k$ logical qubits inside $n$ physical qubits. Its logical Hilbert-space dimension is $D=2^k$. An encoding isometry $V$ maps logical states into the physical codespace; $V^\dagger V=I_D$ and $P=VV^\dagger$ is the projector onto that subspace.

| Concept | Convention in this repository |
|---|---|
| Stabilizer | The actual physical group consists of signed Hermitian Paulis fixing the code. For support counting, $S$ denotes its labels modulo global phases. |
| Centralizer | $N$ is the phase-free set of Paulis commuting with every stabilizer, corresponding to the anchor's symplectic orthogonal space $S^\perp$. |
| Logical Pauli | A label in $N\setminus S$ acts nontrivially on the code. A bar denotes its induced logical operator. |
| Support and weight | Support is the set of nonidentity sites; weight is its size. Equal weight does not imply equal support. |
| Syndrome | $s$ is a binary string recording the eigenvalues of the independent checks. The actual projectors retain the physical signs. |
| Rotation scale | A local rotation is $\exp(-i\theta h_j/2)$. For a Pauli axis, its error amplitude is $\sin(\theta/2)$. |
| Recovery information | Each specified physical unitary may have its own recovery, applied after the prescribed measurement. |
| Error convention | Diamond error means one half of the diamond norm of the difference from the identity channel. |

### Convert the real-Pauli convention before forming observables

Preskill changes convention in §7.10.1: Eq. (7.116) uses the real, anti-Hermitian matrix $`Y_{\mathrm{P}}=ZX=i\sigma_y`$, whose square is $-I$. This repository uses the Hermitian Pauli $Y=\sigma_y$, with $Y^2=I$. The support and symplectic labels survive this conversion, but physical operator signs and phases must be translated.

In the chapter's binary notation, a representative is $Z(\alpha)X(\beta)$. If $\nu$ sites have both bits set, the canonical Hermitian tensor with the same labels is

```math
Q_{\mathrm{H}}=(-i)^\nu Z(\alpha)X(\beta).
```

Every overlapping site contributes $ZX=iY$, so multiplication by $(-i)^\nu$ converts all of them to Hermitian $Y$ factors. For a signed chapter word $\epsilon Z(\alpha)X(\beta)$, the same physical operator is $`\epsilon i^\nu Q_{\mathrm{H}}`$. A chapter stabilizer has even $\nu$, so this factor supplies a real physical sign. This sign must be retained before defining its stabilizer eigenspace. Binary support alone cannot supply the physical sign. Conjugation before conversion is always $E^\dagger B E$; the abbreviated $EBE$ is valid for the Hermitian Pauli errors used here, whose squares are identity.

Preskill's code projector is $\Pi$, whereas this repository writes $P=VV^\dagger$. The symbol $V$ is our encoding isometry; in §7.2 of the chapter it denotes a recovery isometry. The chapter's phase-free space $S^\perp$ is our centralizer $N$, and its quotient $S^\perp/S$ gives the logical Pauli labels $N/S$. The companion gate bridge uses $K$ where this bridge initially uses $D=2^k$ for logical dimension; both use the same Hermitian physical tensors, encoding isometry and code projector.

Discarding phases is legitimate for support and binary commutation calculations because those calculations retain the symplectic commutation relation separately. It is not legitimate when forming the physical projectors. For example, $+XX$ and $-XX$ specify different positive-eigenvalue codespaces even though their phase-free support labels agree. The [model](../science/model.md#1-code-phases-and-syndrome) keeps these two roles separate.

The measurement order is fixed throughout: encode, apply $U=\bigotimes_jU_j$, measure all original checks, then recover using the complete classical record and surviving quantum state. In symbols, the recovery quantifier is

```math
\text{for every specified }U,\quad
\text{there exists a recovery }\mathcal R_U.
```

The quantifier does not supply one universal decoder for an unknown rotation. It also places recovery after readout. Applying a known $U^\dagger$ before readout would answer a different, immediate reversal problem.

### Local derivation of the Clifford and CNOT rules

The support argument needs only a few concrete Pauli conjugations. Let $F=(X+Z)/\sqrt2$ be the Hadamard matrix and $R=\mathrm{diag}(1,i)$. Matrix multiplication gives

```math
FXF^\dagger=Z,\qquad FZF^\dagger=X,
\qquad RXR^\dagger=Y,\qquad RZR^\dagger=Z.
```

Together with $Y=iXZ$, these rules determine the conjugation of every Pauli and allow any nonidentity one-qubit axis to be sent to any other, up to sign. A unitary that maps Pauli operators to Pauli operators under conjugation is called a Clifford unitary. Products of the one-qubit changes above preserve each physical support. This is the only local-Clifford fact needed in the support classification.

For a control $c$ and target $t$, define CNOT by $C\lvert a,b\rangle=\lvert a,a\mathbin\oplus b\rangle$. Its action on the computational basis gives $C=C^\dagger$ and

```math
CX_cC^\dagger=X_cX_t,\qquad CZ_cC^\dagger=Z_c,
\qquad CX_tC^\dagger=X_t,\qquad CZ_tC^\dagger=Z_cZ_t.
```

Unlike a product of one-qubit changes, this two-qubit conjugation can change support. In §11 it conjugates the code and checks to define a second code. In §14 a protected *logical* entangler instead preserves the original codespace and checks. The two operations have different roles. These explicit rules supplement the chapter's stabilizer and encoding examples; no additional educational source is required.

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

for a fixed code with finite $r$ and for $|\theta|\le1$. A zero-probability effect is the zero matrix by positivity, so it is covered as well. This mechanism is developed in [the upper-bound proof](../science/error_exponent.md#2-a-finite-code-collision-polynomial). The estimate gives one explicit code-dependent constant; it does not bound that constant uniformly over larger codes.

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

This positive coefficient and the preceding upper bound establish the sharp worst-case exponent for each fixed code with finite $r$:

```math
\sup_{\boldsymbol h}e_{\mathrm{opt}}(U_{\boldsymbol h}(\theta))
=\Theta(|\theta|^r),
\qquad
\sup_{\boldsymbol h}\varepsilon_{\mathrm{opt}}(U_{\boldsymbol h}(\theta))
=\Theta(|\theta|^r).
```

For $r=\infty$, both errors vanish exactly for every product unitary. At $r=2$ the toy example gives $p_2(q)=q$; at $r=4$ it gives $p_4(q)=3q^2-2q^3$. The repetition mechanism is inherited; the code-wide statement that every code with finite $r$ contains this witness is the new connection. The [complete lower bound](../science/error_exponent.md#3-an-embedded-even-repetition-witness) distinguishes this exact witness equality from the two generally different recovery metrics.

## 11. From the five-qubit code to the one-CNOT separation

Preskill §7.11, Eq. (7.150), supplies the standard five-qubit stabilizer. The [family proof](../science/code_families.md#2-the-five-qubit-building-block) checks its distance and the required parity properties directly. Its stabilizers all have even weight, while every representative of a nontrivial logical Pauli has odd weight. These are properties of complete stabilizer and logical cosets, not only their displayed generators.

Preskill §7.17 defines concatenation by replacing each outer physical qubit with an inner encoded qubit. A nontrivial outer logical must act on at least three inner blocks; each such action costs the inner distance. Minimum representatives attain that product, yielding $[[5^t,1,3^t]]$ for every integer $t\ge1$. The parity property also persists: an even number of odd-weight inner logical factors produces an even-weight outer stabilizer, and odd logical representatives remain odd. These familiar coding ingredients support the new comparison.

Now take two such blocks. A matching stabilizer/logical support would need to match separately inside each block. On at least one block the logical is nontrivial and has odd weight, whereas the stabilizer restriction has even weight. This is impossible, giving $r=\infty$. Total parity alone would be insufficient: a logical acting nontrivially on both blocks can have even total weight.

Conjugate the code by CNOT from physical qubit $0$ in the first block to physical qubit $5^t$ in the second. The construction chooses a minimum logical in the first block that avoids the control, so one logical of weight $3^t$ survives unchanged. A hypothetical shorter logical after CNOT would inverse-conjugate to a logical of weight $3^t$ in the product code touching both blocks. A minimum-weight product-code logical cannot do that, because its nontrivial logical block already uses the entire weight budget. The ordinary distance stays $3^t$.

Meanwhile, a minimum logical $E_B$ with a $Z$ factor at the target transforms into $Z_0E_B$. Its weight is $3^t+1$. A normalizing subword must select the new control factor exactly when it selects the target factor; otherwise inverse CNOT leaves a detectable single-qubit $Z_0$ in the first block. Minimality of $E_B$ then permits only the empty and full subwords. The transformed operator is an even logical circuit. Since $r$ is even and at least the odd distance $3^t$, this circuit forces $r=3^t+1$.

The complete result is

| Physical code | Common parameters | Syndrome reversibility distance |
|---|---|---|
| Two concatenated five-qubit blocks | $[[2\cdot5^t,2,3^t]]$ | $\infty$ |
| Their conjugate under one physical CNOT | $[[2\cdot5^t,2,3^t]]$ | $3^t+1$ |

For $t=1$, the explicit ten-qubit rotation gives the exact error $3q^2-2q^3$. The [verification suite](../verification/README.md) constructs its physical state and all 256 signed syndrome sectors, then checks the recovered channel. The proof supplies the statement for every concatenation level.

The CNOT changes the code and the observables measured afterward. It is not an additional error in the product-unitary noise model. An entangling Clifford can change support relative to physical qubits, so it can change this invariant even when the conventional parameters are preserved.

## 12. Two tasks on the same code

The readout experiment has no Hamiltonian implementing a logical gate: a specified product unitary acts, the original checks are measured, and an informed recovery follows. A protected-gate experiment instead specifies a physical Hamiltonian $H$ with logical generator $h_L=V^\dagger HV$. It requires

```math
HEV=EV(h_L+\lambda_E I)
\quad\text{for every tensor Pauli }E
\text{ with }\mathrm{wt}(E)\le t,
\qquad \lambda_I=0.
```

The equality is in the full physical Hilbert space, so it includes absence of leakage from the protected error image $EV\mathcal H_L$. The real scalars $\lambda_E$ permit consistent energy shifts of error sectors. The constructions below satisfy the stronger, strict identity with every shift zero.

To see the operational role, suppose one of the protected errors occurs at an unknown time $u$ in a gate of duration $\tau$. The identity gives

```math
e^{-iH(\tau-u)}E e^{-iHu}V
=e^{-i\lambda_E(\tau-u)}EVe^{-ih_L\tau}.
```

After ideal correction of $E$, the intended logical evolution remains, apart from an overall phase for that error. For multiple faults and time-dependent controls, the corresponding identity must hold at every instant and the accumulated physical error must stay in the protected set. The complete hypotheses are in the [gate model](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/model.md) and [control account](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/control_robustness.md).

The native resource model expands $H$ in the fixed bare Hermitian tensor-Pauli basis:

```math
H=\sum_Q h_QQ,\qquad
\ell=\max_{h_Q\ne0}\mathrm{wt}(Q),\qquad
J=\sum_{Q\ne I}|h_Q|.
```

Thus $\ell$ counts how many physical qubits one interaction term can involve, while $J$ counts the total absolute strength of nonidentity native interactions. The symbol $\ell$ avoids confusing gate interaction order with the readout distance $r(S)$, also written $r_S$. For a pure logical Pauli generator $h_L=\Omega P+cI$, its logical half-splitting is $\omega=|\Omega|$.

These definitions extend the anchor's language of encoded operators. The Hamiltonian identities and native coefficient bounds are supplied by the scientific sources. Two properties of one code can now be compared without substituting gate faults for the known-product-unitary readout channel.

## 13. Why minimum gate order forces readout loss

The [joint theorem](../science/gate_readout_tradeoff.md#2-joint-theorem) assumes odd detection distance $d\ge3$, a nonscalar protected logical generator, and protection through weight $t$ with $t\ge1$ and $2t+1\le d$. It permits the scalar error-sector shifts above. Its conclusion is

```math
\ell=d+t\quad\Longrightarrow\quad r_S=d+1.
```

The gate theorem supplies the general interaction bound $\ell\ge d+t$. Equality has a consequence beyond resource use: it forces an even logical circuit and therefore a weak-rotation readout witness.

First remove physically fixed qubits. A weight-one stabilizer fixes its site in a Pauli eigenstate. Local Clifford changes and permutations put all such sites into a product state $|a\rangle$, leaving a reduced code with the same distance and $r_S$. In any equal-support pair, the fixed-site factors agree and cancel from $gL$, so they cannot change the readout invariant. Taking the partial matrix element $\langle a|H|a\rangle$ preserves the logical generator and transparency on the retained sites and cannot increase interaction order. The reduced code has no weight-one stabilizer.

Group the reduced Hamiltonian by physical Pauli weight and remove scalar parts after logical compression. If $A_w$ is its traceless logical component at weight $w$ and $A=[h_L]_0\ne0$, the gate's operator moment identities read

```math
A_w=0\quad(w<d),\qquad
\sum_wA_w=A,\qquad
\sum_ww^jA_w=0\quad(1\le j\le t).
```

At order $d+t$, only the $t+1$ weights $d,\ldots,d+t$ can contribute. If at most $t$ of those components were nonzero, the polynomial $q(x)=\prod_{w\in\mathcal W}(1-x/w)$ over their weights would have degree at most $t$ and constant term one. The moments would give $\sum_wq(w)A_w=A$, although every term on the left vanishes. Therefore every allowed weight contributes, including $d+1$. A nonzero traceless compression at that weight requires a nontrivial logical Pauli $E$ of weight $d+1$.

Distance and the removal of weight-one stabilizers make $E$ a logical circuit. A proper stabilizer subword would remove at least two sites and leave a logical below distance. A proper nontrivial logical subword already needs at least $d$ sites and would leave a weight-one normalizer, which the reduced code excludes. No proper normalizer subword remains. Since $d+1$ is even, the circuit theorem gives $r_S\le d+1$; parity and distance give the reverse bound.

Consequently,

```math
r_S>d+1\quad\Longrightarrow\quad\ell\ge d+t+1,
```

including the case $r_S=\infty$. This is a necessary resource constraint. It does not construct a gate at the next order. For even $d$, the readout theorem already gives $r_S=d$ independently of any Hamiltonian. The complete joint account also resolves which logical components of a protected generator admit witnesses; different components can require different product rotations.

## 14. Protected control with the original readout unchanged

### Exact protected axes on the concatenated family

The family $Q_m=[[5^m,1,3^m]]$ from §11 has only odd-weight nontrivial logical representatives. Writing $d=3^m$, this removes alternating interaction weights from the moment problem and strengthens the necessary gate order to $d+2t$. Exact constructions attain the following one-axis resources:

| Protected error weight $t$ | Depth | Minimum order $\ell$ | Minimum $J/\omega$ at that order |
|---:|---|---:|---|
| $1$ | $m\ge1$ | $d+2$ | $J_1(d)=d+1$ |
| $2$ | $m\ge2$ | $d+4$ | $J_2(d)=(4d^2+21d+8)/8$ |
| $3$ | $m\ge2$ | $d+6$ | $J_3(d)=(4d^3+51d^2+170d+24)/24$ |

The sources prove the strict identities and matching resource lower bounds in the declared native model: [one-error gates](../science/minimum_readout_safe_gates.md), [two-error gates](../science/growing_two_error_gates.md), and [three-error gates](../science/growing_three_error_gates.md). The budget optimum is at the specified minimum order. It does not exclude a smaller budget if higher-order interactions are allowed.

A useful common mechanism is the nested hierarchy of strict gates $L_{a,\beta}$ for a logical axis $\beta$ and protection orders $a=0,1,2,3$. Put $L_{-1,\beta}=0$ and

```math
\Delta_{a,\beta}=L_{a,\beta}-L_{a-1,\beta}.
```

On an error image lighter than $a$, the two gates have identical action. Their difference therefore obeys

```math
\Delta_{a,\beta}EV=0
\qquad\text{if }\mathrm{wt}(E)<a.
```

This is a full error-image identity, rather than a scalar expectation-value condition. It makes high total-degree terms vanish when the total physical error weight is small. The three-error construction lifts the complete hierarchy by tensoring increments whose index sum is at most three. The proof tracks coefficient signs to attain the coupling bound after collecting native words; arbitrary independently chosen optimal gates need not supply that sign compatibility.

### Logical entanglement on fixed product blocks

Take $M$ unchanged $Q_m$ blocks with $m\ge2$. A pure logical Pauli string $P$ of logical weight $p\ge1$ has minimum physical representative weight

```math
D=pd.
```

Every active block must carry a nontrivial logical representative of weight at least $d$, and tensor products of minimum representatives attain the bound. The ordinary distance of the product code remains $d$, since a logical can act on only one block. Here $D$ denotes a target-coset weight; earlier sections used $D=2^k$ for logical Hilbert-space dimension. The target-specific weight must not be substituted for the code distance.

For the active encoded blocks $A(P)$, define

```math
H_P=
\sum_{\substack{a_i\ge0\ (i\in A(P))\\\sum_i a_i\le3}}
\bigotimes_{i\in A(P)}\Delta_{a_i,P_i},
```

with identity factors on inactive encoded blocks. The full tensor product of the $L_{3,P_i}$ gates has the desired logical action. An omitted term has $\sum_i a_i>3$, whereas every protected physical error has total weight at most three. At least one factor then has $a_i>\mathrm{wt}(E_i)$ and vanishes. Truncation therefore preserves the full strict identity $H_PEV=EVP$ for every protected error, including entangled logical inputs.

Each retained term has order at most $pd+2\sum_i a_i\le D+6$. Compatible signs make its native budget exactly $J_3(D)$. The [target-coset moment lower bound](../science/protected_entangling_gates.md#4-exact-coupling-lower-bound-at-that-order) permits arbitrary physical Hamiltonians, including other logical cosets and scalar sector shifts, and proves

```math
h_L=\Omega P+cI,\quad\Omega\ne0:
\qquad\ell_{\min}=D+6,\qquad
J_{\min}(\ell=D+6)=|\Omega|J_3(D).
```

For $m=2$ and $P=ZZ$, two blocks form a $[[50,2,9]]$ code. The target minimum is $D=18$, the optimal order is $24$, and $J/\omega=1789$. The logical unitary $e^{-i\pi Z_L\otimes Z_L/4}$ sends $|++\rangle_L$ to a maximally entangled state. Together with the six local logical Pauli directions, this generator yields all fifteen nonidentity two-qubit Pauli directions under commutators, giving $\mathfrak{su}(4)$ control. Protection of a time-dependent combination requires the strict identity at every instant, as in §12.

Every stabilizer still has an even-weight restriction in every encoded block. Every nontrivial product-code logical has an odd-weight restriction in at least one block. Hence their supports cannot match and $r_S=\infty$, even when the target $ZZ$ string has even global weight. A protected logical gate preserves this codespace and its original checks. The physical CNOT in §11 instead changed the code and checks, which is why the two entangling operations have different consequences for the readout invariant.

### Certificates, native resources and operational boundaries

The finite 50-qubit construction has eighteen ordered class pairs, but its dense expansion has 41,032,000,000 distinct native Pauli words. The exact verifier covers all 540,376 physical Pauli errors through weight three using complete response profiles; it does not enumerate that dense interaction list. Compact class descriptions, native term count, coefficient bit length and many-body interaction synthesis are separate resources. The [construction and finite certificate](../science/protected_entangling_gates.md#5-the-first-protected-entangler) state each count and its verification scope.

The pure-string formula is a target-specific optimum. It does not give an optimal weighted sum for an arbitrary mixture of logical Pauli terms at one common interaction cap: components of lower target weight may use the additional allowed locality. The [target-specific scope statement](../science/protected_entangling_gates.md#4-exact-coupling-lower-bound-at-that-order) retains this distinction between axis constructions and a broader optimization problem.

These constructions establish coexistence of exact protected control and $r_S=\infty$ on the same fixed code. The gate statements use their protected-fault and ideal-recovery assumptions. The readout theorem separately uses each specified known product unitary, full ideal original-syndrome measurement and informed post-measurement recovery. They do not establish a joint noisy gate/readout threshold, noisy extraction tolerance, efficient recovery or geometric hardware locality. The native Hamiltonians must satisfy the operator identity during evolution; a large exact interaction sum does not justify replacing it by an arbitrary sequence of native pulses while faults can occur.

## 15. Assemble the complete statement

The main readout chain is now explicit: the prescribed measurement induces logical effects; their scalarity decides exact reversibility; product conjugation restricts nonscalar terms to exact-support orbits; the minimum even circuit identifies the unavoidable weak-rotation exponent; and one physical CNOT separates this property from ordinary code parameters at arbitrarily large distance.

The joint consequence uses that circuit theorem with the gate's operator moments. At odd distance, saturation of the general protected-gate order forces the smallest finite readout distance. The concatenated five-qubit family instead has $r_S=\infty$ and attains its stronger family-specific orders for one, two and three protected errors. Tensor increments then produce optimal pure-string gates, including logical entanglers, while keeping the code and its original syndrome measurement fixed.

The [map's final synthesis](MILESTONES.md#m30-the-complete-scientific-account) and [physical interpretation](../science/interpretation.md) retain the operational boundaries. In the readout theorem, recovery knows the specified product unitary and acts only after the fixed ideal full measurement. The asymptotic theorem fixes the code before taking $\theta\to0$; its constants can depend on that code. The same-parameter family encodes two logical qubits and supplies exact values of $r$, without a constant-rate or uniform-threshold claim. Gate coexistence retains the independent error-transparent Hamiltonian model.

The single anchor supplies the stabilizer language. Established research supplies the outcome-assisted recovery optimum and even-repetition mechanism. The scientific pages prove the support classification, circuit and exponent consequences, equal-parameter family, and shared gate-resource and coexistence results. [Related work](../references/RELATED_WORK.md) records the attribution boundaries; the [learning route](README.md) points directly to the relevant explanations and proofs.
