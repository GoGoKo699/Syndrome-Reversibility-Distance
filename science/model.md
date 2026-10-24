# The measured channel and optimal recovery

The object of study is the information retained after a prescribed syndrome measurement. A product of single-qubit unitaries acts on a stabilizer code, the original full syndrome is measured, and recovery then receives both the physical output and its classical syndrome. The recovery may know the complete unitary. This ordering isolates irreversibility introduced by the measurement itself.

The [support criterion](support_criterion.md) classifies exact reversibility for every product unitary. The [error exponent](error_exponent.md) quantifies the smallest unavoidable loss when that classification fails. This page fixes the channel and proves the recovery facts used by both results.

## 1. Code, phases, and syndrome

Let the physical Hilbert space be $\mathcal H_P=(\mathbb C^2)^{\otimes n}$ and the logical space be $\mathcal H_L=\mathbb C^D$, with $D=2^k$ and $k\geq1$. A stabilizer code is specified by an abelian group of Hermitian Pauli operators that does not contain $-I$. Choose $m=n-k$ independent signed generators $g_1,\ldots,g_m$, and write

```math
g(a)=\prod_{\ell=1}^{m}g_\ell^{a_\ell},
\qquad a\in\mathbb F_2^m.
```

Its codespace is the simultaneous $+1$ eigenspace. Fix an encoding isometry $V:\mathcal H_L\to\mathcal H_P$ and let $P=VV^\dagger$.

Two distinct Pauli conventions serve different purposes:

| Object | Convention |
| --- | --- |
| Physical checks $g(a)$ and projectors | Actual signed Hermitian operators, including the signs that specify the code |
| $S$, the stabilizer label space | The image of the physical stabilizer group after quotienting by all global Pauli phases |
| $N$, its Pauli centralizer | Phase-free Pauli labels commuting with every stabilizer |
| $N\setminus S$ | Nontrivial logical Pauli labels; a change of global phase never creates a new logical operator |
| $\mathrm{supp}(Q)$ and $\mathrm{wt}(Q)$ | Nonidentity sites and their number, independent of global phase |

We use the same letter for a phase-free label and a Hermitian representative when the distinction does not affect an equation. Operator equations involving the code retain their actual signs. In particular, a canonical tensor Pauli whose label belongs to $S$ may compress to either $+I_D$ or $-I_D$.

The full syndrome is $s\in\mathbb F_2^m$, with projector

```math
\Pi_s=2^{-m}\sum_{a\in\mathbb F_2^m}(-1)^{s\cdot a}g(a).
\tag{1}
```

Each syndrome subspace has dimension $D$, and $\sum_s\Pi_s=I_P$. Choose a Pauli isometry $W_s:\mathcal H_L\to\Pi_s\mathcal H_P$ identifying that sector with the logical space. For example, $W_s=T_sV$ for a Pauli $T_s$ having syndrome $s$.

## 2. The prescribed order of operations

For a specified product unitary

```math
U=\bigotimes_{j=1}^{n}U_j,
```

the measured channel is

```math
\mathcal N_U(\rho)
=\sum_s |s\rangle\langle s|\otimes
\Pi_sUV\rho V^\dagger U^\dagger\Pi_s.
\tag{2}
```

The classical register has been dephased: coherent superpositions of different syndrome outcomes are not retained. The recovery acts after (2), may depend on $U$, and is otherwise an arbitrary completely positive trace-preserving map from the joint output to $\mathcal H_L$. Because the input to recovery is block diagonal in $s$, its action is equivalent to a family of trace-preserving recovery channels, one on each syndrome sector. A channel that also manipulates the classical register gives no additional action on these block-diagonal inputs.

Under the identifications $W_s$, the instrument has one logical Kraus operator per syndrome,

```math
K_s=W_s^\dagger\Pi_sUV,
\qquad
M_s=K_s^\dagger K_s=V^\dagger U^\dagger\Pi_sUV.
\tag{3}
```

Thus $M_s\geq0$, $\sum_sM_s=I_D$, and the outcome probability on logical state $\rho$ is $\mathrm{Tr}(M_s\rho)$. Changing $W_s$ only left-multiplies $K_s$ by a unitary and can be absorbed into recovery.

The recovery quantifiers are essential:

```math
\text{for every specified }U\text{, there exists a recovery }\mathcal R_U.
\tag{4}
```

We do not require one recovery to work for an unknown family of rotations. We also do not allow an operation between $U$ and the prescribed syndrome measurement. If the known inverse $U^\dagger$ could be applied there, the physical rotation could be undone before any measurement. Equations (2)–(4) ask which information is still available once the original measurement has already occurred.

## 3. Performance measures

Let $|\Omega_D\rangle=D^{-1/2}\sum_{i=1}^{D}|i\rangle|i\rangle$ be maximally entangled between the logical system and a reference. For a logical channel $\Lambda$, its entanglement fidelity at the maximally mixed state is

```math
F_e(\Lambda)=
\langle\Omega_D|(\Lambda\otimes\mathrm{id})(|\Omega_D\rangle\langle\Omega_D|)|\Omega_D\rangle.
```

Define

```math
\begin{aligned}
F_{\mathrm{opt}}(U)&=\sup_{\mathcal R}F_e(\mathcal R\mathcal N_U),\\
e_{\mathrm{opt}}(U)&=1-F_{\mathrm{opt}}(U),\\
\varepsilon_{\mathrm{opt}}(U)&=
\inf_{\mathcal R}\frac12\Vert \mathcal R\mathcal N_U-\mathrm{id}_L\Vert _\diamond.
\end{aligned}
\tag{5}
```

The diamond norm includes an arbitrary reference system. Both errors lie in $[0,1]$. For any recovery,

```math
1-F_e(\mathcal R\mathcal N_U)
\leq\frac12\Vert \mathcal R\mathcal N_U-\mathrm{id}_L\Vert _\diamond.
\tag{6}
```

Indeed, the diamond norm dominates the trace distance obtained with $|\Omega_D\rangle$, and the two-outcome test onto $|\Omega_D\rangle$ distinguishes that output from the ideal output with probability difference $1-F_e$. In particular, $e_{\mathrm{opt}}(U)\leq\varepsilon_{\mathrm{opt}}(U)$.

## 4. Exact optimal entanglement fidelity

The next proposition is the optimal record-assisted fidelity result of Gregoratti and Werner, applied to the syndrome instrument. Their Proposition 8 and Eqs. (16)–(20) give the same optimum and polar-corrected channel. We include the proof to fix the present conventions and explicitly cover singular branches; see [the source comparison](../references/RELATED_WORK.md#gregorattiwerner-record-assisted-recovery).

**Proposition 1.** For the measured channel (2),

```math
F_{\mathrm{opt}}(U)
=\frac1{D^2}\sum_s\left(\mathrm{Tr}\sqrt{M_s}\right)^2.
\tag{7}
```

A syndrome-dependent polar recovery attains this value and has logical Kraus operators $A_s=\sqrt{M_s}$.

**Proof.** Take a polar decomposition $K_s=Q_sA_s$. Since $K_s$ is square, its polar partial isometry can be extended to a unitary $Q_s$ on the whole $D$-dimensional space. This extension is arbitrary on the kernel of $A_s$ and does not change $K_s=Q_sA_s$. When $K_s=0$, any unitary extension is suitable.

Applying $Q_s^\dagger$ after outcome $s$ produces the channel

```math
\Phi_U(\rho)=\sum_s A_s\rho A_s,
\qquad
\sum_sA_s^2=I_D.
\tag{8}
```

Its entanglement fidelity is the right-hand side of (7). To prove optimality, absorb $Q_s$ into an arbitrary recovery and denote that branch's recovery Kraus operators by $B_{s,\alpha}$, with $\sum_\alpha B_{s,\alpha}^\dagger B_{s,\alpha}=I_D$. Diagonalize $A_s$ with nonnegative eigenvalues $a_i$. The vectors

```math
v_i=\big((B_{s,\alpha})_{ii}\big)_\alpha
```

satisfy

```math
\Vert v_i\Vert _2^2
=\sum_\alpha|(B_{s,\alpha})_{ii}|^2
\leq\sum_{\alpha,j}|(B_{s,\alpha})_{ji}|^2=1.
```

Therefore

```math
\sum_\alpha|\mathrm{Tr}(B_{s,\alpha}A_s)|^2
=\left\Vert \sum_i a_iv_i\right\Vert _2^2
\leq\left(\sum_i a_i\right)^2.
\tag{9}
```

Summing the branch contributions to entanglement fidelity proves (7). The proof makes no invertibility assumption and includes zero-probability branches. $\square$

Polar recovery is optimal for entanglement fidelity. We do not assume it is diamond-optimal for every $U$; it supplies the upper bound used in [the exponent theorem](error_exponent.md). For the explicit repetition witness there, its achievable error can also be realized by a Pauli recovery, and the diamond optimum is determined exactly.

## 5. Exact reversibility and logical information in the syndrome

**Proposition 2.** The following are equivalent for a fixed $U$:

1. A recovery after (2) restores every logical state, including its entanglement with a reference.
2. $F_{\mathrm{opt}}(U)=1$.
3. Every syndrome effect is scalar: $M_s=p_sI_D$ for some $p_s\geq0$.
4. Each nonzero $K_s$ is a scalar multiple of a unitary between the logical space and its syndrome sector.

**Proof.** Exact restoration implies unit entanglement fidelity. Conversely, a channel with unit fidelity on $|\Omega_D\rangle$ has the identity Choi state and hence is the identity channel. For (7), Cauchy–Schwarz gives

```math
\left(\mathrm{Tr}\sqrt{M_s}\right)^2
\leq D\mathrm{Tr}M_s,
```

with equality precisely when all eigenvalues of $\sqrt{M_s}$ coincide. The sum of the upper bounds is $D^2$ because $\sum_sM_s=I_D$. Consequently (7) is one exactly when each $M_s$ is scalar. Polar decomposition then proves the final equivalence. $\square$

Here $p_s=D^{-1}\mathrm{Tr}M_s$ and $\sum_sp_s=1$. The scalar condition states that the syndrome distribution is independent of the logical input. Because each recorded branch has one Kraus operator, this condition is also sufficient for full restoration. A generic instrument with unrecorded Kraus multiplicity would require additional conditions; the one-Kraus structure of ideal projective syndrome readout is used throughout.

This branchwise reversibility criterion is an inherited error-correction fact, not the contribution of the support theorem. Gregoratti and Werner's Proposition 2 gives the record-assisted criterion that each branch effect is scalar. Liu and Chen's Sec. II.4, Proposition 1, gives a recent syndrome-measurement specialization; their Clifford analysis also supplies the logical-projective interpretation. The proofs above make the conventions and metric optimization explicit. See [the related-work comparison](../references/RELATED_WORK.md) for attribution and the distinction between these ingredients and the support classification.

## 6. Scope of the model

The following boundaries are part of the mathematical statements:

- The physical error is a specified product unitary. Averaging over an unknown unitary creates a different channel and generally changes what is recoverable.
- The original complete stabilizer syndrome is measured ideally. Coarse-grained outcomes, noisy extraction circuits, and retaining coherent syndrome amplitudes define different instruments.
- Recovery has full access to the surviving physical output and classical record. The optimization imposes no computational or locality restriction on recovery.
- The small-angle statements are for a fixed finite code, uniformly over bounded single-qubit axes. Their constants can depend on that code.

Within this scope, irreversibility is a property of the code and the prescribed measurement order. The [code-family construction](code_families.md) shows that the conventional parameters $[[n,k,d]]$ do not determine it.
