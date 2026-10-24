# Exact supports and universal syndrome reversibility

For a fixed stabilizer code, ask whether every specified product unitary can be corrected after measuring the original complete syndrome. The recovery is allowed to depend on that unitary. The answer is determined by the exact supports of the stabilizer and logical Pauli operators.

All definitions and the order of operations are fixed in [the measured-channel model](model.md). In particular, $S$ and $N$ below are phase-free Pauli label spaces, and $N\setminus S$ contains precisely the nontrivial logical Pauli labels.

## 1. The support criterion

**Theorem 1.** The following are equivalent:

1. For every $U=\bigotimes_jU_j$, there is an exact recovery after full ideal measurement of the original stabilizer syndrome.
2. For every such $U$ and every syndrome $s$, the effect $M_s=V^\dagger U^\dagger\Pi_sUV$ is proportional to $I_D$.
3. No stabilizer and nontrivial logical Pauli have the same exact support:

```math
\{\mathrm{supp}(g):g\in S\}
\cap
\{\mathrm{supp}(L):L\in N\setminus S\}
=\varnothing.
\tag{1}
```

The set on the stabilizer side contains the **full stabilizer group**, not only a chosen independent generator list. Exact support also matters: neither equal weight nor support containment can replace equality in (1).

### Proof

The equivalence of statements 1 and 2 is Proposition 2 of [the model](model.md). We prove the equivalence with statement 3.

Fourier inversion of the syndrome projectors gives

```math
\begin{aligned}
\Pi_s&=|S|^{-1}\sum_{a\in\mathbb F_2^{n-k}}
(-1)^{s\cdot a}g(a),\\
g(a)&=\sum_s(-1)^{s\cdot a}\Pi_s.
\end{aligned}
\tag{2}
```

For each fixed $U$, all effects $M_s$ are scalar if and only if every compression $V^\dagger U^\dagger g(a)UV$ is scalar.

Fix a physical stabilizer $g$. Product conjugation preserves its exact support. On each nonidentity site, conjugation sends the local Pauli to a real linear combination of $X,Y,Z$; on each identity site it leaves the identity. Thus

```math
U^\dagger gU
\in\mathrm{span}_{\mathbb R}
\{Q:\mathrm{supp}(Q)=\mathrm{supp}(g)\},
\tag{3}
```

where $Q$ ranges over canonical Hermitian tensor Paulis. Conversely, every tensor Pauli with that exact support occurs, up to its global sign, as a product-Clifford conjugate of $g$. The span of the entire conjugacy orbit is therefore exactly the space in (3).

Compression of a physical Pauli $Q$ has three possibilities:

```math
V^\dagger QV=
\begin{cases}
0,&Q\notin N,\\
\pm I_D,&Q\in S,\\
\overline Q,&Q\in N\setminus S,
\end{cases}
\tag{4}
```

where $\overline Q$ in the last line is a nontrivial logical Pauli, up to a sign, and has zero trace. The first line follows by inserting a stabilizer that anticommutes with $Q$ on either side of the compression. The second is the stabilizer eigenvalue condition. For the third, the quotient $N/S$ is the logical Pauli space; a nonzero element has a logical anticommuting partner, so its action is traceless and cannot be scalar.

If (1) holds, every term in (3) compresses to zero or a scalar. Hence all compressed stabilizers and all $M_s$ are scalar. If (1) fails, choose $g\in S$ and $L\in N\setminus S$ with equal support. A product Clifford can map the physical $g$ to a Hermitian representative of $L$, up to sign. Its compression is a nontrivial logical Pauli by (4), so at least one syndrome effect is nonscalar. This proves the theorem. $\square$

## 2. A collision exposes a logical observable

The necessity proof gives a finite-angle operational witness, not only a failure of a formal condition.

**Corollary 1.** If (1) fails, there is a product Clifford $U$ for which

```math
F_{\mathrm{opt}}(U)\leq\frac12,
\qquad
\varepsilon_{\mathrm{opt}}(U)\geq\frac12.
\tag{5}
```

**Proof.** Choose $U$ so that $U^\dagger gU$ acts as a nontrivial Hermitian logical Pauli $\overline L$ on the code. Absorb a possible sign into $\overline L$. Since each syndrome sector is an eigenspace of $g$, with eigenvalue $\chi_s(g)\in\{+1,-1\}$,

```math
K_s\overline L=\chi_s(g)K_s.
\tag{6}
```

Thus $K_s$ annihilates the opposite logical eigenspace. Both eigenspaces of $\overline L$ have dimension $D/2$, and $\mathrm{rank}M_s\leq D/2$. Cauchy–Schwarz on the support of $M_s$ yields

```math
\left(\mathrm{Tr}\sqrt{M_s}\right)^2
\leq\frac D2\mathrm{Tr}M_s.
```

Summing and applying the optimal-fidelity formula gives $F_{\mathrm{opt}}\leq1/2$. The diamond lower bound follows from $\varepsilon_{\mathrm{opt}}\geq1-F_{\mathrm{opt}}$. $\square$

The syndrome has recorded an eigenvalue of a logical Pauli. Even with perfect knowledge of $U$, subsequent recovery cannot reconstruct the coherence destroyed between its eigenspaces. The interpretation of Clifford branches as logical projections is an established ingredient; the new role of (1) is to characterize exactly which stabilizer codes admit such a product-unitary witness. Attribution is detailed in [related work](../references/RELATED_WORK.md).

## 3. Syndrome reversibility distance

Define

```math
r(S)=
\min_{\substack{g\in S,\ L\in N\setminus S\\
\mathrm{supp}(g)=\mathrm{supp}(L)}}
\mathrm{wt}(gL),
\tag{7}
```

with $r(S)=\infty$ when the set is empty. Multiplication in (7) is phase-free. Because $g$ and $L$ commute globally, a Hermitian representative of their product is available, but its sign has no effect on the weight.

Equivalently, $r$ is the least weight of a nontrivial logical Pauli $E$ for which there exists a stabilizer $g$ satisfying

```math
g_jE_j=-E_jg_j
\quad\text{at every }j\in\mathrm{supp}(E).
\tag{8}
```

To see the equivalence, set $E=gL$ in a support collision. Wherever $g_j$ and $L_j$ differ, they are distinct nonidentity Paulis; their product is the third axis and anticommutes with $g_j$. Wherever they agree, $E_j=I$. Conversely, if (8) holds, $L=gE$ is a nontrivial logical label with the same support as $g$: its active factors remain nonidentity, and its other factors are unchanged.

Since $E\in N$, it commutes globally with $g$. Equation (8) then forces $\mathrm{wt}(E)$ to be even. Also $E$ is a nontrivial logical Pauli, so its weight is at least the ordinary distance $d=\min_{E\in N\setminus S}\mathrm{wt}(E)$. Therefore

```math
r<\infty\quad\Longrightarrow\quad r\text{ is even and }r\geq d.
\tag{9}
```

The [error-exponent theorem](error_exponent.md) gives (7) its operational meaning. Its use does not assume that an optimal recovery can be implemented efficiently.

## 4. Equivalent characterization by even logical circuits

For a fixed Hermitian Pauli word $E=\bigotimes_jE_j$ and a subset $A\subseteq\mathrm{supp}(E)$, let $E_A$ retain the same local factor $E_j$ on $A$ and be identity elsewhere. A **logical circuit** means a nontrivial logical Pauli $E\in N\setminus S$ with

```math
E_A\notin N
\quad\text{for every proper nonempty }
A\subset\mathrm{supp}(E).
\tag{10}
```

The word “circuit” refers to a minimal binary linear dependence among the local syndrome columns. It does not mean a gate circuit. Condition (10) concerns the fixed local labels of $E$; it does not exclude all other Pauli operators with smaller supports.

**Theorem 2.** The reversibility distance is precisely the minimum even weight of a logical circuit:

```math
r(S)=\min\{\mathrm{wt}(E): E\text{ is a logical circuit of even weight}\},
\tag{11}
```

where an empty minimum is infinity.

### Proof: a minimum local witness is a circuit

Suppose $r<\infty$, and choose $E,g$ attaining the minimum in (8). If a proper nonempty subword $E_A$ were in $N$, it would commute globally with $g$. Local anticommutation at all its sites implies $|A|$ is even. If $E_A\in N\setminus S$, it is itself a smaller witness with the same $g$. If $E_A\in S$, then $EE_A$ is a smaller nontrivial logical Pauli supported on the complementary sites, again satisfying (8) with $g$. Both contradict minimality. Hence $E$ is an even logical circuit.

### Proof: an even circuit supplies a stabilizer witness

Let $E$ be a logical circuit of even weight $w$, with active sites $J$. Form the binary matrix

```math
H_E=(v_j)_{j\in J}\in\mathbb F_2^{(n-k)\times w},
\qquad
(v_j)_\ell=
\begin{cases}
1,&g_\ell E_j=-E_jg_\ell,\\
0,&g_\ell E_j=E_jg_\ell,
\end{cases}
\tag{12}
```

where $E_j$ in the commutation relation means the physical single-site Pauli. A subset of columns sums to zero exactly when the corresponding subword normalizes the code. As $E\in N$ and has no proper nonempty normalizer subword,

```math
\ker H_E=\mathrm{span}\{\mathbf1_J\}.
\tag{13}
```

Because $w$ is even, $\mathbf1_J\cdot\mathbf1_J=0$ in $\mathbb F_2$. Therefore

```math
\mathbf1_J\in(\ker H_E)^\perp
=\mathrm{rowspace}(H_E).
```

Choose a row combination equal to $\mathbf1_J$ and take the same product of stabilizer generators. The resulting $g$ anticommutes locally with every $E_j$ on $J$, establishing (8). Each even logical circuit is therefore an admissible witness, proving equality of the two minima and the equivalence of their emptiness. $\square$

## 5. Relation to ordinary distance

**Corollary 2.** Every stabilizer code of even distance has $r=d$.

**Proof.** Let $E$ be a minimum-weight nontrivial logical Pauli. It has no proper nonempty normalizer subword: a logical subword would already be shorter, and a stabilizer subword could be removed to make a shorter logical. Thus $E$ is a circuit. If $d$ is even, Theorem 2 gives $r\leq d$, while (9) gives the reverse inequality. $\square$

For odd $d$, the smallest possible finite value is $d+1$, but $r$ may be larger or infinite. The invariant is not the minimum weight of an arbitrary even logical Pauli: the circuit condition, or equivalently the local stabilizer witness (8), is essential.

Physical-qubit permutations and local Clifford transformations preserve $r$, since they preserve exact supports and the number of sites where two Pauli labels differ. For tensor products,

```math
r(S_1\otimes S_2)=\min\{r(S_1),r(S_2)\},
\tag{14}
```

with the usual infinity convention. To prove the lower bound, restrict a support collision to a block where the logical factor is nontrivial; that block is itself a collision and its product weight cannot exceed the total. For the upper bound, extend any minimizing collision in one block by identity in the other. If neither block admits a collision, the same restriction argument excludes one in the product.

The [one-CNOT family](code_families.md) uses these facts to produce equal-$[[n,k,d]]$ codes with $r=\infty$ and $r=d+1$, respectively, at unbounded distance.
