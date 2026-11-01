# Three-error gates at minimum interaction order

The twice-concatenated five-qubit code $Q_2=[[25,1,9]]$ admits a strict error-transparent logical Pauli Hamiltonian for every physical Pauli error of weight at most three, with

```math
\boxed{\ell=15,\qquad J/\omega=2867/8,\qquad r_S=\infty.}
```

Both the interaction order and the coupling budget at that order are optimal for this fixed code. The lower bound permits arbitrary Hermitian Hamiltonians, including terms outside the supplied construction and consistent scalar energy shifts between error sectors. The construction supplies each logical axis and determines the exact budget for every logical direction in the fixed canonical Pauli frame:

```math
h_L=\Omega_X X+\Omega_Y Y+\Omega_Z Z+cI
\quad\Longrightarrow\quad
J_{\min}=\frac{2867}{8}
\bigl(|\Omega_X|+|\Omega_Y|+|\Omega_Z|\bigr)
\qquad(\ell=15).
```

The new ingredient is an outer weight-five logical class. Third-order transparency forces this class to have nonzero logical contribution. It was absent from the [two-error optimum](two_error_gates.md), whose interaction order was thirteen. This finite theorem supplies the seed for the [all-depth three-error construction](growing_three_error_gates.md), which lifts the nested hierarchy by a truncated tensor expansion.

## 1. Model and minimum interaction order

Use the [fixed-code model](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/model.md). Let $V$ encode one logical qubit into $Q_2$, let $h_L$ be its desired Hermitian generator, and require

```math
HEV=EV(h_L+c_E I)
\qquad\text{for every tensor Pauli }E\text{ with }\mathrm{wt}(E)\le3,
```

with consistent real sector shifts $c_E$. The explicit Hamiltonians below satisfy the stronger equations with every $c_E=0$. Let $\ell$ be the largest nonzero native Pauli weight in $H$, and let $J$ be the sum of absolute nonidentity bare Hermitian Pauli coefficients. Write $A=[h_L]_0$ for the traceless logical generator and $\omega=\lVert A\rVert_\infty$ for its half-splitting.

The [odd-weight obstruction for this code family](minimum_readout_safe_gates.md#6-scope-and-verification) gives $\ell\ge d+2t$. Here $d=9$ and $t=3$, so

```math
\ell\ge15.
```

In detail, all nontrivial logical representatives have odd weight. At order below fifteen only weights $9,11,13$ could contribute to $A$. The three vanishing weight moments required by three-error transparency, together with normalization, need at least four available weights. The construction below reaches order fifteen.

## 2. Seven positive logical-coset classes

Use five inner blocks, each containing five physical qubits. The inner code has positive stabilizer generators

```math
XZZXI,\quad IXZZX,\quad XIXZZ,\quad ZXIXZ,
```

and logical representatives $\overline\beta=\beta^{\otimes5}$ for $\beta\in\{X,Y,Z\}$. A signed positive logical representative is $P_Q=s_Q Q$, where $Q$ is a bare Hermitian Pauli, $s_Q\in\{\pm1\}$, and $P_QV_1=V_1\beta$.

Let $G_3^\beta$ and $G_5^\beta$ be the normalized positive signed-coset averages at weights three and five. They contain ten and six words, respectively. Let $K_4$ be the normalized average of the fifteen nonidentity positive signed stabilizers. Identity is a separate inner choice. These are the same phase-aware ingredients used in the [two-error construction](two_error_gates.md).

For a fixed desired logical axis $\alpha$, choose one outer positive signed logical representative and replace each nonidentity outer factor $\beta$ by $G_3^\beta$ or $G_5^\beta$. A selected identity factor can be replaced by $K_4$. The outer representative's Hermitian sign is retained. Thus every expanded atom acts as positive logical $\alpha$ on the concatenated code, regardless of its bare Pauli sign.

Define seven normalized averages. For the first six classes, average over the ten outer weight-three representatives. Their three active children receive the indicated inner logical weights. Promotions and padding locations are chosen uniformly. For the final class, average over the six outer weight-five representatives and use $G_3$ in every child.

| Class | Inner block weights, up to permutation | Definition | Native words |
|---|---|---|---:|
| $A_0$ | $(3,3,3,0,0)$ | No promotions or padding | 10,000 |
| $A_1$ | $(5,3,3,0,0)$ | Promote one of three active children to $G_5$ | 18,000 |
| $A_2$ | $(5,5,3,0,0)$ | Promote two active children | 10,800 |
| $A_3$ | $(5,5,5,0,0)$ | Promote all three active children | 2,160 |
| $F_0$ | $(3,3,3,4,0)$ | Put $K_4$ in one of two inactive children | 300,000 |
| $F_1$ | $(5,3,3,4,0)$ | Promote one active child and pad one inactive child | 540,000 |
| $B_5$ | $(3,3,3,3,3)$ | Outer weight-five logical; all five children use $G_3$ | 600,000 |

The displayed block weights identify the class, up to permutation. Within a class, the actual inner logical labels and weights identify the root representative and its promotion/padding choices. The complete inner cosets are disjoint, so different choices produce distinct native words. The seven classes are disjoint as well. Their total support is 1,480,960 words per axis.

These counts describe the explicit dense construction. They are not lower bounds on necessary term count.

## 3. Exact construction and all protected error images

Set

```math
H_\alpha=
\frac{715}{16}A_0
-\frac{1755}{16}A_1
+\frac{3645}{32}A_2
-\frac{1485}{32}A_3
-\frac{675}{32}F_0
+\frac{675}{32}F_1
-\frac32 B_5.
\tag{T1}
```

The seven signed masses sum to one. Every signed atom $P_Q$ belongs to the same positive logical $\alpha$ coset, so for every physical Pauli $E$,

```math
P_QEV=\chi(E,Q)EV\alpha,
\qquad\chi(E,Q)\in\{+1,-1\}.
```

This is a full physical-space identity. It preserves the error image rather than merely fixing its projected energy. Consequently, a normalized class $T$ has response $f_T(E)$ defined by

```math
TEV=f_T(E)EV\alpha,
\qquad
f_T(E)=\mathbb E_{Q\sim T}[\chi(E,Q)].
```

Inner tensor factors give exact products of conjugation characters. Enumerating the complete signed five-qubit cosets and the physical Pauli errors through weight three gives the following seven response profiles for each logical axis. Representative errors use one-based site indices, with sites 1--5 in the first bottom block and sites 6--10 in the second.

| Representative/profile | $A_0$ | $A_1$ | $A_2$ | $A_3$ | $F_0$ | $F_1$ | $B_5$ |
|---|---:|---:|---:|---:|---:|---:|---:|
| $I$ | $1$ | $1$ | $1$ | $1$ | $1$ | $1$ | $1$ |
| $X_1$ | $13/25$ | $31/75$ | $23/75$ | $1/5$ | $23/75$ | $1/5$ | $1/5$ |
| $X_1X_2$ | $9/25$ | $89/225$ | $97/225$ | $7/15$ | $11/75$ | $41/225$ | $-1/15$ |
| $X_1X_6$ | $29/125$ | $13/125$ | $1/225$ | $-1/15$ | $23/375$ | $-11/1125$ | $1/25$ |
| $X_1Y_2X_3$ | $1/5$ | $1/5$ | $1/5$ | $1/5$ | $1/5$ | $1/5$ | $-1/3$ |
| $X_1X_2X_6$ | $17/125$ | $121/1125$ | $47/675$ | $1/45$ | $1/125$ | $-1/675$ | $-1/75$ |
| $X_1X_6X_{11}$ | $53/625$ | $-1/75$ | $-11/225$ | $-1/27$ | $1/125$ | $-59/5625$ | $1/125$ |

Each row contracted with the seven masses in (T1) is exactly one. The profiles cover the entire protected error set, not only the listed witnesses:

| Response profile | Number of physical errors through weight three |
|---|---:|
| Identity | 1 |
| Single-error profile | 375 |
| Two-in-one-block profile | 1,350 |
| Two-in-distinct-blocks profile | 2,250 |
| Inner logical triple profile | 150 |
| Two errors in one block and one in another | 27,000 |
| One error in each of three blocks | 33,750 |
| **Total** | **64,876** |

For clarity, the single-error profile includes all 75 weight-one errors and 300 weight-three errors. The two-in-one-block profile includes 450 weight-two errors and 900 weight-three errors. Of the 1,350 triples supported in a single inner block, 150 are nontrivial inner logical words, 300 have the single-error profile, and 900 have the two-in-one-block profile. The finite verifier checks these classifications directly for every physical label and every logical axis.

It follows that

```math
H_\alpha EV=EV\alpha
\qquad(\mathrm{wt}(E)\le3).
\tag{T2}
```

All expanded words lie in one logical coset, whose signed representatives commute. Their maximum physical weight is fifteen. Since the class supports are disjoint, their absolute coefficient sums add:

```math
J/\omega=
\frac{715}{16}+\frac{1755}{16}
+\frac{3645}{32}+\frac{1485}{32}
+\frac{675}{32}+\frac{675}{32}+\frac32
=\frac{2867}{8}.
```

The model requires the transparent Hamiltonian throughout the gate. Commutativity alone does not justify replacing it by a sequence of individually applied terms during which faults can occur; see the [control account](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/control_robustness.md).

### Exact native coefficients

The normalized class averages are uniform, so the relative coefficient of each positive signed-coset atom is simple:

| Class | Relative coefficient per atom |
|---|---:|
| $A_0$ | $143/32000$ |
| $A_1$ | $-39/6400$ |
| $A_2$ | $27/2560$ |
| $A_3$ | $-11/512$ |
| $F_0$ | $-9/128000$ |
| $F_1$ | $1/25600$ |
| $B_5$ | $-1/400000$ |

The bare Hermitian Pauli coefficient is this number multiplied by the atom's logical-coset sign $s_Q$. The common denominator is 3,200,000, whose binary length is 22 bits. The construction trades a large term count for these small exact rational descriptions; no term-count optimum or experimental precision claim is made.

## 4. Why the seven patterns exhaust arbitrary order-fifteen Hamiltonians

The coupling lower bound does not assume the uniform class ansatz in (T1). Expand an arbitrary Hermitian $H$ in the bare tensor-Pauli basis. Only nontrivial logical normalizer terms can contribute to its traceless logical compression. Detectable Paulis have zero compression, and stabilizer Paulis contribute only scalars.

A normalizer word of the concatenated code must be an inner normalizer on each of the five bottom blocks. An inner stabilizer has weight zero or four. An inner nontrivial logical representative has weight three or five. Replacing each inner word by its logical label gives an outer normalizer word. If the whole physical word is a nontrivial logical operator, its outer label is nontrivial and therefore has weight three or five.

If the outer logical has weight three, its three active children start at total weight nine. Promotions cost two each, and padding an inactive child costs four. At total weight at most fifteen, the possibilities are precisely zero through three promotions with no padding, or zero or one promotion with one padded child. Two padded children would already cost seventeen. If the outer logical has weight five, its minimum physical weight is fifteen, forcing every child to have weight three. These are exactly the seven patterns in Section 2.

The argument applies to every logical axis and every native Pauli coefficient, without any symmetry assumption on $H$.

## 5. Block moments and the exact coupling lower bound

For physical site $i$, let

```math
\mathcal R_i(O)=\frac14(O+X_iOX_i+Y_iOY_i+Z_iOZ_i),
\qquad
\mathcal W_b=\sum_{i\text{ in bottom block }b}(\mathrm{Id}-\mathcal R_i).
```

The commuting operators $\mathcal W_b$ multiply a tensor-Pauli word by its weight $w_b$ in block $b$. Products of at most three such factors expand into linear combinations of conjugations by physical Paulis of support at most three. The coefficient sum of each positive-degree monomial is zero. Three-error transparency therefore makes their traceless logical compressions vanish, even when the error sectors carry consistent scalar shifts.

In particular, define

```math
W=\sum_b w_b,\qquad R=\sum_b w_b^2,\qquad S=\sum_b w_b^3.
```

The six moment rows $W,W^2,W^3,R,WR,S$ all have zero traceless logical compression. Let $C_j$ be the total traceless logical contribution from pattern $j$, in the order $(A_0,A_1,A_2,A_3,F_0,F_1,B_5)$. Normalization and these moments give

```math
\begin{pmatrix}
1&1&1&1&1&1&1\\
9&11&13&15&13&15&15\\
81&121&169&225&169&225&225\\
729&1331&2197&3375&2197&3375&3375\\
27&43&59&75&43&59&45\\
243&473&767&1125&559&885&675\\
81&179&277&375&145&243&135
\end{pmatrix}
\begin{pmatrix}
C_{A_0}\\C_{A_1}\\C_{A_2}\\C_{A_3}\\C_{F_0}\\C_{F_1}\\C_{B_5}
\end{pmatrix}
=
\begin{pmatrix}A\\0\\0\\0\\0\\0\\0\end{pmatrix}.
\tag{T3}
```

This integer matrix has determinant 2,949,120. Its unique solution is

```math
(C_{A_0},C_{A_1},C_{A_2},C_{A_3},C_{F_0},C_{F_1},C_{B_5})
=\left(\frac{715}{16},-\frac{1755}{16},\frac{3645}{32},
-\frac{1485}{32},-\frac{675}{32},\frac{675}{32},-\frac32\right)A.
```

Let $J_j$ be the absolute bare coefficient sum of all nontrivial logical words in pattern $j$. Each compressed nontrivial logical Pauli has norm one, so $J_j\ge\lVert C_j\rVert_\infty$. Other native terms only add nonnegative budget. Hence every allowed $H$ satisfies

```math
\boxed{J\ge\sum_j\|C_j\|_\infty
=\frac{2867}{8}\,\omega.}
\tag{T4}
```

An exact scalar dual row, in the same row order as (T3), is

```math
y=\left(\frac{2867}{8},-\frac{2375}{24},7,-\frac16,
\frac{47}{8},-\frac18,-\frac12\right).
```

Its values on the seven columns are $(+1,-1,+1,-1,-1,+1,-1)$, and its right-hand-side value is $2867/8$. The construction (T1) attains this bound, so there is no gap between a restricted uniform ansatz and the arbitrary-Hamiltonian optimum.

If only total weights $9,11,13,15$ are retained, the unique scalar shell masses are $(715,-1755,1485,-429)/16$, giving the weaker relaxation $J/\omega\ge274$. The block refinement increases the lower bound by $675/8$. Different signs within weights thirteen and fifteen are invisible to the total-weight relaxation.

### Exact budget for arbitrary logical directions

The native Pauli basis fixes a canonical logical Pauli frame: the three nontrivial normalizer cosets act as $X,Y,Z$, represented here by $X^{\otimes25},Y^{\otimes25},Z^{\otimes25}$. Write $A=\sum_\beta\Omega_\beta\beta$. Further partition each pattern's terms by these three logical cosets. Every native Pauli in coset $\beta$ compresses to $+\beta$ or $-\beta$. Triangle inequalities on its three logical coefficients therefore give the stronger estimate

```math
J_j\ge\sum_{\beta=X,Y,Z}|[C_j]_\beta|
=|a_j|\sum_{\beta=X,Y,Z}|\Omega_\beta|,
```

where $a_j$ is the corresponding forced mass from (T3). Stabilizer, detectable and scalar terms cannot reduce these absolute coefficient sums. Summing over the seven disjoint patterns gives

```math
J\ge\frac{2867}{8}\sum_\beta|\Omega_\beta|.
```

Conversely, $H=\sum_\beta\Omega_\beta H_\beta+cI$ is strictly transparent by (T2). The three axes use disjoint nontrivial logical cosets, so its native absolute budget is exactly the right-hand side. This proves the stated all-direction optimum at order fifteen. Its basis dependence is physical: the budget is defined in the prescribed native tensor-Pauli basis, whose normalizer cosets select this logical frame. In terms of half-splitting $\omega=\sqrt{\Omega_X^2+\Omega_Y^2+\Omega_Z^2}$, the ratio is $J_{\min}/\omega=(2867/8)\lVert \Omega\rVert_1/\lVert \Omega\rVert_2$. The axis value in (T4) is the case with one nonzero component.

### Why the outer weight-five branch is unavoidable

On the allowed inner normalizer weights $w\in\{0,3,4,5\}$, the polynomial

```math
p(w)=\frac{4w^3-33w^2+68w}{15}
```

is zero on stabilizer weights zero and four, and one on logical weights three and five. Thus $\sum_b p(w_b)$ counts the number of active outer logical factors. It is three on the first six patterns and five on $B_5$.

Its traceless compression vanishes because it is a linear combination of the permitted degree-one, -two and -three block moments. Together with $\sum_j C_j=A$, this gives

```math
3\sum_{j\ne B_5}C_j+5C_{B_5}=0,
\qquad
\boxed{C_{B_5}=-\frac32 A.}
```

Operationally, weight-three errors can be nontrivial logical Paulis of an inner block. Averaging conjugation over the thirty such inner minimum logical words acts trivially on an inner stabilizer and has average character $-1/3$ on any nontrivial inner logical label. Their transparency therefore creates an outer active-count moment. Three-error protection forces a contribution from the outer weight-five branch; omitting it cannot be repaired by changing coefficients among the first six patterns.

## 6. Reproduction and scope

The [exact verifier](../verification/verify_three_error_gates.py) reconstructs the signed five-qubit group and cosets, exhausts the allowed inner weight signatures, and checks the seven normalized class responses for every one of

```math
1+3(25)+9\binom{25}{2}+27\binom{25}{3}=64,876
```

physical Pauli errors, separately for $X,Y,Z$. It also checks the pattern counts, moment matrix, forced masses and exact dual certificate. Tensor-factor multiplication and signed logical replacement establish the action of every class atom. The verifier uses this factorization rather than multiplying a matrix of dimension $2^{25}$ or claiming a separate full native-word expansion for every error. All certificate arithmetic is integer or rational; no numerical optimization is needed.

The code itself has not changed. Its $r_S=\infty$ property still concerns each specified known product unitary followed by full ideal original-syndrome measurement and an informed recovery, as proved for the [underlying concatenated family](minimum_readout_safe_gates.md). The protected-gate guarantee concerns errors through weight three during the declared Hamiltonian evolution with ideal final recovery. These remain distinct noise models; their coexistence does not supply a joint fault-tolerance threshold.

The exact optimum is for this code at interaction order fifteen, with the axis value $2867/8$ and the all-direction native budget above. The [growing-distance extension](growing_three_error_gates.md) attains the corresponding optima at every depth $m\ge2$. The [protected logical-string extension](protected_entangling_gates.md) supplies optimal entanglers on fixed product codes. Geometric locality, physical synthesis and noisy recovery remain further questions.
