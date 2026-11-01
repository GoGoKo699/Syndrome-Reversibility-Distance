# Three-error optima at every concatenation depth

The depth $m$ concatenated five-qubit code $Q_m$ has parameters $[[5^m,1,3^m]]$. For every $m\ge2$, each logical Pauli axis admits a strict Hamiltonian transparent to every physical Pauli error of weight at most three, with

```math
\boxed{n=5^m,\qquad d=3^m,\qquad
\ell=d+6,\qquad
\frac{J}{\omega}=J_3(d)=\frac{4d^3+51d^2+170d+24}{24}.}
```

Both interaction order and coupling budget at that order are optimal for this code family. The code retains $r_S=\infty$ under its separate known-product-unitary, full ideal original-syndrome readout model. This is an equality statement in the native Pauli resource model, not a combined noisy gate/readout threshold.

The construction starts from the [twenty-five-qubit three-error gate](three_error_gates.md). A truncated tensor expansion lifts the complete hierarchy of zero-, one-, two- and three-error gates together. Its increments have compatible native coefficient signs, so the lifting attains the stronger block-moment coupling bound without cancellation. It does not require extending the earlier uniform-promotion/frontier-padding ansatz.

## 1. A strict transparency lifting theorem

Suppose a child code with distance $D$ has strict logical Pauli Hamiltonians $L_{s,\beta}$ for $s=0,\ldots,T$ and $\beta=X,Y,Z$, satisfying

```math
L_{s,\beta}EV=EV\beta
\quad\text{whenever }\mathrm{wt}(E)\le s,
\qquad
\ell(L_{s,\beta})\le D+2s.
```

Put $L_{-1,\beta}=0$ and define increments

```math
\Delta_{a,\beta}=L_{a,\beta}-L_{a-1,\beta}.
```

For $a>\mathrm{wt}(E)$, both Hamiltonians have the same full action on $EV$, and therefore

```math
\Delta_{a,\beta}EV=0.
\qquad\text{(T1)}
```

Encode five child logical qubits into the five-qubit code. Let $\sigma_RR$ run over its ten positive signed weight-three representatives of the desired outer logical axis. If $I(R)$ is the three-element active child set and $\beta_i$ is the Pauli at child $i$, define

```math
L'_{t,\alpha}
=\frac1{10}\sum_{R}\sigma_R
\sum_{\substack{a_i\ge0\ (i\in I(R))\\\sum_i a_i\le t}}
\bigotimes_{i\in I(R)}\Delta_{a_i,\beta_i},
\qquad 0\le t\le T,
\qquad\text{(T2)}
```

with identity operators on the two inactive children. The child indices in (T2) need not exceed $t$. Equivalently, with a formal variable $z$, the parent increment is the coefficient of $z^t$ in the product of the three child increment polynomials:

```math
\Delta'_{t,\alpha}
=\frac1{10}\sum_R\sigma_R\,[z^t]
\bigotimes_{i\in I(R)}
\left(\sum_{a=0}^{T}\Delta_{a,\beta_i}z^a\right).
\qquad\text{(T3)}
```

These are finite operator identities. The word “truncated” refers to the increment indices, not to an approximation to physical time evolution.

### Full error-image proof

Let a physical Pauli $E$ have total weight at most $t$, and write $E_i$ for its restriction to child $i$, with $s_i=\mathrm{wt}(E_i)$. For a fixed root representative, first consider the full tensor product of the three $L_{t,\beta_i}$ factors. Each child factor acts strictly as $\beta_i$ on its error image. Errors in inactive children commute with their identity factors. The full product therefore implements the root logical representative on the complete parent error image.

Expand that product using $L_{t,\beta_i}=\sum_{a=0}^t\Delta_{a,\beta_i}$. Every omitted term in (T2) has $\sum_i a_i>t$. It cannot have $a_i\le s_i$ for all active children, since $\sum_i s_i\le t$. Thus at least one factor has $a_i>s_i$ and vanishes on $E_iV$ by (T1). Dropping the omitted terms changes no protected error image. Averaging the positive signed outer representatives then gives

```math
L'_{t,\alpha}EV'=EV'\alpha.
```

This proves the full physical-space intertwiner, including preservation of the error image. It is not merely a compressed energy equality. The maximum term order is at most

```math
\sum_{i\in I(R)}(D+2a_i)\le3D+2t.
```

Consequently the lifting theorem gives strict gates through order $T$ on every larger concatenation depth once one seed hierarchy is available. No coefficient sign or coupling-cost assumption is needed for this feasibility statement. Those extra properties are used below to prove native coupling optimality for $T=3$.

## 2. The seed increments and their signs

At depth two let $L_{0,\alpha}=A_0$, the positive normalized minimum-weight logical class. Use the known strict one- and two-error gates as $L_1,L_2$, and the [seven-class three-error construction](three_error_gates.md) as $L_3$. All their classes use the same positive signed-coset convention. The following table gives the coefficients of each normalized class in the increments. A blank entry is zero.

| Class | $\Delta_0$ | $\Delta_1$ | $\Delta_2$ | $\Delta_3$ |
|---|---:|---:|---:|---:|
| $A_0$ | $1$ | $9/2$ | $99/8$ | $429/16$ |
| $A_1$ | | $-9/2$ | $-99/4$ | $-1287/16$ |
| $A_2$ | | | $243/16$ | $3159/32$ |
| $A_3$ | | | | $-1485/32$ |
| $F_0$ | | | $-45/16$ | $-585/32$ |
| $F_1$ | | | | $675/32$ |
| $B$ | | | | $-3/2$ |

Here $A_j$ has $j$ weight-five active bottom factors and all other active factors of weight three. The $F_0,F_1$ classes have one inactive weight-four stabilizer factor and zero or one weight-five promotion. The $B$ class has five active weight-three factors at the seed. Each class atom acts as the positive prescribed logical axis. Its actual bare Hermitian Pauli coefficient is its listed relative coefficient multiplied by its logical-coset sign.

All nonzero entries in one row have the same sign. Moreover, the minimum increment index required to contain an atom is at least its native weight excess over the minimum, divided by two: $A_j$ needs index at least $j$, $F_0$ needs at least two, and $F_1,B$ need index three.

### Sign preservation under lifting

For total increment index at most three, tensor products of child increments have only the following bottom-block patterns:

- Promotions without padding have class sign $(-1)^j$, where $j$ is the total number of weight-five factors.
- A single weight-four factor with no promotion has negative sign. Two such factors would require total index at least four and cannot occur.
- A weight-four factor plus one promotion has positive sign: it is either a child $F_1$ atom or a negative child $F_0$ atom times a negative $A_1$ atom in another child.
- A $B$ atom requires index three; the other child factors then have index zero and positive sign, so its sign remains negative.

These statements apply both to the parent increments (T3) and to their cumulative gates (T2). Every term belongs to the same positive logical coset after including its outer and child coset signs. If different expansions produce the same bare word, that word has one fixed coset sign and the same relative coefficient polarity in every contribution. The contributions add without cancellation.

Induction from the seed therefore preserves the seven class signs

```math
(+,-,+,-,-,+,-)
\quad\text{for }(A_0,A_1,A_2,A_3,F_0,F_1,B)
\qquad\text{(T4)}
```

at every depth. The construction keeps a separate signed increment hierarchy. It does not assume that normalized class distributions in the cumulative gates are identical for different protection orders.

## 3. Native coupling cost of the recursion

Let $K_a(D)$ be the absolute sum of native relative coefficients in the child increment $\Delta_a$ for one logical axis. Sign preservation makes the root average and tensor products in (T3) exact for this cost: collecting coincident Pauli words does not reduce the sum. Thus

```math
K'_0=1,\qquad K'_1=3K_1,\qquad
K'_2=3K_2+3K_1^2,\qquad
K'_3=3K_3+6K_1K_2+K_1^3.
\qquad\text{(T5)}
```

The seed table and exact polynomial substitution give

```math
K_0(D)=1,\qquad K_1(D)=D,\qquad
K_2(D)=\frac{D(4D+13)}8,\qquad
K_3(D)=\frac{D(4D^2+39D+107)}{24}.
\qquad\text{(T6)}
```

Applying (T5) to these four expressions returns $K_a(3D)$ for each $a$. Consequently the cumulative three-error gate has

```math
\frac{J}{\omega}=\sum_{a=0}^3K_a(d)
=\frac{4d^3+51d^2+170d+24}{24}.
```

The lower cumulative orders recover $1$, $d+1$, and $(4d^2+21d+8)/8$. The construction supplies another optimal realization of the earlier two-error family. Its promotion and padding distributions differ from the uniform-promotion/frontier-padding realization in the [earlier account](growing_two_error_gates.md). In particular, the present lifting needs no padding in an inactive root child. The earlier padding correlation condition was a condition on that particular ensemble, not a uniqueness theorem for all optimal Hamiltonians.

## 4. The arbitrary-Hamiltonian lower bound

The logical cosets of $Q_m$ contain only odd-weight representatives. The available-weight moment obstruction therefore gives $\ell\ge d+2t$ for any nontrivial Hamiltonian transparent to all Pauli errors through weight $t$, when $2t+1\le d$. At $t=3$ and $m\ge2$, this proves $\ell\ge d+6$.

For the coupling refinement, divide the code into $N=5^{m-1}$ bottom five-qubit blocks and set $k=3^{m-1}=d/3$. A physical Pauli with nonzero nontrivial code compression must commute with every inner stabilizer. Each bottom factor is consequently an inner stabilizer of weight zero or four, or a nontrivial inner logical of weight three or five. The induced outer word is a nontrivial logical representative of $Q_{m-1}$, of odd weight at least $k$.

If that induced weight exceeds $k+2$, the physical weight is at least $3(k+4)>d+6$. If it is $k+2$, every active bottom factor must have weight three and no other nonidentity factors can occur. If it is $k$, the budget permits at most three weight-five promotions, or one weight-four inactive factor and at most one promotion. Thus exactly seven contributing bottom-weight patterns are possible.

For a pattern write $W=\sum_b w_b$, $R=\sum_b w_b^2$, and $S=\sum_b w_b^3$:

| Class | Bottom-block weights | $W$ | $R$ | $S$ |
|---|---|---:|---:|---:|
| $A_j$, $0\le j\le3$ | $k-j$ blocks of weight three and $j$ of weight five | $d+2j$ | $3d+16j$ | $9d+98j$ |
| $F_0$ | $k$ blocks of weight three and one of weight four | $d+4$ | $3d+16$ | $9d+64$ |
| $F_1$ | $k-1$ blocks of weight three, one of weight five and one of weight four | $d+6$ | $3d+32$ | $9d+162$ |
| $B$ | $k+2$ blocks of weight three | $d+6$ | $3d+18$ | $9d+54$ |

All remaining bottom factors are identity. This classification uses normalizer constraints, not an assumption that an arbitrary Hamiltonian is a uniform class average.

### Degree-three block moments

For a physical site $i$, let

```math
\mathcal R_i(O)=\tfrac14(O+X_iOX_i+Y_iOY_i+Z_iOZ_i),
\qquad
\mathcal W_b=\sum_{i\in b}(\mathrm{Id}-\mathcal R_i).
```

The operator $\mathcal W_b$ multiplies a Pauli term by its weight $w_b$ in block $b$. A polynomial in these block twirls of degree at most three and with zero constant term expands into conjugations by physical Pauli errors of weight at most three, with coefficients summing to zero. Full error transparency makes every such conjugated compression equal to the same logical generator modulo scalars. Its scalar-free polynomial moment therefore vanishes.

This gives the six zero moments represented by

```math
W,\quad W^2,\quad W^3,\quad R,\quad WR,\quad S.
```

Let $A=[V^\dagger HV]_0$ and let $C_j$ be the scalar-free contribution from physical Pauli terms in pattern $j$. The seven rows $(1,W,W^2,W^3,R,WR,S)$ give an exact operator linear system with right-hand side $(A,0,0,0,0,0,0)$. Its determinant is the nonzero constant $2,949,120$. To see why it is independent of $d$, put $a=W-d$, $b=R-3d$ and $c=S-9d$. Unit-triangular row operations replace the seven rows by $(1,a,a^2,a^3,b,ab,c)$, whose columns are constant; in particular $ab=WR-dR-3dW+3d^2$. The seed determinant then gives the displayed value. The unique solution is $C_j=m_j(d)A$, where

| Class | Forced signed mass $m_j(d)$ |
|---|---:|
| $A_0$ | $(d+2)(d+4)(d+6)/48$ |
| $A_1$ | $-d(d+4)(d+6)/16$ |
| $A_2$ | $d(2d+9)(d+6)/32$ |
| $A_3$ | $-d(2d+15)(d+6)/96$ |
| $F_0$ | $-5d(d+6)/32$ |
| $F_1$ | $5d(d+6)/32$ |
| $B$ | $-d/6$ |

Detectable Pauli terms have zero code compression, and stabilizer terms contribute only scalars. They do not invalidate the system and add nonnegative native budget. Consistent scalar shifts between protected error sectors likewise disappear on taking the scalar-free part.

For one logical Pauli axis, the absolute native coefficient sum in a class is at least $\lVert C_j\rVert_\infty=|m_j|\omega$. Hence every allowed Hamiltonian satisfies

```math
J\ge\sum_j|m_j(d)|\,\omega=J_3(d)\omega.
```

The recursion has precisely the signs in (T4), so its class masses attain this lower bound. An exact dual row, in the displayed moment order, is

```math
y(d)=\left(
J_3(d),
-\frac{12d^2+111d+404}{24},
\frac{d+5}{2},
-\frac16,
\frac{d+38}{8},
-\frac18,
-\frac12
\right).
```

It takes the values $(+1,-1,+1,-1,-1,+1,-1)$ on the seven pattern columns and has objective $J_3(d)$. The [exact verifier](../verification/verify_growing_three_error.py) checks both the primal and dual identities coefficient by coefficient in $d$, rather than inferring them from finitely many numerical samples.

Keeping only the four total weights $d,d+2,d+4,d+6$ gives the weaker scalar relaxation

```math
\frac{J_{\rm scalar}}{\omega}
=\frac{d^3+9d^2+20d+6}{6},\qquad
J_3(d)-\frac{J_{\rm scalar}}{\omega}
=\frac{5d(d+6)}8.
```

The scalar relaxation combines opposite-sign contributions at weights $d+4$ and $d+6$. Resolving their bottom-block patterns exposes the additional absolute coefficient cost. At $d=9$, the scalar value is $274$ and the exact code-specific value is $2867/8$.

### A specified logical direction in the canonical Pauli basis

There is also a sharper direction-resolved statement. Fix the canonical logical axes represented by the transversal $X$, $Y$ and $Z$ operators, and write

```math
A=\Omega_X X+\Omega_Y Y+\Omega_Z Z.
```

Every physical Pauli with nontrivial code compression acts as one of these three logical Paulis, up to its signed-coset sign. Thus the absolute coefficient sum in pattern $j$ is at least $|m_j|\sum_\alpha|\Omega_\alpha|$, not merely its spectral-norm lower bound. Summing the three constructed axes attains

```math
\boxed{J_{\min}(\Omega_X,\Omega_Y,\Omega_Z)
=J_3(d)\left(|\Omega_X|+|\Omega_Y|+|\Omega_Z|\right).}
```

The logical half-splitting is $\omega=\sqrt{\Omega_X^2+\Omega_Y^2+\Omega_Z^2}$. The ratio to $\omega$ therefore depends on direction. A change of logical coordinate basis is not a free physical rotation of the native Pauli interactions; the specified canonical basis is part of this direction-resolved statement.

## 5. Examples, description size and verification

| Depth $m$ | Physical qubits $n$ | Distance $d$ | Optimal order $\ell$ | Optimal axis budget $J/\omega$ |
|---:|---:|---:|---:|---:|
| $2$ | $25$ | $9$ | $15$ | $2867/8$ |
| $3$ | $125$ | $27$ | $33$ | $40175/8$ |
| $4$ | $625$ | $81$ | $87$ | $824723/8$ |

The [recursive implementation](../verification/three_error_recursion.py) records the short increment rule. It does not assert that a dense native expansion is small. At the seed, the seven positive classes already contain 1,480,960 distinct terms per axis; tensor expansion grows rapidly with depth. The recursion is an exact description of the Hamiltonian, separate from the cost of expanding or physically synthesizing its individual many-body terms.

The direct recursion nevertheless has a controlled exact coefficient size. The seed gates and increments have a common native coefficient denominator dividing $q_2=3,200,000$. Three child products and the ten-term outer average give the safe recurrence

```math
q_m=10q_{m-1}^3
=q_2^{\,3^{m-2}}10^{(3^{m-2}-1)/2}.
```

Every collected native coefficient has denominator dividing $q_m$; summing repeated words does not increase this denominator. Its bit length is 22, 69, 208 and 627 at depths two through five, respectively, and is $O(d)$ in general. These are safe upper bounds, not minimal common denominators or hardware precision requirements. Since the absolute coefficient sum is polynomial in $d$, the corresponding numerator bit length is also $O(d)$. The [implementation](../verification/three_error_recursion.py), through `coefficient_denominator_bound(level)`, computes this bound without expanding the native list.

A separate existence bound follows from response-preserving positive support reduction. With

```math
D_3(n)=1+3n+9\binom n2+27\binom n3,
```

the full transparency equations for an axis have $D_3(n)$ character coordinates. Use the prescribed class signs and nonnegative absolute coefficients. The construction is feasible; eliminating dependent response columns leaves at most $D_3(n)=O(n^3)$ native Pauli terms. The same block moments force the same class masses and therefore preserve the optimum. Integer sign columns give polynomial-bit rational coefficients for an independent support by Cramer's rule and the Hadamard bound. This is a support-existence statement. A practical native sparse compiler or explicit sparse 125-qubit three-error word list is not supplied here.

The [growing verifier](../verification/verify_growing_three_error.py) checks the exact polynomial moments and dual, the increment cost recurrence, the forced pattern masses at finite depths, and the signed increment base. At depth three, eleven complete child-error response profiles cover all 8,649,376 physical Pauli errors of weight at most three on 125 qubits. The computation evaluates these exact profiles and their coverage; it does not multiply a matrix of dimension $2^{125}$ or claim to enumerate a dense native Hamiltonian term by term. The arbitrary-depth physical-space conclusion follows from the lifting proof.

All three axes have the same construction. The gate guarantee remains fixed at three physical Pauli errors even when the code can correct more. The infinite syndrome-reversibility claim still allows recovery to know the specified product unitary and assumes ideal full original-syndrome measurement. State preparation, noisy recovery, geometric range, many-body interaction synthesis and a joint noisy gate/readout threshold remain outside the declared model.
