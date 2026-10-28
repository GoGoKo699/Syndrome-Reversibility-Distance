# Minimum-locality protected gates at growing distance

The iterated five-qubit code admits exact single-error-transparent logical gates at the smallest interaction order compatible with infinite syndrome reversibility distance. At concatenation depth $m\ge1$, the construction has

```math
\boxed{n=5^m,\qquad d=3^m,\qquad
\ell=d+2,\qquad J=(d+1)\omega,\qquad r_S=\infty.}
```

The locality and coupling equalities hold for each of the three logical Pauli axes. Linear combinations give arbitrary continuous control of the one encoded qubit while retaining exact transparency to every single-qubit error. The number of native Pauli terms grows exponentially with $d$; these resource equalities do not establish an economical hardware implementation.

The five-qubit code, its canonical protected Hamiltonians, and the concatenated code family are inherited constructions. The result developed here is the recursive choice of two Pauli-weight classes and its exact gate-resource analysis. It extends the distance-three comparison in the [gate/readout tradeoff](gate_readout_tradeoff.md). Throughout this page, the protected error weight is **$t=1$**, including when the code distance is larger than three.

## 1. Code family and readout property

Use the standard five-qubit stabilizer code with positive generators

```math
XZZXI,\qquad IXZZX,\qquad XIXZZ,\qquad ZXIXZ,
```

and logical representatives $XXXXX$, $YYYYY$, and $ZZZZZ$. Denote its encoding by $V_1$. Define $V_m$ recursively by encoding each of the five outer physical qubits with $V_{m-1}$, so the code at depth $m$ has $n=5^m$ physical qubits and encodes one logical qubit.

Its distance is $d=3^m$. A physical normalizer induces a logical Pauli label on each inner block. For a nontrivial logical of the concatenated code, these five labels form a nontrivial logical of the outer code, so at least three inner blocks have nontrivial logical labels. Induction gives weight at least $3\cdot3^{m-1}$. Choosing a weight-three outer representative and minimum-weight representatives in its three active blocks attains this bound.

Every stabilizer of the concatenated code has even weight, and every nontrivial logical Pauli has odd weight. At the base level, nonidentity stabilizers have weight four and nontrivial logical representatives have weights three or five. Inductively, the inner stabilizers are even, and lifting a weight-four outer stabilizer gives a product of four odd-weight logical representatives on disjoint blocks, again of even weight. Products of commuting Paulis add their weight parities, so the full stabilizer group remains even. A chosen nontrivial logical lift is odd, and multiplication by a commuting even stabilizer preserves that parity.

Consequently no stabilizer and nontrivial logical have equal physical support. The [syndrome support criterion](https://github.com/GoGoKo699/Syndrome-Reversibility-Distance/blob/main/science/support_criterion.md) therefore gives $r_S=\infty$ for every $m$. This means exact recoverability of every specified product unitary after full ideal measurement of the original syndrome, when recovery knows the unitary. It does not include unknown rotations or faults in syndrome extraction. The [code-family account](https://github.com/GoGoKo699/Syndrome-Reversibility-Distance/blob/main/science/code_families.md) establishes the same concatenated-family readout property.

## 2. The five-qubit ingredient

Fix a logical axis $\alpha\in\{X,Y,Z\}$. Let $\mathcal C_{\alpha,w}$ be the **signed Hermitian** representatives of that logical axis having physical weight $w$. Each $Q\in\mathcal C_{\alpha,w}$ is signed so that

```math
QV_1=V_1\alpha.
```

The five-qubit logical coset contains ten representatives of weight three and six of weight five. In either class, at every physical site, each nonidentity local label $X$, $Y$, or $Z$ occurs exactly twice. This count is independent of the logical axis $\alpha$. Thus the two normalized sums

```math
T_{1,\alpha}^{(0)}=\frac1{10}\sum_{Q\in\mathcal C_{\alpha,3}}Q,
\qquad
T_{1,\alpha}^{(1)}=\frac1{6}\sum_{Q\in\mathcal C_{\alpha,5}}Q
```

have local-label probabilities $1/5$ and $1/3$, respectively, in their positive signed-coset coordinates. Both compress to $\alpha$. These probabilities refer to the coefficients of the displayed normalized sums, not to a physical randomization protocol.

At $m=1$, the Hamiltonian constructed below is

```math
H_{1,\alpha}=\frac52T_{1,\alpha}^{(0)}
-\frac32T_{1,\alpha}^{(1)}.
```

Its signed-coset coefficients are $+1/4$ at weight three and $-1/4$ at weight five. It is therefore exactly the [canonical five-qubit protected Hamiltonian](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/coupling_cost.md#the-five-qubit-benchmark), not a new base-code construction.

## 3. Recursive normalized sums

We construct two normalized signed-coset sums $T_{m,\alpha}^{(0)}$ and $T_{m,\alpha}^{(1)}$. Every term in the first has weight $d=3^m$, and every term in the second has weight $d+2$. Both implement $\alpha$ on the code.

For $m\ge2$, write an outer representative as $Q=s_Q Q_1\otimes\cdots\otimes Q_5\in\mathcal C_{\alpha,3}$, where $s_Q\in\{+1,-1\}$ and each $Q_j$ is a literal Pauli label. Define the block operators

```math
B_j^{(0)}(Q)=
\begin{cases}
I^{\otimes5^{m-1}},&Q_j=I,\\
T_{m-1,Q_j}^{(0)},&Q_j\ne I.
\end{cases}
```

The identity acts on an entire block of $5^{m-1}$ qubits. Then set

```math
T_{m,\alpha}^{(0)}
=\frac1{10}\sum_{Q\in\mathcal C_{\alpha,3}}
s_Q\bigotimes_{j=1}^{5} B_j^{(0)}(Q).
```

For the second sum, choose one of the three active child blocks and use its second sum instead:

```math
T_{m,\alpha}^{(1)}
=\frac1{30}\sum_{Q\in\mathcal C_{\alpha,3}}s_Q
\sum_{a\in\operatorname{supp}(Q)}
\bigotimes_{j=1}^{5} B_j^{(a)}(Q),
```

where $B_a^{(a)}(Q)=T_{m-1,Q_a}^{(1)}$ and $B_j^{(a)}(Q)=B_j^{(0)}(Q)$ for $j\ne a$.

Equivalently, the first sum chooses a weight-three logical representative at every active node of the concatenation tree. The second chooses one active bottom block uniformly and replaces its weight-three representative by a weight-five representative. This raises the total weight by exactly two.

Each expansion has positive coefficients in signed-coset coordinates, with coefficient sum one. A sampled signed Pauli acts as the prescribed outer logical axis because every inner representative acts as its own block label. Hence, as full physical-space identities,

```math
T_{m,\alpha}^{(0)}V_m=V_m\alpha,
\qquad
T_{m,\alpha}^{(1)}V_m=V_m\alpha.
```

### One-site label counts

Let $\mu_m^{(a)}$ be the coefficient mass of any specified nonidentity local label at any specified physical site in $T_{m,\alpha}^{(a)}$. It is independent of the site, the local label, and the logical axis. The base counts above give $\mu_1^{(0)}=1/5$ and $\mu_1^{(1)}=1/3$.

A fixed outer child is active with probability $3/5$. When it is active in the second construction, it receives the second sum with probability $1/3$. Therefore

```math
\mu_m^{(0)}=\frac35\mu_{m-1}^{(0)},
\qquad
\mu_m^{(1)}=\frac15\mu_{m-1}^{(1)}
+\frac25\mu_{m-1}^{(0)}.
```

Induction yields the exact isotropy identities

```math
\boxed{
\mu_m^{(0)}=\frac{d}{3n},\qquad
\mu_m^{(1)}=\frac{d+2}{3n}.
}
```

Only these one-site identities are needed for single-error transparency. They do not establish protection against two or more faults.

## 4. Full error transparency and optimal resources

Define the unit-logical-strength Hamiltonian

```math
\boxed{
H_{m,\alpha}
=\frac{d+2}{2}T_{m,\alpha}^{(0)}
-\frac d2T_{m,\alpha}^{(1)}.
}
```

The coefficient sum in signed-coset coordinates is one, so $H_{m,\alpha}V_m=V_m\alpha$. For each site and each nonidentity local label, its signed coefficient mass is

```math
\frac{d+2}{2}\frac{d}{3n}
-\frac d2\frac{d+2}{3n}=0.
```

To check full transparency, let $E$ be any single-site Pauli. Every signed-coset representative $Q$ obeys

```math
QE V_m=\chi(E,Q)E V_m\alpha,
\qquad \chi(E,Q)\in\{+1,-1\}.
```

Conjugation by $E$ reverses precisely the two nonidentity labels anticommuting with it at its site. Both have zero total signed coefficient mass. The total coefficient therefore remains one, proving

```math
\boxed{
H_{m,\alpha}E V_m=E V_m\alpha
\qquad(E=I,X_i,Y_i,Z_i).
}
```

This is an equality into the full physical Hilbert space, including absence of leakage from the protected error images. Linearity covers arbitrary single-qubit error operators. The [gate model](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/model.md) then gives the intended logical evolution for a fault at an unknown time, followed by ideal recovery.

The two sums have disjoint weights $d$ and $d+2$. Within each signed coset, repeated contributions, if any, have the same sign, so each normalized sum has absolute native coefficient sum one. Consequently

```math
\ell=d+2,\qquad
J(H_{m,\alpha})=\frac{d+2}{2}+\frac d2=d+1.
```

Scaling by $\Omega$ gives logical half-splitting $\omega=|\Omega|$ and $J=(d+1)\omega$.

These resources are optimal within the declared model. Since $r_S=\infty$, the [joint locality theorem](gate_readout_tradeoff.md) excludes $\ell=d+1$ for a nontrivial single-error-transparent generator; the general gate bound excludes smaller orders. Thus $d+2$ is the minimum possible order. At that order, the [single-error coupling bound](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/coupling_cost.md#one-protected-error) gives

```math
J\ge\frac{\ell+d}{\ell-d}\,\omega=(d+1)\omega,
```

which the construction attains. This is coupling optimality at the minimum permitted order, not an assertion that larger-order gates cannot use a smaller coupling budget.

## 5. Term count and logical control

Let $a_m$ and $b_m$ count the distinct terms in $T_{m,\alpha}^{(0)}$ and $T_{m,\alpha}^{(1)}$. Their recurrences are

```math
a_1=10,\quad b_1=6,\qquad
a_m=10a_{m-1}^3,\quad
b_m=30b_{m-1}a_{m-1}^2.
```

The outer logical labels distinguish different outer choices. In the second sum, the unique child block with weight $3^{m-1}+2$ identifies the chosen child, and induction distinguishes its inner choices. Thus the counts have no duplicate physical words. Solving the recurrences gives

```math
\boxed{
a_m=10^{(d-1)/2},\qquad
b_m=\frac d5a_m,\qquad
\#\operatorname{terms}(H_{m,\alpha})
=\left(1+\frac d5\right)10^{(d-1)/2}.
}
```

| Depth $m$ | Physical qubits $n$ | Distance $d$ | Interaction order $\ell$ | $J/\omega$ | Terms per logical axis |
|---:|---:|---:|---:|---:|---:|
| 1 | 5 | 3 | 5 | 4 | 16 |
| 2 | 25 | 9 | 11 | 10 | 28,000 |
| 3 | 125 | 27 | 29 | 28 | 64,000,000,000,000 |

The term count is exponential in $d$, even though the total absolute coupling sum is linear in $d$. The construction does not bound geometric range, control complexity, or the cost of synthesizing these native interactions.

All terms for a fixed logical axis commute: they belong to one stabilizer coset. This does not make their sequential application transparent to faults during the sequence; the required cancellations must hold throughout the gate, as explained in [control robustness](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/control_robustness.md).

Each of the three axes has the same construction. Therefore any pointwise combination with real coefficients $\Omega_\alpha(t)$

```math
H(t)=\sum_{\alpha=X,Y,Z}\Omega_\alpha(t)H_{m,\alpha}
```

satisfies the full transparency equations with logical generator $\sum_\alpha\Omega_\alpha(t)\alpha$. This gives arbitrary continuous control of the single encoded qubit. It does not supply entangling gates between encoded qubits or establish multiqubit computational universality.

The three logical cosets are disjoint, so the combined native budget is

```math
J(t)=(d+1)\sum_{\alpha=X,Y,Z}|\Omega_\alpha(t)|,
\qquad
\omega(t)=\sqrt{\Omega_X(t)^2+\Omega_Y(t)^2+\Omega_Z(t)^2}.
```

Thus the coupling optimum stated above is attained separately for each Pauli axis; the displayed combination is not claimed coupling optimal for every arbitrary logical direction.

## 6. Scope and verification

The construction protects one arbitrary single-qubit fault at an unknown time. Although the depth-$m$ code corrects all errors of weight at most $(3^m-1)/2$, this Hamiltonian construction proves only $t=1$ transparency. A higher-error gate must satisfy additional multisite identities. For $m=2$ and $E=X_1X_2$ within the first bottom block (one-based indices), direct signed-term enumeration gives the two class responses $9/25$ and $89/225$. Therefore

```math
H_{2,Z}EV_2=\left(\frac{11}{2}\frac9{25}-\frac92\frac{89}{225}\right)EV_2Z
=\frac15EV_2Z,
```

which exhibits failure of transparency to this correctable weight-two fault.

**Higher-error obstruction for this family.** If a nontrivial Hamiltonian on the depth-$m$ code is transparent to every Pauli error of weight at most $t$, where $1\le t$ and $2t+1\le d$, then

```math
\boxed{\ell\ge d+2t.}
```

Indeed, every nontrivial logical representative on this code has odd weight, so the scalar-free logical components $A_w$ vanish at every even weight, as well as below $d$. If $\ell<d+2t$, only the $t$ weights $d,d+2,\ldots,d+2t-2$ can contribute. The [available-weight moment argument](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/locality.md#5-available-weights-can-strengthen-the-bound) requires at least $t+1$ contributing weights, a contradiction. For example, on the depth-two code, $d=9$ and $t=2$ require $\ell\ge13$; no saturating two-error construction is supplied here. This family-specific argument uses odd logical-coset weights and is not asserted for every code with $r_S=\infty$.

Exact syndrome reversibility and exact gate transparency remain separate operational statements about the same fixed code. The former permits recovery tailored to a known product unitary after ideal full syndrome measurement; the latter assumes the specified native Hamiltonian and ideal final recovery after a protected fault. Their coexistence does not automatically establish a combined noise threshold or a complete implementation advantage.

The [verification guide](../verification/README.md) records the finite certificates. Exact enumeration of the base cosets checks the local-label counts and signs. At depth two the verifier expands all 28,000 terms for each of the three axes, independently checks signed logical-coset membership, and certifies all 76 no-error/single-Pauli intertwiners per axis. It also checks the two-error counterexample and recursive integer/rational identities through depth four. The recurrences and proofs above establish the arbitrary-depth result without enumerating exponentially many terms.
