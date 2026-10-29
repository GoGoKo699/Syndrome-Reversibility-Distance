# Two-error-transparent gates on the twenty-five-qubit code

The twice-concatenated five-qubit code admits strict transparency to every Pauli error of weight at most two at interaction order thirteen. Its native coupling cost at that order is also exactly determined. For each logical Pauli axis,

```math
\boxed{[[25,1,9]],\qquad t=2,\qquad
\ell=13,\qquad J=\frac{521}{8}\,\omega,\qquad r_S=\infty.}
```

Both interaction order and coupling budget are optimal under the declared fixed-code model. The coupling optimum is specific to this code at interaction order thirteen. It is stronger than the scalar weight-moment relaxation because it distinguishes two physically different kinds of weight-thirteen terms.

The [growing two-error family](growing_two_error_gates.md) extends the optimal resource formula to every concatenation depth $m\ge2$, using a support-tree-correlated padding rule. The present account retains the explicit twenty-five-qubit construction and its sparse certificate.

The code, its distance and readout property, and the canonical five-qubit Hamiltonian are inherited ingredients. The new construction combines four normalized Pauli classes; its optimality proof adds a moment of the interaction weights inside the five inner blocks. The [single-error family](minimum_readout_safe_gates.md) supplies the earlier construction and notation, while the [sparse gate account](sparse_protected_gates.md) treats a separate single-error description problem.

## 1. Model and the first lower bound

Let $V$ encode the twice-concatenated five-qubit code. Regard its physical qubits as five inner blocks of five qubits each. The code has distance nine, encodes one qubit, and has even-weight stabilizers and odd-weight nontrivial logical Pauli representatives. The [readout support criterion](https://github.com/GoGoKo699/Syndrome-Reversibility-Distance/blob/main/science/support_criterion.md) gives $r_S=\infty$.

We require a Hermitian code-preserving Hamiltonian with nonscalar logical generator $h_L$ and full identities

```math
HEV=EV(h_L+\lambda_E I)
\qquad (\operatorname{wt}(E)\le2),
```

where the consistent real scalar sector shifts permitted by the [gate model](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/model.md) are allowed in the lower bound. The construction below has every $\lambda_E=0$. Interaction order $\ell$ is maximum native Pauli weight, and $J$ is the absolute sum of nonidentity Pauli coefficients. Write $\omega$ for half the spectral width of $h_L$.

Only odd interaction weights can contribute nonscalar logical compression. Protecting errors through weight two requires at least three contributing weights by the [available-weight moment argument](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/locality.md#5-available-weights-can-strengthen-the-bound). At distance nine, the first three are $9,11,13$. Hence every nontrivial two-error-transparent generator on this code has

```math
\ell\ge13.
```

## 2. Four normalized signed Pauli classes

For an inner five-qubit block and a logical axis $\beta\in\{X,Y,Z\}$, let $G_{3,\beta}$ and $G_{5,\beta}$ be the normalized sums of its ten weight-three and six weight-five representatives. Each representative retains the Hermitian sign that makes it implement $+\beta$ on the inner code. Thus both normalized sums have logical action $\beta$.

Also define the normalized nonidentity stabilizer sum

```math
K_4=\frac1{15}\sum_{g\in S_5\setminus\{I\}}g
=\frac{16P_5-I}{15},
```

where the signed stabilizers specify the positive code eigenspace and $P_5$ is its projector. All fifteen terms have weight four, and $K_4$ acts as identity on the code.

Fix an outer logical axis $\alpha$. Average over its ten signed weight-three representatives $Q=s_Q Q_1\otimes\cdots\otimes Q_5$. For each such representative, exactly three blocks are active. Replace an active label $Q_j=\beta$ by the corresponding $G_{w,\beta}$ and an inactive label by physical identity, then include the outer sign $s_Q$. Define four normalized sums:

| Class | Active inner factors | Inactive inner factors | Physical weight | Distinct native terms in the dense expansion |
|---|---|---|---:|---:|
| $A_0$ | Three $G_3$ factors | Two identities | 9 | 10,000 |
| $A_1$ | One $G_5$ factor, chosen uniformly among the three active blocks; two $G_3$ factors | Two identities | 11 | 18,000 |
| $A_2$ | Two $G_5$ factors, chosen uniformly among the three active pairs; one $G_3$ factor | Two identities | 13 | 10,800 |
| $B$ | Three $G_3$ factors | One $K_4$, chosen uniformly among the two inactive blocks; one identity | 13 | 300,000 |

Every physical term in each class is a signed positive representative of the same full-code logical axis $\alpha$. Each normalized class therefore obeys $A_jV=V\alpha$ or $BV=V\alpha$. Its absolute native coefficient sum is one.

The inner block logical labels identify the outer representative. Within $A_1$ and $A_2$, the weight-five child blocks identify the promotion choices; within $B$, the nonidentity stabilizer block identifies the padding choice. Thus the printed counts have no duplicate physical words. The four classes are mutually disjoint, including the two weight-thirteen classes, whose inner logical/stabilizer block patterns differ.

## 3. Exact responses to all errors of weight at most two

For a normalized class $T$ and a Pauli error $E$, define its response $f_T(E)$ by the full physical-space identity

```math
TEV=f_T(E)EV\alpha.
```

This form follows term by term because all terms belong to the same signed logical coset. It remains to calculate their Pauli conjugation signs.

### Inner five-qubit responses

For any single-site inner Pauli error $F$,

```math
G_{3,\beta}F V_5=\frac15 FV_5\beta,\qquad
G_{5,\beta}F V_5=-\frac13 FV_5\beta,\qquad
K_4FV_5=-\frac1{15}FV_5.
```

The first two identities follow from the exact signed five-qubit coset counts. The last follows from $P_5FV_5=0$.

If an inner Pauli $E$ has weight two, it also has nonzero syndrome. Since the five-qubit code is perfect, exactly one single-site Pauli $F$ has the same syndrome. The Pauli $F^\dagger E$ is a nontrivial logical operator $\tau$: it normalizes the inner code, has weight at most three, and cannot be a nonidentity stabilizer because those have weight four. Its logical class is one of $X,Y,Z$.

Consequently the two logical-class responses acquire the additional sign $\chi(\tau,\beta)$, equal to $+1$ when the logical Paulis commute and $-1$ when they anticommute:

```math
f_{G_{3,\beta}}(E)=\frac{\chi(\tau,\beta)}5,
\qquad
f_{G_{5,\beta}}(E)=-\frac{\chi(\tau,\beta)}3,
\qquad f_{K_4}(E)=-\frac1{15}.
```

For any nontrivial $\tau$, the average of $\chi(\tau,\beta)$ over $\beta=X,Y,Z$ is $-1/3$.

### Outer averaging and the response table

The ten minimum outer representatives have two elementary properties: their supports run once over all ten three-site subsets, and at each outer site every nonidentity Pauli label occurs twice. These finite five-qubit facts hold for each logical axis. They imply that every physical error of weight at most two has one of the four response rows below:

| Physical error type | Number of errors | $f_{A_0}$ | $f_{A_1}$ | $f_{A_2}$ | $f_B$ |
|---|---:|---:|---:|---:|---:|
| Identity | 1 | $1$ | $1$ | $1$ | $1$ |
| One nonidentity Pauli | 75 | $13/25$ | $31/75$ | $23/75$ | $23/75$ |
| Two Paulis in the same inner block | 450 | $9/25$ | $89/225$ | $97/225$ | $11/75$ |
| Two Paulis in different inner blocks | 2,250 | $29/125$ | $13/125$ | $1/225$ | $23/375$ |

For example, a specified block is inactive with probability $2/5$ in $A_0$, and otherwise its inner factor is $G_3$. Its single-error response is therefore $2/5+(3/5)(1/5)=13/25$. For a weight-two error inside that block, averaging the nontrivial logical class gives $2/5+(3/5)(-1/15)=9/25$.

For errors in different inner blocks, the probability that zero, one, or both of the specified blocks are active is $1/10$, $6/10$, or $3/10$. Single-error inner responses are independent of the inner logical axis, so these support probabilities, together with the uniform promotion or padding choices, determine the last row. For two errors in one block, the uniform outer local-axis counts give the logical-sign average above. Thus the table covers each physical Pauli error, rather than an average over a chosen error distribution.

## 4. Strict transparency and the attained resources

Define

```math
\boxed{
H_\alpha=\frac{143}{8}A_0
-\frac{117}{4}A_1
+\frac{243}{16}A_2
-\frac{45}{16}B.
}
```

The four masses sum to one, and their dot product with every row of the response table is one. Therefore

```math
\boxed{H_\alpha EV=EV\alpha
\quad\text{for every Pauli }E\text{ with }\operatorname{wt}(E)\le2.}
```

These are full physical-space identities, not only compressed energies. They include preservation of each protected error image. Linearity covers arbitrary operators supported on at most two physical qubits. Pointwise transparency also gives the intended logical evolution through up to two physical faults at unknown times, followed by ideal recovery, since every cumulative Pauli error during the trajectory has weight at most two.

The four physical classes are disjoint and their signs are uniform in positive logical-coset coordinates. Hence

```math
\ell=13,\qquad
J(H_\alpha)=\frac{143}{8}+\frac{117}{4}
+\frac{243}{16}+\frac{45}{16}
=\frac{521}{8}.
```

Scaling by $\Omega$ gives logical half-splitting $\omega=|\Omega|$ and $J=(521/8)\omega$. The dense expansion contains 338,800 terms per logical axis. This is a description of that expansion, not a lower bound on the number of terms.

The construction applies separately to $X,Y,Z$. Real linear combinations preserve transparency and allow arbitrary continuous control of the encoded qubit. The budget of the displayed combination is $(521/8)\sum_\alpha|\Omega_\alpha|$; coupling optimality below is attained separately for each Pauli axis, not asserted for every arbitrary logical direction.

### An exact sparse representation

The dense expansion is not needed to specify the gate. The [stored sparse certificate](../verification/data/two_error_sparse_orbits.json) gives 103 signed logical-$Z$ orbit seeds under independent cyclic shifts of the outer block index and the shared within-block index. This is the same transitive 25-element coordinate-translation group used in the [single-error sparse account](sparse_protected_gates.md). Each stored orbit has 25 distinct words, and the 103 orbits are disjoint, giving 2,575 native Pauli terms. Cycling physical labels $X\mapsto Y\mapsto Z\mapsto X$ produces the other two logical axes, whose signs and identities are independently checked.

The sparse relative coefficients have the same four class signs and exact total masses as the displayed construction. Their absolute sum is $521/8$, their interaction order is thirteen, and their conjugation response is one for each of the 2,776 protected errors. A numerical search selected the orbit support; the stored rational coefficients and every full error-image identity are checked independently without an optimizer. The collected native coefficients have a 125-bit common denominator. Thus this representation reduces support at the cost of more complicated coefficients; it is not asserted to minimize term count or denominator size.

There is also a general support guarantee for this fixed problem. Use nonnegative absolute coefficients with the four class signs and one equality for each protected error response, including identity. The dense construction proves feasibility. Exact nonnegative support reduction leaves at most 2,776 independent columns. The block-moment proof below fixes the class masses for every feasible solution, so this reduction retains the optimal budget. The explicit 2,575-term witness improves on that elementary support bound. No arbitrary-depth two-error sparse recursion is established here.

## 5. A code-specific coupling lower bound

The scalar weight moments alone cannot distinguish $A_2$ and $B$, since both have weight thirteen. A second moment of the weights inside the inner blocks supplies the missing constraint.

### Classification of all contributing native Pauli words

A Pauli normalizing the full concatenated code must normalize each inner block, because the inner checks are among the full stabilizers. Its five induced inner logical labels form an outer normalizer. A nontrivial full logical has outer logical weight three or five.

An outer weight-five logical needs at least fifteen physical sites and is excluded at interaction order thirteen. For outer weight three, each active inner factor has weight three or five. Each inactive block is identity or a weight-four inner stabilizer. Total weight at most thirteen therefore leaves precisely four possibilities, up to block permutation:

| Type | Inner-block weights | Total weight $W$ | Sum of squared block weights $R$ |
|---|---|---:|---:|
| $0$ | $(3,3,3,0,0)$ | 9 | 27 |
| $1$ | $(5,3,3,0,0)$ | 11 | 43 |
| $2$ | $(5,5,3,0,0)$ | 13 | 59 |
| $B$ | $(3,3,3,4,0)$ | 13 | 43 |

This classification concerns every nontrivial logical normalizer word, not only terms chosen by the construction. Paulis outside the full normalizer compress to zero, while stabilizers compress to scalars. They cannot contribute to the following scalar-free moments.

### Four operator equations

Let $A=[h_L]_0$ and let $C_0,C_1,C_2,C_B$ be the scalar-free compressed parts of the Hamiltonian terms in the four types. No separate preservation assumption is made for the four sub-Hamiltonians. The normalization and the first two global [weight moments](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/locality.md) give

```math
C_0+C_1+C_2+C_B=A,
```

```math
9C_0+11C_1+13C_2+13C_B=0,
```

```math
81C_0+121C_1+169C_2+169C_B=0.
```

For an inner block $b$, define the Pauli-weight superoperator

```math
\mathcal W_b=\sum_{i\in b}(\mathrm{Id}-\mathcal R_i),
\qquad
\mathcal R_i(O)=\frac14(O+X_iOX_i+Y_iOY_i+Z_iOZ_i).
```

Every native Pauli is an eigenoperator of $\mathcal W_b$ with eigenvalue its weight in that block. Expanding $\mathcal W_b^2$ gives a linear combination of conjugations by errors on at most two physical sites, with coefficients summing to zero. Full two-error transparency, including allowed scalar shifts, thus makes its scalar-free compression vanish. Summing over the five blocks gives the fourth equation

```math
27C_0+43C_1+59C_2+43C_B=0.
```

Pauli conjugation changes only the sign of a native Pauli word, so words outside the normalizer still compress to zero in every moment. Scalar parts and scalar sector shifts disappear from the equations. Thus the argument applies to an arbitrary physical Hamiltonian satisfying the model, even if it contains terms outside the chosen logical coset.

Solving the four linear operator equations gives the unique masses

```math
\boxed{
C_0=\frac{143}{8}A,\quad
C_1=-\frac{117}{4}A,\quad
C_2=\frac{243}{16}A,\quad
C_B=-\frac{45}{16}A.
}
```

### From operator masses to coupling budget

For this one-logical-qubit code, the traceless generator $A$ has operator norm $\omega$. If $J_j$ is the sum of absolute native coefficients in type $j$, then $\|C_j\|_\infty\le J_j$: every contributing compressed Hermitian logical Pauli has norm one, and the triangle inequality applies. Terms in the four types are disjoint. Any additional stabilizer or detectable terms can only add to the native budget. Consequently

```math
J\ge\sum_jJ_j
\ge\sum_j\|C_j\|_\infty
=\left(\frac{143}{8}+\frac{117}{4}
+\frac{243}{16}+\frac{45}{16}\right)\omega
=\frac{521}{8}\omega.
```

The displayed construction attains equality. This proves exact code-specific coupling optimality at interaction order thirteen. If only total weights $9,11,13$ were used in the scalar moment relaxation, its bound would be $J/\omega\ge119/2$. The stronger value $521/8$ comes from the necessarily negative weight-thirteen stabilizer-padding contribution, which the total-weight relaxation merges with the positive two-promotion contribution.

An exact scalar dual certificate makes the last comparison directly reproducible. For the four constraint rows (normalization, $W$, $W^2$, $R$), choose

```math
y=\left(\frac{521}{8},-12,\frac12,\frac18\right).
```

Its values on the four class columns are $(+1,-1,+1,-1)$, precisely the signs of their forced masses. Its value on the right-hand side is $521/8$. The [verification program](../verification/verify_two_error_gates.py) solves the primal and dual systems in exact rational arithmetic.

## 6. Scope and verification

The result establishes exact two-error transparency for this fixed twenty-five-qubit code with native Pauli interactions, and retains its infinite known-product-unitary syndrome reversibility. The code can correct errors through weight four, but this Hamiltonian is proved transparent only through weight two. The resource theorem does not apply to arbitrary distance-nine codes or arbitrary moving-code circuits. No claim is made about geometric locality, efficient many-body synthesis, noisy recovery, or a joint noise threshold.

The [verification guide](../verification/README.md) records the exact certificates. Independent factorized signed-Pauli calculations evaluate all 2,776 errors of weight at most two for each logical axis without constructing a $2^{25}$-dimensional matrix. A separate sparse verifier expands all 2,575 stored native terms per axis, checks signed logical membership, and evaluates every protected error response by integer contraction. Signed logical membership ensures that exact conjugation responses certify full error-image identities. The finite checks support the supplied code and constructions; the all-Hamiltonian lower bound is the operator-moment proof above.
