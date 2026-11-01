# Optimal protected entanglers on fixed product codes

Two unchanged concatenated five-qubit blocks admit a logical entangling Hamiltonian transparent to every physical Pauli error of weight at most three. At concatenation depth $m\ge2$, put $n_0=5^m$ and $d=3^m$. The product code has parameters $[[2n_0,2,d]]$, and a logical $Z\otimes Z$ generator attains

```math
\boxed{\ell=2d+6,\qquad
J/\omega=J_3(2d),\qquad
J_3(x)=\frac{4x^3+51x^2+170x+24}{24}.}
```

Both interaction order and coupling budget at that order are optimal for the prescribed logical Pauli generator. The first instance is a $[[50,2,9]]$ code with order $24$ and budget $1789$. Its code and original stabilizer measurement remain fixed, and its syndrome reversibility distance remains $r_S=\infty$.

The same construction gives a more general pure logical Pauli-string law. On $M$ product blocks, let $P\in\{I,X,Y,Z\}^{\otimes M}$ have logical weight $p\ge1$, and set $D=pd$. For logical generator $h_L=\Omega P+cI$ with $\Omega\ne0$,

```math
\boxed{\ell_{\min}=D+6,\qquad
J_{\min}(\ell=D+6)=|\Omega|J_3(D).}
\tag{E1}
```

Here $\omega=|\Omega|$. The value $D$ is the minimum physical weight of the specified logical coset. The ordinary distance of the product code is still $d$, since a logical operator can act on just one encoded block. This distinction is necessary for the entangling-gate lower bound.

## 1. Fixed code, target and resource model

Let $V_m$ encode $Q_m=[[5^m,1,3^m]]$, and use the product encoding $V=V_m^{\otimes M}$. The target logical axes are the canonical Pauli frame induced by the transversal $X$, $Y$ and $Z$ representatives on each block. The gate identity is

```math
H_PEV=EVP
\qquad\text{for every physical tensor Pauli }E
\text{ with }\operatorname{wt}(E)\le3.
\tag{E2}
```

This is an equality into the full physical Hilbert space. It preserves every protected error image. Multiplying $H_P$ by $\Omega$ and adding $cI$ gives the stated logical generator. The [native resource model](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/model.md) counts the maximum support of one physical Pauli interaction and the absolute sum of all nonidentity bare Hermitian Pauli coefficients.

The lower bound permits arbitrary Hermitian Hamiltonians satisfying the transparency condition with consistent scalar error-sector shifts. Their Pauli terms may belong to any logical coset, be stabilizers, or be detectable. The proof below projects onto the specified target coefficient, so it does not assume that an arbitrary competing Hamiltonian already has the construction's tensor form.

For this product code, a physical Pauli compresses to a signed logical Pauli string or to zero. A word in the target $P$ coset restricts on each active encoded block to a representative of the corresponding nontrivial logical axis. Its weight there is odd and at least $d$. On an inactive encoded block its restriction is a stabilizer of even weight. Consequently every target representative has weight at least $D=pd$ and parity $p\bmod2$. Tensor products of minimum representatives attain $D$.

## 2. The protected tensor construction

The [single-block hierarchy](growing_three_error_gates.md) supplies strict protected logical Hamiltonians $L_{t,\beta}$ for every $t=0,1,2,3$ and $\beta=X,Y,Z$. Their interaction order is at most $d+2t$. Put $L_{-1,\beta}=0$ and

```math
\Delta_{a,\beta}=L_{a,\beta}-L_{a-1,\beta}.
```

If a physical Pauli $E_i$ in that block has weight smaller than $a$, both terms have identical full action on its error image. Hence

```math
\Delta_{a,\beta}E_iV_m=0
\qquad(a>\operatorname{wt}(E_i)).
\tag{E3}
```

Let $A(P)$ be the active encoded blocks of the target string. Define

```math
H_P=
\sum_{\substack{a_i\ge0\ (i\in A(P))\\\sum_i a_i\le3}}
\bigotimes_{i\in A(P)}\Delta_{a_i,P_i},
\tag{E4}
```

with physical identity factors on inactive encoded blocks. Only increment indices zero through three occur. This finite truncation specifies an exact Hamiltonian.

### Strict error-image identity

Start from the full tensor product of the active $L_{3,P_i}$ factors. For any physical Pauli of total weight at most three, every active block restriction has weight at most three. The full product therefore has the desired strict logical action $P$ on the whole product-code error image, including entangled logical inputs. Inactive encoded blocks have identity factors and their physical errors pass through unchanged.

Expand each active $L_{3,P_i}$ as $\sum_{a=0}^3\Delta_{a,P_i}$. An omitted term from (E4) has $\sum_i a_i>3$. Since $\sum_i\operatorname{wt}(E_i)\le3$, at least one active factor has $a_i>\operatorname{wt}(E_i)$ and vanishes by (E3). Removing every omitted term changes no protected error image. Thus (E4) proves (E2) as a physical-space operator identity.

Each retained tensor term has interaction order at most

```math
\sum_{i\in A(P)}(d+2a_i)
\le pd+6=D+6.
```

The construction uses the complete nested increment hierarchy. The normalized positive components of different protection orders need not have identical higher-error response distributions.

### Native coupling cost

The single-block increment budgets are

```math
b_0(d)=1,\qquad b_1(d)=d,\qquad
b_2(d)=\frac{d(4d+13)}8,\qquad
b_3(d)=\frac{d(4d^2+39d+107)}{24}.
```

Introduce their formal polynomial

```math
B(x,z)=1+xz+\frac{x(4x+13)}8z^2
+\frac{x(4x^2+39x+107)}{24}z^3.
```

Direct coefficient identities give

```math
B(x,z)B(y,z)=B(x+y,z)\pmod{z^4}.
\tag{E5}
```

For example, the degree-two coefficient is $b_2(x)+b_2(y)+xy=b_2(x+y)$, and the degree-three identity adds $x\,b_2(y)+y\,b_2(x)$ to $b_3(x)+b_3(y)$. Therefore $B(d,z)^p=B(pd,z)\pmod{z^4}$.

The triangle inequality and multiplicativity of native Pauli coefficient sums under tensor products now bound the cost of (E4) by

```math
J(H_P)\le\sum_{a=0}^3[z^a]B(d,z)^p
=\sum_{a=0}^3b_a(D)=J_3(D).
\tag{E6}
```

The matching arbitrary-Hamiltonian lower bound below proves equality. Thus the cost proof does not require assuming in advance that tensor contributions have disjoint native support. In the supplied hierarchy their compatible logical-coset signs also ensure that coincident native words do not cancel.

## 3. Minimum order for the target logical coset

For a logical operator $A$, let

```math
\pi_P(A)=2^{-M}\operatorname{Tr}(PA).
```

This selects its coefficient along the target Pauli $P$. Detectable physical Paulis have zero compression; stabilizers and other logical cosets have zero $\pi_P$. A physical Pauli in the target coset has $\pi_P=+1$ or $-1$.

Suppose an arbitrary competing Hamiltonian has target logical coefficient $\Omega\ne0$. Exact transparency, including consistent scalar sector shifts, gives the target-projected weight moments

```math
\sum_w c_w=\Omega,\qquad
\sum_w w^j c_w=0\quad(j=1,2,3).
```

Only weights $D,D+2,D+4,\ldots$ can have nonzero target contributions. If $\ell<D+6$, at most the three weights $D,D+2,D+4$ are available. A polynomial of degree at most three with value one at zero and zero on all available weights then contradicts the four displayed equations. Therefore

```math
\ell\ge D+6.
```

For the 50-qubit $ZZ$ target, this uses the target minimum $D=18$. Applying only the product code's ordinary distance nine would not give this sharp statement.

## 4. Exact coupling lower bound at that order

Partition the physical qubits into bottom five-qubit blocks. Every physical Pauli with nonzero target compression is an inner normalizer on every such block. Its bottom weight is zero or four for an inner stabilizer and three or five for an inner nontrivial logical.

After replacing each bottom word by its logical label, the induced word lies in the target $P$ coset of the product of the depth-$(m-1)$ codes. Its total induced weight is at least $q=pd/3=D/3$ and has the same parity as $q$. Since every nonidentity induced label costs at least three physical qubits, total physical order $D+6$ permits only induced weight $q$ or $q+2$. This reasoning allows arbitrary stabilizer factors on encoded blocks where the target is identity.

At induced weight $q$, weight-five promotions cost two and inactive weight-four stabilizers cost four. At induced weight $q+2$, every active bottom factor must have weight three and no additional factor can occur. Thus exactly the same seven global bottom-weight patterns arise as in the [single-block coupling proof](growing_three_error_gates.md#4-the-arbitrary-hamiltonian-lower-bound), with $d$ there replaced by $D$:

| Pattern | Bottom weights above the minimum target support | $W$ | $R=\sum_b w_b^2$ | $S=\sum_b w_b^3$ |
|---|---|---:|---:|---:|
| $A_j$, $0\le j\le3$ | $j$ weight-five promotions among $q$ active factors | $D+2j$ | $3D+16j$ | $9D+98j$ |
| $F_0$ | One inactive weight-four factor | $D+4$ | $3D+16$ | $9D+64$ |
| $F_1$ | One promotion and one inactive weight-four factor | $D+6$ | $3D+32$ | $9D+162$ |
| $B$ | $q+2$ active factors, all of weight three | $D+6$ | $3D+18$ | $9D+54$ |

The target-projected degree-three block moments are the seven rows

```math
(1,W,W^2,W^3,R,WR,S).
```

Their right-hand side is $(\Omega,0,0,0,0,0,0)$. Products of at most three single-site Pauli twirls give the six zero rows; applying $\pi_P$ removes scalar shifts and every competing logical component. The moment matrix has constant determinant $2,949,120$. Its unique pattern contributions are $m_j(D)\Omega$ with

```math
\begin{aligned}
m_{A_0}&=(D+2)(D+4)(D+6)/48,\\
m_{A_1}&=-D(D+4)(D+6)/16,\\
m_{A_2}&=D(2D+9)(D+6)/32,\\
m_{A_3}&=-D(2D+15)(D+6)/96,\\
m_{F_0}&=-5D(D+6)/32,\qquad
m_{F_1}=5D(D+6)/32,\\
m_B&=-D/6.
\end{aligned}
```

Each target Pauli compresses to a signed $P$, so a pattern's absolute native coefficient sum is at least the absolute value of its target contribution. Every other native term adds nonnegative cost. Consequently

```math
J\ge|\Omega|\sum_j|m_j(D)|
=|\Omega|J_3(D).
\tag{E7}
```

An exact dual in the displayed moment order is

```math
y(D)=\left(J_3(D),
-\frac{12D^2+111D+404}{24},
\frac{D+5}{2},-\frac16,
\frac{D+38}{8},-\frac18,-\frac12\right).
```

It evaluates to $(+1,-1,+1,-1,-1,+1,-1)$ on the seven pattern columns and has objective $J_3(D)$. This lower bound covers arbitrary native terms and consistent scalar sector shifts. Equations (E6) and (E7) establish the exact optimum (E1).

This theorem concerns one specified logical Pauli string. It does not assign a globally optimal weighted sum of these pure-string budgets to an arbitrary logical Hamiltonian at one common interaction-order cap. A lower-weight logical component can use the extra order available under a cap chosen for a higher-weight component.

## 5. The first protected entangler

At $m=2$ and $P=ZZ$, the [seven seed classes](three_error_gates.md) give a compact explicit description of the 50-qubit gate. Write $T_i\otimes T_j$ for the products of their normalized positive signed logical-$Z$ classes. The following table specifies all nonzero class-pair masses. An orientation count of two includes both $T_i\otimes T_j$ and $T_j\otimes T_i$, each with the listed mass.

| Class pair | Mass per orientation | Orientations |
|---|---:|---:|
| $A_0,A_0$ | $220$ | 1 |
| $A_0,A_1$ | $-297$ | 2 |
| $A_1,A_1$ | $243$ | 1 |
| $A_0,A_2$ | $729/4$ | 2 |
| $A_0,A_3$ | $-1485/32$ | 2 |
| $A_1,A_2$ | $-2187/32$ | 2 |
| $A_0,F_0$ | $-135/4$ | 2 |
| $A_0,F_1$ | $675/32$ | 2 |
| $A_1,F_0$ | $405/32$ | 2 |
| $A_0,B_5$ | $-3/2$ | 2 |

There are eighteen ordered class pairs. Their signed masses sum to one and their absolute masses sum to $1789$. The classes are disjoint within each 25-qubit block, so the complete dense specification has 41,032,000,000 distinct native words. Their collected relative per-word coefficients have common denominator 32,000,000,000, of bit length 35. These counts and exact coefficients can be computed from the small class-pair table without enumerating the native list. No native term-count or minimal-denominator optimum is claimed.

For every physical error of total weight at most three, the two blocks restrict it to the already certified child response profiles. Tracking their actual physical weights gives ten complete two-block profiles with multiplicities totaling

```math
1+3(50)+9\binom{50}{2}+27\binom{50}{3}
=540,376.
```

The [exact verifier](../verification/verify_entangling_gates.py) checks the finite error certificate, class-pair coefficients, moment and budget identities, and universal-control statements. Its coverage calculation is separate from the large dense native expansion. The [constructor](../verification/logical_pauli_gates.py) specifies the arbitrary-product gate using the existing exact increment hierarchy.

| Depth $m$ | Active logical blocks $p$ | Physical qubits with $M=p$ | Ordinary distance $d$ | Target minimum $D$ | Optimal order | Optimal axis budget |
|---:|---:|---:|---:|---:|---:|---:|
| $2$ | $2$ | $50$ | $9$ | $18$ | $24$ | $1789$ |
| $2$ | $3$ | $75$ | $9$ | $27$ | $33$ | $40175/8$ |
| $3$ | $2$ | $250$ | $27$ | $54$ | $60$ | $32824$ |

For the pure axis gates and increments, let $q_m$ be the [safe single-block coefficient denominator](growing_three_error_gates.md#5-examples-description-size-and-verification). Formula (E4) has a safe common denominator dividing $q_m^p$. This is a general upper bound; the smaller collected denominator above is an exact property of the finite two-block construction. Coefficient description size, native term count and physical interaction synthesis remain different resources.

## 6. Entanglement and universal encoded control

The code space remains $Q_m^{\otimes M}$ throughout these logical evolutions. For two encoded qubits, the logical unitary at angle $\pi/4$ is entangling:

```math
e^{-i\pi Z_{L,1}Z_{L,2}/4}|++\rangle_L
=\frac{(1-i)|00\rangle+(1+i)|01\rangle
+(1+i)|10\rangle+(1-i)|11\rangle}{\sqrt8}.
```

Either one-qubit reduced state is $I/2$. The six single-logical-qubit Pauli generators together with $Z_L\otimes Z_L$ generate all fifteen nonidentity two-qubit Pauli directions under $i[A,B]$. Their anti-Hermitian multiples therefore generate $\mathfrak{su}(4)$ and give universal continuous control on the two-logical-qubit space. The larger pure-string family supplies every nonidentity logical Pauli direction on an arbitrary fixed product register.

Each generator obeys the same strict error-image identity. Thus any piecewise or time-dependent combination that satisfies that identity at every instant has the corresponding protected time-ordered logical evolution on the fixed code. Multiple faults are covered while every accumulated physical error stays in the protected weight-at-most-three set. These control statements do not assume that large native terms can be applied sequentially during a fault-sensitive pulse; see the [control account](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/control_robustness.md).

## 7. Syndrome readout and operational scope

Every stabilizer of $Q_m$ has even physical weight, and every nontrivial logical coset has odd weight. A product-code stabilizer has an even-weight restriction in every encoded block. A nontrivial product-code logical Pauli has an odd-weight restriction in at least one block. Their physical supports therefore cannot coincide. The [support criterion](https://github.com/GoGoKo699/Syndrome-Reversibility-Distance/blob/main/science/support_criterion.md) gives $r_S=\infty$ for the complete product code, even when a particular target string has even global weight.

The statement holds for entangled logical inputs. The protected logical entangler preserves the code subspace and leaves its stabilizer group and original syndrome measurement unchanged. The readout theorem still concerns each specified known product unitary, followed by full ideal original-syndrome measurement and recovery informed of that unitary. Gate faults and the readout experiment retain their separate operational assumptions.

The results establish exact protected encoded control in the declared native Hamiltonian model. They do not establish a joint noisy gate/readout threshold, efficient decoding, small native expansion, geometric locality, recovery-device fault tolerance, or a physical synthesis advantage for many-body interactions.
