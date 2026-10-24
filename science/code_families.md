# Equal code parameters with different syndrome reversibility

The ordinary parameters of a stabilizer code do not determine whether its syndrome measurement preserves information after coherent local rotations. This page constructs pairs of codes with identical parameters and unbounded distance. One physical CNOT, applied to the definition of the encoding and stabilizers, changes the syndrome reversibility invariant from infinity to the smallest finite value compatible with the distance.

The operational order is the one in [Model and recovery](model.md): a known product unitary acts, the original full syndrome of the chosen code is measured ideally, and a recovery depending on that unitary acts afterward. The [support criterion](support_criterion.md) and [sharp error exponent](error_exponent.md) supply the two general results used below.

## 1. The family separation

Let $Q_t$ be the $t$-fold concatenated five-qubit code, with $t\geq 1$, and put

```math
N_t=5^t,\qquad m_t=3^t.
```

Use two blocks, $A$ and $B$, with physical indices $0,\ldots,N_t-1$ and $N_t,\ldots,2N_t-1$. Define

```math
\mathcal C_t=Q_t\otimes Q_t,
\qquad
\widetilde{\mathcal C}_t=C_t\mathcal C_t,
\qquad
C_t=\mathrm{CNOT}(0\longrightarrow N_t).
```

Here conjugation by $C_t$ transforms both the code subspace and its full stabilizer group. Each code is subsequently measured in its own prescribed stabilizer syndrome.

**Theorem.** These codes have the same parameters,

```math
\mathcal C_t,\widetilde{\mathcal C}_t:
[[2\cdot5^t,\,2,\,3^t]],
```

but their syndrome reversibility invariants satisfy

```math
r(\mathcal C_t)=\infty,
\qquad
r(\widetilde{\mathcal C}_t)=3^t+1.
```

In particular:

| Concatenation level | Common parameters | Product code $r$ | Code after one CNOT $r$ |
|---|---|---|---|
| $t=1$ | $[[10,2,3]]$ | $\infty$ | $4$ |
| $t=2$ | $[[50,2,9]]$ | $\infty$ | $10$ |
| $t=3$ | $[[250,2,27]]$ | $\infty$ | $28$ |

The product code admits exact recovery after syndrome measurement for every product unitary. For the transformed code, the optimal worst-case weak-rotation error has exponent $3^t+1$. The proof also gives explicit rotations with the exact error function $p_{3^t+1}(q)$ from [Error exponent](error_exponent.md), where $q=\sin^2(\theta/2)$ and $0\leq q\leq 1/2$.

## 2. The five-qubit building block

Use the commuting stabilizer generators

```text
XZZXI
IXZZX
XIXZZ
ZXIXZ
```

with their positive eigenvalues defining the codespace. The canonical logical operators are $\overline X=XXXXX$ and $\overline Z=ZZZZZ$. Pauli strings in algebraic support statements below are understood modulo global phase. Physical syndrome projectors retain the signed stabilizers.

### Distance and explicit minimum logical representatives

The four syndrome bits record anticommutation with the generators in the displayed order. The syndromes of the fifteen single-qubit Pauli errors are:

| Physical qubit | $X$ syndrome | $Y$ syndrome | $Z$ syndrome |
|---|---|---|---|
| 0 | 0001 | 1011 | 1010 |
| 1 | 1000 | 1101 | 0101 |
| 2 | 1100 | 1110 | 0010 |
| 3 | 0110 | 1111 | 1001 |
| 4 | 0011 | 0111 | 0100 |

These are all fifteen distinct nonzero four-bit strings. Consequently, no weight-one Pauli normalizes the code. A weight-two Pauli on distinct sites is the product of two distinct single-site errors. Its syndrome is their binary sum, which is nonzero. There is therefore no normalizer element of weight one or two.

Each nontrivial logical class has a representative of weight three. A useful fixed choice is

```text
logical X: IYYIX
logical Z: ZYIIY
logical Y: YXIIX
```

Each word differs from the corresponding canonical logical Pauli by a stabilizer, up to phase. These statements can be checked by multiplying the four generators, and are verified by complete enumeration in [the exact-parameter check](../verification/verify_code_parameters.py). Thus the code has parameters $[[5,1,3]]$. The chosen logical $X$ avoids physical qubit 0, while the chosen logical $Z$ has a $Z$ factor on that qubit. Both features will be preserved recursively.

### Parity of every stabilizer and logical representative

For phase-free Paulis $a,b$, let $\langle a,b\rangle$ be their binary symplectic pairing. A sitewise check gives

```math
\mathrm{wt}(ab)
\equiv \mathrm{wt}(a)+\mathrm{wt}(b)
+\langle a,b\rangle \pmod 2.
```

At a site where the local factors are different nonidentity Paulis, their product remains nonidentity and they anticommute. The other three cases give the displayed identity immediately. In particular, the parity of weight is additive for commuting Paulis.

All four stabilizer generators have even weight and commute with one another, so every stabilizer has even weight. Each of $XXXXX$, $ZZZZZ$, and $YYYYY$ has odd weight and commutes with the stabilizers. Multiplying any of them by a stabilizer preserves odd parity. Hence **every representative of every nontrivial logical Pauli class has odd weight**. Exact enumeration gives one identity and fifteen weight-four stabilizers, and ten weight-three plus six weight-five representatives in each nontrivial logical class; these counts are useful finite checks, but the parity proof does not rely on them.

## 3. Concatenation preserves the required properties

Concatenate by replacing each physical qubit of the outer five-qubit code with a copy of the preceding inner code. The stabilizer consists of the five inner stabilizer groups together with encoded representatives of the four outer generators. Standard signed logical operators may be used for this definition; changing a representative within its stabilizer coset leaves the same stabilizer group and code subspace after choosing the corresponding sign.

We prove inductively that $Q_t$ has the following properties:

1. Its parameters are $[[5^t,1,3^t]]$.
2. Every stabilizer has even weight, and every nontrivial logical Pauli representative has odd weight.
3. Each of the three nontrivial logical classes has a representative of weight $3^t$.
4. There is such a logical $X$ representative avoiding qubit 0, and such a logical $Z$ representative with local factor $Z$ on qubit 0.

The preceding section proves the base case.

For the distance lower bound, any Pauli normalizing the concatenated code must first normalize each inner code. Replace its action in each inner block by its logical Pauli class. This produces a Pauli on the five outer qubits that normalizes the outer code. If this outer class were a stabilizer, the original physical Pauli would be a product of an encoded outer stabilizer and inner stabilizers, and hence would be a stabilizer of the concatenated code. Thus a nontrivial logical Pauli induces a nontrivial outer logical Pauli. It acts nontrivially on at least three inner logical qubits, each of which requires physical weight at least $3^{t-1}$. Its total weight is at least $3^t$.

To attain that weight in each logical class, take the appropriate weight-three outer representative displayed above and replace each nonidentity factor by a minimum-weight representative of that inner logical class. This has weight $3\cdot3^{t-1}$. It proves the distance upper bound and property 3. The number of encoded qubits is one: the $5(5^{t-1}-1)$ inner checks and four independent outer checks give $5^t-1$ independent stabilizer generators.

For parity, every inner stabilizer remains even. Each outer stabilizer generator acts on four blocks. Replacing each nonidentity factor by an odd-weight inner logical representative gives an even-weight physical representative of that generator. All generators commute, so the parity identity shows that every concatenated stabilizer is even. The canonical logical $X$, $Z$, and $Y$ act on five inner logical blocks; their physical lifts have odd weight. Every representative in these three cosets remains odd after multiplication by a stabilizer.

Finally, recursively lift `IYYIX` to obtain the minimum logical $X$. Its entire first inner block is the identity, so it avoids physical qubit 0. Recursively lift `ZYIIY` to obtain the minimum logical $Z$. Its first inner block is the chosen inner logical $Z$, whose first physical factor is $Z$ by induction. This proves property 4 without requiring an extra local Clifford or an assumption about spatial symmetry.

## 4. The product code is completely reversible

Let $S_t$ and $N_t^{\mathrm P}$ denote the phase-free stabilizer and Pauli centralizer of one block. The superscript distinguishes the centralizer from the physical block length $N_t=5^t$.

A stabilizer of $\mathcal C_t$ is a product $g_Ag_B$, with $g_A,g_B\in S_t$. A nontrivial logical Pauli of the product code is a product $L_AL_B$, with each factor in $N_t^{\mathrm P}$ and at least one factor outside $S_t$.

Suppose their exact physical supports were equal. On a block where the logical factor is nontrivial, its support has odd cardinality, whereas the corresponding stabilizer factor has even support cardinality. Equal supports on that block are impossible. Therefore the two full support sets cannot agree.

The support criterion gives $r(\mathcal C_t)=\infty$. This argument allows logical operators nontrivial on either or both blocks; total weight parity alone would not distinguish all of those cases. The blockwise support comparison is essential.

## 5. One CNOT preserves the distance

Write $m=3^t$ and $N=5^t$. A two-qubit CNOT preserves the number of physical and logical qubits. It changes the weight of any Pauli by at most one because only its two endpoint factors change, and it fixes the identity pair.

Assume that the transformed code has a nontrivial logical Pauli $F$ of weight at most $m-1$. Its inverse conjugate

```math
F_0=C_t^\dagger F C_t
```

is a nontrivial logical Pauli of the product code, and has weight at most $m$. The product code has distance $m$, so its weight must equal $m$. The inverse conjugation must have increased the weight by exactly one. On the endpoint pair this is necessarily a change from one nonidentity factor to two. Thus $F_0$ acts nontrivially in both blocks.

However, any nontrivial logical Pauli of the product code must have a nontrivial logical restriction to at least one block. That restriction alone has weight at least $m$. If the total weight is exactly $m$, the other block must be the identity. This contradicts the conclusion that $F_0$ touches both blocks. The transformed distance is therefore at least $m$.

The minimum logical $X$ in block $A$ constructed above avoids the control qubit. It is also the identity on the target block. The CNOT consequently leaves it unchanged. It remains a nontrivial logical Pauli of weight $m$ in the transformed code, proving the reverse inequality. The distance is exactly $m$.

## 6. The transformed invariant is exactly $m+1$

Let $E_B$ be the chosen minimum logical $Z$ in block $B$. Its weight is $m$ and its local factor on the target qubit $N$ is $Z$. CNOT conjugation gives

```math
E'=C_t E_B C_t^\dagger=Z_0E_B.
```

This is a nontrivial logical Pauli of the transformed code, of even weight $m+1$. We now show that it is a logical circuit in the precise sense of [Support criterion](support_criterion.md): no proper nonempty subword obtained by deleting some of its nonidentity local factors lies in the centralizer.

First, a minimum-weight nontrivial logical $E_B$ has no proper nonempty normalizer subword. If such a subword were itself a nontrivial logical, it would contradict minimum distance. If it were a stabilizer, removing it would give a smaller representative of the same nontrivial logical class, again a contradiction.

Consider any subword $F$ of $E'$. It selects a subset $J$ of the factors of $E_B$, together with a bit $f\in\lbrace 0,1\rbrace $ specifying whether the control factor $Z_0$ is selected. Let $h\in\lbrace 0,1\rbrace $ specify whether the target factor is selected. Since that target factor is $Z$, inverse conjugation gives

```math
C_t^\dagger F C_t=Z_0^{f\oplus h}(E_B)_J.
```

If $F$ normalizes the transformed code, the right-hand side must normalize the product code. A single $Z_0$ does not normalize block $A$: it is not a stabilizer because all stabilizers have even weight, and it is not a nontrivial logical because the distance is at least three. Therefore $f=h$. The remaining factor $(E_B)_J$ must normalize block $B$. By the minimum-weight argument, $J$ is either empty or the entire support of $E_B$. In the first case $f=h=0$ and $F=I$; in the second $f=h=1$ and $F=E'$. These are the only possibilities. Hence $E'$ is an even logical circuit.

The circuit characterization gives $r(\widetilde{\mathcal C}_t)\leq m+1$. The invariant is even and is at least the ordinary distance $m$. Since $m=3^t$ is odd, every finite value of $r$ is at least $m+1$. Combining the bounds proves

```math
r(\widetilde{\mathcal C}_t)=m+1.
```

The axis rotations on the support of $E'$ give an exact witness, rather than only an asymptotic lower bound. Their recovered channel is a logical Pauli dephasing channel with rate $p_{m+1}(\sin^2(\theta/2))$, as proved in [Error exponent](error_exponent.md).

## 7. The ten-qubit example and its physical channel

For $t=1$, the product code consists of two five-qubit blocks. Conjugate it by CNOT $0\to5$. A concrete support collision is

```text
g = XIYXXYXIXZ
L = YIYXXXZIXX
E = ZIIIIZYIIY
```

Here $g$ is a stabilizer, $L$ is a nontrivial logical Pauli, and $gL$ equals $E$ up to phase. The strings $g$ and $L$ have the same support. On each of the four occupied sites of $E$, the factor of $g$ anticommutes with that of $E$.

Apply rotations of angle $\theta$ about $Z_0$, $Z_5$, $Y_6$, and $Y_9$, and leave the other qubits unchanged. For $q=\sin^2(\theta/2)\leq1/2$, the exact optimal entanglement infidelity and exact optimal half-diamond error are

```math
p_4(q)=3q^2-2q^3
=\frac{3}{16}\theta^4+O(\theta^6).
```

The full physical calculation in [verify_physical_channel.py](../verification/verify_physical_channel.py) starts from a $1024\times4$ encoding isometry, applies the rotations, and projects all 256 signed syndrome sectors. It computes the singular values of each physical branch map to obtain the optimal fidelity. It also applies actual Pauli recovery to the projected states and checks that the resulting logical channel is the predicted dephasing channel. Thus the channel construction does not assume the repetition formula it is testing.

At angles $0.1$, $0.2$, $0.4$, $0.8$, and $\pi/2$, the calculated optimal infidelity agrees with $3q^2-2q^3$ within $5\times10^{-16}$ in the recorded run. Exactly eight syndrome sectors have nonzero probability. The maximum Frobenius residual of the recovered normalized Choi matrix from the predicted dephasing channel is below $1.3\times10^{-15}$. A separate, generic product of nonidentical rotations about three-component axes on the original product code gives scalar syndrome effects to a maximum Frobenius residual below $4.8\times10^{-16}$.

The [exact-parameter check](../verification/verify_code_parameters.py) separately enumerates every stabilizer and centralizer element of both ten-qubit codes. It finds $d=3$ for both, no support collisions for the product code, and $r=4$ after CNOT, computed directly from the definition. The [family check](../verification/verify_concatenated_family.py) verifies the recursive generator algebra and the circuit witnesses at $t=1,2,3$. The large-code distance lower bound is the proof above; the script does not claim an exhaustive search at those larger sizes. Commands, dependencies, and saved outputs are described in [Verification](../verification/README.md).

## 8. Interpretation and scope

The CNOT in this construction changes the code and the observables subsequently measured. It is not an additional noise event in the product-unitary model. The comparison identifies a sensitivity of syndrome readout to the physical arrangement of logical information that the conventional parameters $n,k,d$ do not determine.

The two families have unbounded distance and encode two logical qubits. They are not asserted to have a nonvanishing asymptotic rate. The exact recovery allowed here can depend on the full physical unitary and can be computationally or physically expensive. The construction therefore establishes an information-preservation separation under ideal prescribed readout; efficiency, noisy measurement circuits, and code-family threshold constants are separate questions.

The original five-qubit code, ordinary concatenation, and parity mechanisms are established ingredients. Their role and the comparison with coherent-error and syndrome-measurement literature are documented in [Related work](../references/RELATED_WORK.md). The added conclusion here is the explicit equal-parameter separation and its exact syndrome reversibility exponent under the stated operational model.
