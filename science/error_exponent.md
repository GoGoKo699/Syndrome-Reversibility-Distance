# The exact weak-rotation error exponent

The support classification has a quantitative consequence. For a fixed stabilizer code, its syndrome reversibility distance $r$ is the exact worst-case exponent of the error that remains after optimal recovery from the original full syndrome measurement. The recovery knows the physical rotation and is unrestricted after that measurement.

The [model](model.md) defines $\mathcal N_U$, the optimal entanglement infidelity $e_{\mathrm{opt}}(U)$, and the optimal half diamond error $\varepsilon_{\mathrm{opt}}(U)$. The [support theorem](support_criterion.md) defines $r$ and proves its equivalent characterization as the minimum even weight of a logical circuit.

## 1. Statement and quantifiers

For Hermitian single-qubit Hamiltonians satisfying $\Vert h_j\Vert _\infty\leq1$, set

```math
U_{\boldsymbol h}(\theta)
=\bigotimes_{j=1}^{n}\exp(-i\theta h_j/2).
\qquad\text{(1)}
```

Identity components of $h_j$ only contribute global phases. Define the worst-case optimized errors

```math
\begin{aligned}
\mathcal E_S(\theta)
&=\sup_{\Vert h_j\Vert _\infty\leq1}
e_{\mathrm{opt}}(U_{\boldsymbol h}(\theta)),\\
\mathcal D_S(\theta)
&=\sup_{\Vert h_j\Vert _\infty\leq1}
\varepsilon_{\mathrm{opt}}(U_{\boldsymbol h}(\theta)).
\end{aligned}
\qquad\text{(2)}
```

The recovery is optimized separately for each Hamiltonian tuple. These are not minimax errors of a common decoder for unknown axes.

**Theorem 3.** If $r(S)=\infty$, both errors in (2) vanish for every $\theta$. If $r(S)<\infty$, then, for the fixed code as $\theta\to0$,

```math
\mathcal E_S(\theta)=\Theta(|\theta|^r),
\qquad
\mathcal D_S(\theta)=\Theta(|\theta|^r).
\qquad\text{(3)}
```

The upper bound is uniform over all Hamiltonians in (1). One fixed tuple of Pauli axes, with $h_j=0$ outside a minimum witness, attains the lower bound. Constants in (3) can depend on the code.

The proof below supplies an explicit finite-code upper bound and a witness with an exact error formula at every $0\leq\theta\leq\pi/2$.

## 2. A finite-code collision polynomial

Let $\mathcal C_S$ contain the pairs $(g,L)$ with $g\in S$, $L\in N\setminus S$, and exactly equal support. Define the polynomial

```math
B_S(x)=
\sum_{(g,L)\in\mathcal C_S}
x^{\mathrm{wt}(gL)}.
\qquad\text{(4)}
```

Each phase-free pair is counted once; physical signs do not change the count. If $r$ is finite, its smallest nonzero power is $r$. All powers are even, because $g$ and $L$ commute globally and differ on an even number of active sites. If there are no collisions, $B_S$ is the zero polynomial.

**Proposition 3.** For every tuple in (1) and $|\theta|\leq1$,

```math
e_{\mathrm{opt}}(U_{\boldsymbol h}(\theta))
\leq\varepsilon_{\mathrm{opt}}(U_{\boldsymbol h}(\theta))
\leq\min\{1,B_S(|\theta|)\}.
\qquad\text{(5)}
```

In particular, if $r<\infty$,

```math
\mathcal D_S(\theta)\leq B_S(1)|\theta|^r
\qquad (|\theta|\leq1).
\qquad\text{(6)}
```

The polynomial is a quantitative certificate, not an asserted efficient way to compute $r$ or implement recovery. Its number of terms can grow exponentially with code size.

### 2.1 Suppression of changed local axes

Let $P_j$ be a nonidentity single-qubit Pauli. Its conjugate by $U_j(\theta)$ is a real combination of the three nonidentity Paulis. The coefficient of $P_j$ has absolute value at most one. For a different axis $Q_j$, its coefficient has absolute value at most $|\theta|$. Indeed,

```math
\left\Vert \frac{d}{d\theta}
\left(U_j(\theta)^\dagger P_jU_j(\theta)\right)\right\Vert _\infty
=\frac12\Vert [h_j,P_j]\Vert _\infty
\leq1.
\qquad\text{(7)}
```

At $\theta=0$ the $Q_j$ coefficient is zero, and its magnitude is bounded by the operator-norm difference from $P_j$.

For a physical stabilizer $g$, a canonical Pauli word $L$ in the expansion of $U^\dagger gU$ has the same exact support as $g$. Its coefficient has absolute value at most

```math
|\theta|^{\mathrm{wt}(gL)},
\qquad\text{(8)}
```

since the differing local axes occur precisely at the nonidentity sites of the phase-free product $gL$.

### 2.2 Bounding the nonscalar syndrome effects

Write

```math
p_s=D^{-1}\mathrm{Tr}M_s,
\qquad
\Delta_s=M_s-p_sI_D.
\qquad\text{(9)}
```

In the Pauli expansion of a compressed conjugated stabilizer, nonnormalizing words vanish, stabilizer words contribute only to the scalar part, and nontrivial logical words are traceless with operator norm one. Combining (8) with the syndrome Fourier factor $1/|S|$ therefore gives, for every $s$,

```math
\Vert \Delta_s\Vert _\infty\leq|S|^{-1}B_S(|\theta|).
\qquad\text{(10)}
```

There are $|S|$ syndrome labels, so

```math
\sum_s\Vert \Delta_s\Vert _\infty\leq B_S(|\theta|).
\qquad\text{(11)}
```

This bound does not divide by $p_s$. In particular, it remains valid when an outcome becomes very rare or vanishes as $\theta\to0$.

### 2.3 Polar recovery without loss of exponent

The polar recovery of [Proposition 1](model.md) has positive Kraus operators $A_s=\sqrt{M_s}$ and channel $\Phi(\rho)=\sum_sA_s\rho A_s$. Since $\sum_sA_s^2=I_D$,

```math
\Phi-\mathrm{id}
=-\frac12\sum_s[A_s,[A_s,\,\cdot\,]].
\qquad\text{(12)}
```

Set $C_s=A_s-\sqrt{p_s}I_D$. Scalar shifts do not change commutators. Left or right multiplication by an operator $C$ has diamond norm at most $\Vert C\Vert _\infty$, so $\Vert [C,\,\cdot\,]\Vert _\diamond\leq2\Vert C\Vert _\infty$. It follows that

```math
\frac12\Vert \Phi-\mathrm{id}\Vert _\diamond
\leq\sum_s\Vert C_s\Vert _\infty^2.
\qquad\text{(13)}
```

The scalar matrix $p_sI_D$ commutes with $M_s$. For every eigenvalue $\lambda\geq0$ of $M_s$,

```math
(\sqrt\lambda-\sqrt{p_s})^2\leq|\lambda-p_s|.
```

Consequently

```math
\Vert C_s\Vert _\infty^2
\leq\Vert M_s-p_sI_D\Vert _\infty.
\qquad\text{(14)}
```

This remains true for $p_s=0$; in that case positivity and zero trace actually imply $M_s=0$. Equations (11)–(14) prove the second inequality of (5). The first follows from the general fidelity–diamond comparison, and the bound by one follows because both channels are trace preserving. As $B_S$ starts at power $r$, (6) follows. $\square$

The squared commutator bound is important: estimating square roots alone and then using a generic channel-distance estimate can lose half the exponent. Here the channel identity (12) retains the correct order without any lower bound on syndrome probabilities.

## 3. An embedded even repetition witness

Assume $r<\infty$. By [Theorem 2](support_criterion.md), choose an even logical circuit

```math
E=\prod_{j\in J}E_j,
\qquad |J|=r,
```

where the local factors $E_j$ are Hermitian nonidentity Paulis and $E$ is a nontrivial logical operator. Take the fixed Hamiltonian tuple $h_j=E_j$ on $J$ and $h_j=0$ elsewhere. For $0\leq\theta\leq\pi/2$, let

```math
c=\cos(\theta/2),\qquad
s=\sin(\theta/2),\qquad
q=s^2\leq\frac12.
```

The physical rotation has the exact expansion

```math
U_\theta
=\sum_{A\subseteq J}c^{r-|A|}(-is)^{|A|}E_A.
\qquad\text{(15)}
```

### 3.1 Only complementary subsets share a syndrome

The words $E_A$ and $E_B$ have the same syndrome precisely when $E_{A\triangle B}\in N$. The circuit property says this happens exactly when $A\triangle B$ is empty or all of $J$. Thus the error words in (15) fall into $2^{r-1}$ distinct complementary pairs $\{A,J\setminus A\}$, one per syndrome label reached by these words. At $\theta=0$, some such branches have zero amplitude; this does not change the pairing statement.

For a pair represented by $A$ of weight $w$, use $E_A$ as the Pauli return from that syndrome sector. Since disjoint local factors commute and square to identity, the corresponding logical Kraus operator is

```math
K_A
=c^{r-w}(-is)^w I_D
+c^w(-is)^{r-w}\overline E,
\qquad\text{(16)}
```

where $\overline E=V^\dagger EV$ is a nontrivial Hermitian logical Pauli. Its two eigenspaces have dimension $D/2$.

Because $r$ is even, the relative phase of the coefficients in (16) is real. After removing a common phase, this is $aI_D+\eta b\overline E$, with $a,b\geq0$ and $\eta\in\{+1,-1\}$. Thus

```math
\frac1D\mathrm{Tr}\sqrt{K_A^\dagger K_A}
=\frac{|a+b|+|a-b|}{2}
=\max\{a,b\}.
\qquad\text{(17)}
```

The branch contributes $\max\{a^2,b^2\}$ to optimal entanglement fidelity. Since $\sum_{\{A,A^c\}}(a^2+b^2)=1$, the optimal infidelity is the sum of the smaller probability in each pair.

### 3.2 Exact optimal error

For $q\leq1/2$, the smaller probability belongs to the larger-weight subset. In the equal-weight case, a complementary pair is counted twice among subsets of size $r/2$. Therefore

```math
\begin{aligned}
p_r(q)
&=\sum_{w>r/2}\binom rwq^w(1-q)^{r-w}\\
&\quad+\frac12\binom r{r/2}q^{r/2}(1-q)^{r/2}
\end{aligned}
\qquad\text{(18)}
```

satisfies $e_{\mathrm{opt}}(U_\theta)=p_r(q)$. This is the majority-decision error of an even-length classical repetition code, including a half-weight tie.

In fact the optimal diamond error is exactly the same:

```math
\boxed{
\varepsilon_{\mathrm{opt}}(U_\theta)
=e_{\mathrm{opt}}(U_\theta)
=p_r\!\left(\sin^2(\theta/2)\right).
}
\qquad\text{(19)}
```

**Proof of diamond equality.** Choose the lower-weight word as the return for every complementary pair, breaking equal-weight ties arbitrarily. The corrected branch has the form $\alpha_A I_D+\beta_A\overline E$, with $|\alpha_A|\geq|\beta_A|$ and $\alpha_A\beta_A^*$ real. Its unwanted squared amplitude is $|\beta_A|^2$, and their sum is $p_r(q)$.

Trace preservation gives

```math
\sum_A(|\alpha_A|^2+|\beta_A|^2)=1,
\qquad
\sum_A\mathrm{Re}(\alpha_A\beta_A^*)=0.
```

Because every product in the second sum is real, the cross terms cancel entirely. The recovered logical channel is consequently

```math
\Lambda_\theta(\rho)
=(1-p_r(q))\rho+p_r(q)\overline E\rho\overline E.
\qquad\text{(20)}
```

A nontrivial Hermitian Pauli has $\Vert \mathrm{Ad}_{\overline E}-\mathrm{id}\Vert _\diamond=2$. The upper bound is the general distance between channels; equality is attained by a balanced superposition of one $+1$ and one $-1$ eigenvector, which $\overline E$ sends to an orthogonal state. Hence (20) has half diamond error $p_r(q)$. Every recovery has half diamond error at least its entanglement infidelity and thus at least the optimal infidelity already calculated. This proves (19). $\square$

The repetition-code weak-measurement mechanism is established in the coherent-error literature, notably Huang, Doherty, and Flammia. Its role here is the universal embedding: every finite-$r$ stabilizer code contains a minimum circuit that realizes (18), making the upper exponent sharp even against arbitrary post-syndrome recovery. See [the attribution comparison](../references/RELATED_WORK.md).

## 4. Small-angle and finite-angle consequences

The tie term is the leading contribution in (18). Since $q=\theta^2/4+O(\theta^4)$,

```math
p_r\!\left(\sin^2(\theta/2)\right)
=\frac12\binom r{r/2}\left(\frac\theta2\right)^r
+O(\theta^{r+2}).
\qquad\text{(21)}
```

The coefficient is positive. Equations (6) and (21), with $e_{\mathrm{opt}}\leq\varepsilon_{\mathrm{opt}}$, prove both statements in (3). Reversing the sign of $\theta$ leaves the witness errors unchanged, so the result is stated in $|\theta|$.

At the finite angle $\theta=\pi/2$, binomial symmetry gives $p_r(1/2)=1/2$. The witness then reduces the logical channel to complete dephasing in the $\overline E$ eigenspaces. For $r=4$,

```math
p_4(q)=3q^2-2q^3,
\qquad
p_4\!\left(\sin^2(\theta/2)\right)
=\frac3{16}\theta^4+O(\theta^6).
\qquad\text{(22)}
```

This is the ten-qubit example verified directly in [the code-family account](code_families.md).

## 5. What the exponent establishes

The exact exponent is a statement about information irreversibly removed by a fixed measurement instrument. It distinguishes codes even when their conventional parameters agree, and it determines the first possible order of optimal recovery error under weak bounded local rotations.

The finite-code upper constant $B_S(1)$ is allowed to grow with the code. Consequently, (3) by itself supplies no uniform threshold or fixed-angle suppression for a growing family. The exact witness formula (18) describes one adversarial axis choice and does not supply a general efficient decoder. Noise in the syndrome-extraction circuit, uncertainty about the rotation, and restrictions on available recovery gates change the problem specified in [the model](model.md).
