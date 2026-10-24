# Syndrome Reversibility Distance

**When does an error-syndrome measurement irreversibly reveal encoded quantum information?** This project gives an exact answer for stabilizer codes exposed to known products of single-qubit rotations. It identifies a code invariant that controls the best possible recovery after the prescribed measurement and is not determined by the usual code parameters.

The operation order is part of the question: encode, apply the physical product unitary, measure the full original stabilizer syndrome, then recover. Recovery may depend on the unitary and has access to the complete syndrome and remaining quantum state.

## Main results

Let $S$ be the full stabilizer group, and let $N$ be its Pauli centralizer. Group membership below ignores overall Pauli phases. A nontrivial logical Pauli is an element of $N\setminus S$.

**Exact reversibility.** Every product unitary is perfectly recoverable after syndrome measurement if and only if no stabilizer and nontrivial logical Pauli have exactly the same physical support:

```math
\{\operatorname{supp}(g):g\in S\}
\cap
\{\operatorname{supp}(L):L\in N\setminus S\}
=\varnothing.
```

**Sharp error exponent.** Define the syndrome reversibility distance by

```math
r(S)=\min_{\substack{g\in S,\ L\in N\setminus S\\
\operatorname{supp}(g)=\operatorname{supp}(L)}}
\operatorname{wt}(gL),
\qquad \min\varnothing=\infty.
```

For a fixed code with finite $r$, the worst-case optimal half-diamond error after bounded local rotations of angle scale $\theta$ is $\Theta(\theta^r)$. The matching converse permits every syndrome-dependent quantum recovery. An equivalent description is the minimum even weight of a nontrivial logical Pauli whose proper nonempty Pauli subwords all have nonzero syndrome. In particular, finite $r$ is even and satisfies $r\ge d$; every even-distance code has $r=d$.

**Identical code parameters, different reversibility.** For every integer $t\ge1$, two families have the same parameters $[[2\cdot5^t,2,3^t]]$:

| Code | Syndrome reversibility distance |
|---|---:|
| Two independent, $t$-fold concatenated five-qubit blocks | $\infty$ |
| The same code conjugated by one physical CNOT between the blocks | $3^t+1$ |

The CNOT changes the code's position relative to the physical tensor-product structure. The number of physical qubits, number of encoded qubits, and ordinary distance remain unchanged.

The smallest pair consists of two $[[10,2,3]]$ codes. For an explicit four-site rotation on the CNOT-conjugated code, both optimal entanglement infidelity and optimal half-diamond error are exactly

```math
p_4(q)=3q^2-2q^3,
\qquad q=\sin^2(\theta/2),\quad 0\le q\le\tfrac12.
```

## Read the science

| Question | Authoritative account |
|---|---|
| What is retained, what can recovery do, and which errors are optimized? | [Model and recovery](science/model.md) |
| Why is equality of exact supports the criterion? | [Support classification](science/support_criterion.md) |
| Why does the invariant give the optimal exponent, including rare syndromes? | [Sharp error exponent](science/error_exponent.md) |
| How does one CNOT preserve distance while changing reversibility? | [Code families and the ten-qubit example](science/code_families.md) |
| What does the result say about syndrome measurement and coherent errors? | [Physical interpretation](science/interpretation.md) |
| Which ingredients are inherited, and how do nearby results compare? | [Related work](references/RELATED_WORK.md) and [bibliography](references/references.bib) |
| How are the finite examples reproduced? | [Verification](verification/README.md) |

The general outcome-conditioned recovery criterion and optimal-fidelity formula are established results of Gregoratti and Werner. The even-repetition error mechanism is established in work of Huang, Doherty and Flammia. Their roles, and the connection to more recent coherent-error studies, are identified in the [attribution account](references/RELATED_WORK.md). The proofs here connect that recovery framework to the physical supports of arbitrary stabilizer codes.

## Reproduce the checks

From the repository root:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_all.py
```

The checks include exact binary Pauli calculations and a direct physical-Hilbert-space calculation for all syndrome sectors of the ten-qubit example. The universal classification and growing-distance family are established by the written proofs; finite computations check their specified instances.

## Scope

The syndrome measurement is full, ideal, and fixed to the original stabilizer code. Recovery occurs afterward, may depend on the full physical unitary, and has unrestricted computational and quantum-processing resources. The small-angle estimates hold for each fixed code; their constants are not asserted to remain uniform over a growing family. The [interpretation](science/interpretation.md) explains how these assumptions enter the result.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).

