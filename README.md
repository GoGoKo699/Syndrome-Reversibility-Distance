# Syndrome Reversibility Distance

**When does an error-syndrome measurement irreversibly reveal encoded quantum information?** This project gives an exact answer for stabilizer codes exposed to known products of single-qubit rotations. A support invariant determines exact reversibility and the sharp small-angle recovery exponent. Codes with identical ordinary parameters can have different values of this invariant.

The operation order is fixed: encode, apply the physical product unitary, measure the full **original** stabilizer syndrome, then recover. Recovery may know the entire unitary and receives both the syndrome and the remaining quantum state.

## Main results

Let $S$ be the full stabilizer group and $N$ its Pauli centralizer, with membership understood modulo global Pauli phases. A nontrivial logical Pauli belongs to $N\setminus S$.

**Exact reversibility.** Every product unitary is perfectly recoverable after the prescribed measurement if and only if no stabilizer and nontrivial logical Pauli have exactly the same physical support:

```math
\{\mathrm{supp}(g):g\in S\}
\cap
\{\mathrm{supp}(L):L\in N\setminus S\}
=\varnothing.
```

**Sharp recovery exponent.** Define the syndrome reversibility distance by

```math
r(S)=\min_{\substack{g\in S,\ L\in N\setminus S\\
\mathrm{supp}(g)=\mathrm{supp}(L)}}
\mathrm{wt}(gL),
\qquad \min\varnothing=\infty.
```

For each fixed code with finite $r$, the worst-case optimal half-diamond error and the worst-case optimal entanglement infidelity under bounded local rotations of angle scale $\theta$ are both $\Theta(\theta^r)$. The converse allows every syndrome-dependent quantum recovery. Equivalently, $r$ is the minimum even weight of a nontrivial logical Pauli whose proper nonempty Pauli subwords all have nonzero syndrome. Thus finite $r$ is even, $r\ge d$, and even ordinary distance implies $r=d$. An arbitrary even logical word need not be such a circuit.

**Same code parameters, different reversibility.** For every concatenation depth $m\ge1$, the following codes have parameters $[[2\cdot5^m,2,3^m]]$:

| Code | Syndrome reversibility distance |
|---|---:|
| Two independent concatenated five-qubit blocks | $\infty$ |
| The same code conjugated by one physical CNOT between blocks | $3^m+1$ |

The CNOT changes the code's position relative to the physical tensor-product structure while retaining its length, encoded dimension and ordinary distance. In the smallest ten-qubit example, an explicit four-site rotation has optimal entanglement infidelity and optimal half-diamond error

```math
p_4(q)=3q^2-2q^3,
\qquad q=\sin^2(\theta/2),\quad 0\le\theta\le\pi/2.
```

## Protected gates with the same readout

The sister project [Error-Transparent Locality](https://github.com/GoGoKo699/Error-Transparent-Locality) asks how much native interaction order and coupling an exact protected Hamiltonian requires. Their [joint theorem](science/gate_readout_tradeoff.md) states that, for odd $d\ge3$, $t\ge1$ and $d\ge2t+1$, a nontrivial gate transparent to every Pauli error of weight at most $t$ at interaction order $\ell=d+t$ forces $r(S)=d+1$. Hence $r(S)>d+1$, including infinity, requires $\ell\ge d+t+1$.

The concatenated five-qubit family $Q_m=[[5^m,1,3^m]]$ has $r(S)=\infty$ and nevertheless supports exact protected local gates and entanglers. With $d=3^m$, its family-specific minimum orders and native coupling optima are:

| Protected physical error weight | Minimum interaction order | Optimal axis budget $J/\omega$ | Depth |
|---|---:|---:|---|
| One | $d+2$ | $d+1$ | $m\ge1$ |
| Two | $d+4$ | $(4d^2+21d+8)/8$ | $m\ge2$ |
| Three | $d+6$ | $J_3(d)=(4d^3+51d^2+170d+24)/24$ | $m\ge2$ |

Here $\omega$ is half the logical spectral width and $J$ is the sum of absolute nonidentity bare Hermitian Pauli coefficients. These statements use exact sector intertwining. They concern the stated orders; they do not claim optimality at every larger order. The [science guide](science/README.md#protected-control-with-infinite-reversibility-distance) separates the dense constructions, sparse synthesis, coupling certificates and finite checks.

On a fixed product of $Q_m$ blocks with $m\ge2$, prescribe $h_L=\Omega P+cI$ with $\Omega\ne0$ and a nonidentity logical Pauli string $P$ of weight $p\ge1$. It has a three-error-transparent realization at minimum order $pd+6$ with optimal budget $|\Omega|J_3(pd)$. The target coset minimum $pd$ differs from ordinary code distance $d$. In particular, a protected logical $ZZ$ gate on two 25-qubit blocks attains order 24 and budget $1789|\Omega|$. Protected single-block controls and this entangler give universal encoded control in the native Hamiltonian model, while the product code and original syndrome measurement retain $r(S)=\infty$.

The physical CNOT used for the code-family separation **conjugates the code**. The protected logical entangler **keeps the product code fixed**. The readout experiment and protected-gate noise model have separate assumptions; no joint threshold or hardware-efficiency claim is made.

## Read and learn

Start with the [science guide](science/README.md) for a topic-based proof route, the full result catalog and the limits of each construction.

| Question | Authoritative account |
|---|---|
| What channel is measured, and what may recovery use? | [Model and recovery](science/model.md) |
| Why do exact supports and minimal logical circuits classify reversibility? | [Support criterion](science/support_criterion.md) |
| Why is the optimized exponent sharp, including rare syndromes? | [Error exponent](science/error_exponent.md) |
| How can one physical CNOT preserve distance but change reversibility? | [Code families](science/code_families.md) |
| What does the ordering mean physically? | [Interpretation](science/interpretation.md) |
| Which ingredients come from prior research? | [Related work](references/RELATED_WORK.md) and [bibliography](references/references.bib) |

The selected educational anchor is **Daniel Gottesman, _An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation_**, arXiv:0904.2557v1 (2009; published 2010). The [Gottesman bridge](learning/GOTTESMAN_BRIDGE.md) maps the source sections and notation to this project and supplies the additional channel-recovery arguments. The [30-milestone map](learning/MILESTONES.md) has 12 background/preliminaries milestones followed by 18 narrative/innovation milestones, each with an understanding target, earlier dependencies, a conceptual connection and source pointers. Use the [learning guide](learning/README.md) to enter this common route. The [LLM guide](llms.txt) identifies the authoritative files and distinctions to preserve.

The general outcome-conditioned recovery criterion and optimal-fidelity formula are established results of Gregoratti and Werner. The even-repetition mechanism is established in work of Huang, Doherty and Flammia. The [attribution account](references/RELATED_WORK.md) distinguishes these ingredients from the support/circuit classification and its operational consequences proved here.

## Reproduce the checks

From the repository root:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_all.py
```

The [verification guide](verification/README.md) describes all 12 programs, their exact algebraic certificates, numerical physical-channel checks and limits. These include direct physical-Hilbert-space recovery calculations for the ten-qubit example, binary Pauli calculations, protected-gate response checks and optimality certificates. The universal classification and growing families rely on the written proofs.

## Scope

The syndrome measurement is full, ideal and fixed to the original stabilizer code. Recovery occurs afterward, may depend on the complete physical unitary, and has unrestricted computational and quantum-processing resources. Small-angle constants may depend on the fixed code. Exact protected Hamiltonians use native Pauli order and coefficient cost, without a geometric layout, efficient synthesis or fault-tolerant threshold assumption.

## Purpose and contact

This repository serves as a record of the work and a guide for the author’s self-directed learning. For discussion or potential collaboration, please contact Ruge Lin at [gogoko699@gmail.com](mailto:gogoko699@gmail.com).
