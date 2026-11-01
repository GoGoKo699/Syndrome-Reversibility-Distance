# Science guide

The central object is the channel produced when a known product of single-qubit unitaries is followed by the full original stabilizer-syndrome measurement. The project classifies when that measurement loses logical information, quantifies the smallest unavoidable loss, and separates this information from ordinary code distance. Protected Hamiltonians then show how infinite reversibility distance can coexist with encoded control on the same code.

## The readout proof route

| Step | What to establish | Source |
|---|---|---|
| Specify the retained information | Recovery receives the classical syndrome and the remaining quantum state, knows the unitary, and acts only after measurement. Branch effects determine logical information in the outcomes. | [Model and optimal recovery](model.md) |
| Classify exact reversibility | Every branch effect is scalar for every product unitary exactly when stabilizer and nontrivial logical support families are disjoint. Equal support also gives a constructive witness. | [Support criterion](support_criterion.md) |
| Identify the invariant | Minimize the weight of the product of an equal-support stabilizer/logical pair. Equivalently, minimize the even weight of a logical circuit with no proper nonempty zero-syndrome subword. | [Circuit characterization](support_criterion.md) |
| Prove the sharp exponent | Collision polynomials and polar recovery give the upper bound; an exact even-repetition witness gives a converse against every recovery. | [Error exponent](error_exponent.md) |
| Separate it from ordinary distance | Concatenated five-qubit blocks have disjoint support parities; one physical CNOT creates a minimal collision without changing the code parameters. | [Code families](code_families.md) |
| Interpret the conclusion | The invariant depends on physical tensor factors, operation order and the prescribed measurement. Fixed-code asymptotics do not provide uniform family constants. | [Physical interpretation](interpretation.md) |

Write $r(S)$, or $r_S$, for the same invariant throughout. It is infinite precisely when every known product unitary is exactly recoverable after the prescribed measurement. When finite it is even and at least ordinary distance $d$; for even $d$ it equals $d$. Circuit minimality refers to subwords of fixed Pauli labels. An arbitrary even logical word is insufficient.

The family comparison uses $Q_m=[[5^m,1,3^m]]$. Two independent blocks have $r_S=\infty$; a single physical CNOT conjugation gives $r_S=3^m+1$ at the same parameters $[[2\cdot5^m,2,3^m]]$. The ten-qubit rotation witness has both optimized errors exactly $3q^2-2q^3$, with $q=\sin^2(\theta/2)$ and $0\le\theta\le\pi/2$.

## Protected control with infinite reversibility distance

[Error-Transparent Locality](https://github.com/GoGoKo699/Error-Transparent-Locality) studies exact protected Hamiltonians on a fixed code. Its native coupling cost is $J=\sum_{Q\ne I}|h_Q|$ in the bare Hermitian tensor-Pauli expansion. Write $\omega$ for half the logical spectral width. The budget differs from operator norm and does not include geometry or implementation overhead. The shared accounts below are included here so the coexistence result can be read independently.

For odd $d\ge3$, $t\ge1$ and $d\ge2t+1$, the [joint theorem](gate_readout_tradeoff.md) proves that attaining the general interaction-order bound $\ell=d+t$ forces $r_S=d+1$. Requiring $r_S>d+1$ raises the necessary order to at least $d+t+1$. The nine-qubit four-body gate has $r_S=4$ and an explicit four-site readout witness. Removing fixed weight-one ancillas preserves $d$, $k$, $r_S$, the logical generator and transparency without increasing interaction order; merely padding a code with fixed qubits does not make an even logical word into a minimal circuit.

On $Q_m$, all representatives of a nontrivial logical coset have odd weight. For $2t+1\le d$, this stronger structure forces $\ell\ge d+2t$. The following constructions attain that bound, with matching coupling certificates at that order.

| Topic | Exact result and certificate | Source |
|---|---|---|
| One protected error | All three axes attain $\ell=d+2$ and $J/\omega=d+1$ for $m\ge1$ using two positive logical shells with signed masses. | [Minimum readout-safe gates](minimum_readout_safe_gates.md) |
| Sparse one-error gates | An explicit construction uses at most $4n$ native terms per axis, with a separate $3n+1$ existence bound; the 25-qubit example uses 75 terms. | [Sparse protected gates](sparse_protected_gates.md) |
| Two protected errors | The 25-qubit gate has $\ell=13$ and $J/\omega=521/8$; a four-pattern moment certificate proves the optimum. | [Two-error gates](two_error_gates.md) |
| Growing two-error family | For every $m\ge2$, $\ell=d+4$ and $J/\omega=(4d^2+21d+8)/8$. Correlated padding follows the concatenation tree. | [Growing two-error gates](growing_two_error_gates.md) |
| Exact sparse two-error synthesis | Positive component certificates and rational pairwise coupling give a general polynomial reference compiler. A separate support search supplies the 125-qubit witness with 69,125 terms per axis, order 31 and budget $3491/8$; the full depth-three reference compiler was not used to generate it. Its native common denominator has 1,458 bits. | [Sparse two-error synthesis](sparse_two_error_synthesis.md) |
| Three protected errors | The 25-qubit gate has $\ell=15$ and $J/\omega=2867/8$, with seven-pattern primal and dual certificates. It has 1,480,960 terms per axis and a 22-bit common denominator. | [Three-error gates](three_error_gates.md) |
| Growing three-error family | Strict truncated-tensor lifting attains $\ell=d+6$ and $J/\omega=J_3(d)=(4d^3+51d^2+170d+24)/24$ for all $m\ge2$. In the fixed canonical logical frame, the exact single-block direction cost is $J_3(d)\sum_\alpha\lvert\Omega_\alpha\rvert$. | [Growing three-error gates](growing_three_error_gates.md) |
| Protected Pauli strings and entanglers | On fixed product $Q_m$ blocks with $m\ge2$, a target $h_L=\Omega P+cI$ with $\Omega\ne0$ and nonidentity logical Pauli string $P$ of weight $p\ge1$ attains minimum order $pd+6$ and budget $\lvert\Omega\rvert J_3(pd)$. Local controls plus a protected $ZZ$ entangler generate universal encoded control. | [Protected entangling gates](protected_entangling_gates.md) |

Every code in this construction route retains $r_S=\infty$. The three-error lifting combines nested zero- through three-error maps using their increments and a total-degree truncation. It does not preserve every old uniform promotion ensemble. The pure-string coupling optimum is not an optimal weighted budget formula for an arbitrary mixed logical Hamiltonian under one common locality cap.

For the two-block 50-qubit $ZZ$ example, ordinary distance is 9, the target coset minimum is 18, minimum order is 24 and optimal budget is $1789|\Omega|$. The compact dense description specifies 41,032,000,000 distinct native terms with a 35-bit common denominator. No claim of efficient native expansion follows from that compact description.

## What the finite checks certify

The [verification guide](../verification/README.md) distinguishes exact matrix checks, binary Pauli enumeration, compressed response equations, expanded membership checks and analytic certificates. In particular:

| Instance | Coverage | What is not inferred |
|---|---|---|
| Ten-qubit CNOT-conjugated code | All measured syndrome sectors and the direct physical channel for the stated witness | A numerical substitute for the all-code support theorem |
| 125-qubit two-error gate | All 70,126 protected errors through 562 exact response orbits, with all three expanded axes checked | Small coefficient bit complexity |
| 25-qubit three-error gate | All 64,876 protected errors on each axis | A practical sparse three-error compiler |
| 125-qubit three-error gate | Eleven complete profiles cover all 8,649,376 protected errors | Brute-force expansion of the dense Hamiltonian |
| 50-qubit protected entangler | All 540,376 protected errors, with child-axis independence covering all nine logical axis products | Expansion of all 41,032,000,000 native terms |

The all-code and all-depth conclusions come from proofs. Finite checks validate the specified constructions and certificates. Known-product-unitary readout and error-transparent gate evolution remain separate operational models, even when they use the same stabilizer code.

## Sources and learning

The [related-work account](../references/RELATED_WORK.md) separates established recovery results, coherent-error mechanisms and protected-control methods from the claims proved here. Use the [bibliography](../references/references.bib) for source identities and versions.

The common educational anchor is John Preskill, _Quantum Information_, Chapter 7, “Quantum Error Correction”, in the author-hosted draft updated March 2026. This project's [Preskill bridge](../learning/PRESKILL_BRIDGE.md) and [30 milestones](../learning/MILESTONES.md) build the readout proof and then its protected-control connection. The [learning guide](../learning/README.md) explains how to use both. [Error-Transparent Locality has its own Preskill bridge](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/learning/PRESKILL_BRIDGE.md) and milestone map for the protected-gate argument. Return to the [repository overview](../README.md) for the main statements and reproduction commands.
