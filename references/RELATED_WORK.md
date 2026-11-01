# Related work and attribution

This project asks whether a prescribed stabilizer syndrome measurement preserves the information in a code after a known product unitary. The product of physical single-qubit unitaries acts first, the original full syndrome is measured ideally next, and arbitrary recovery follows with knowledge of the unitary and the outcome. The recovery is optimized separately for each specified unitary. This is the channel and ordering fixed in [the model](../science/model.md).

## Results and inherited ingredients

Let $S$ be the phase-free stabilizer label space and $N$ its Pauli centralizer. Define

```math
r_S=\min_{\substack{g\in S,\ L\in N\setminus S\\
\mathrm{supp}(g)=\mathrm{supp}(L)}}\mathrm{wt}(gL),
\qquad \min\varnothing=\infty.
```

The science pages denote this invariant by $r(S)$ or $r$. It uses exact physical support and local Pauli labels, not ordinary code distance alone.

| Result proved here | Inherited ingredients and additional conclusion |
|---|---|
| [Exact support criterion](../science/support_criterion.md) | Record-assisted recovery is established. The classification identifies precisely which stabilizer codes are recoverable after the prescribed measurement for every product unitary: no stabilizer and nontrivial logical Pauli may have exactly equal support, equivalently $r_S=\infty$. It applies to arbitrary numbers of logical qubits. |
| [Logical-circuit characterization and sharp exponent](../science/error_exponent.md) | The even-repetition dephasing mechanism is established. Here $r_S$ equals the minimum even weight of a nontrivial logical Pauli with no proper nonempty normalizer subword. For a fixed code with finite $r_S$, both worst-case optimal entanglement infidelity and optimal half diamond error are $\Theta(|\theta|^{r_S})$ for bounded single-qubit generators. Constants may depend on the code. |
| [Equal-parameter families](../science/code_families.md) | The five-qubit code and concatenation are established. Two concatenated blocks and the same code conjugated by one physical CNOT have identical $[[2\cdot5^m,2,3^m]]$ parameters but $r_S=\infty$ and $r_S=3^m+1$, respectively. Each is measured in its own prescribed original syndrome. |

Finite $r_S$ is even and at least $d$; at even distance it equals $d$. These statements concern irreversible information loss after the fixed measurement, with the known unitary available to recovery. They do not assert that a common decoder corrects unknown rotations or that code-family constants give a threshold. The source comparisons below identify the recovery, parity, coherent-noise, and code ingredients. Bibliographic records are in [references.bib](references.bib).

## Gregoratti–Werner: record-assisted recovery

M. Gregoratti and R. F. Werner, *Quantum lost and found*, Journal of Modern Optics **50**, 915–933 (2003). [Article](https://doi.org/10.1080/09500340308234541) · [arXiv version](https://arxiv.org/abs/quant-ph/0209025).

For a fixed observed Kraus decomposition, Proposition 2 gives exact recovery precisely when each Kraus effect is a scalar identity. Proposition 8, Section V, gives the maximum entanglement fidelity and recovery by polar decomposition, including completion on the unused output subspace. Eqs. (16)–(20) specialize directly to the fidelity formula and positive-Kraus channel in [the model](../science/model.md).

These recovery facts are inherited. The support theorem classifies when their hypothesis holds for every product unitary on a stabilizer code; the exponent theorem quantifies failure. Polar recovery is an entanglement-fidelity optimizer and a constructive diamond-error upper bound. It is not asserted to minimize diamond error for every instrument.

## Cheng et al.: symmetry-protected logical unitarity

Z. Cheng, E. Huang, V. Khemani, M. J. Gullans, and M. Ippoliti, *Emergent unitary designs for encoded qubits from coherent errors and syndrome measurements*, PRX Quantum **6**, 030333 (2025). [Article](https://doi.org/10.1103/bnld-2chd) · [arXiv v2](https://arxiv.org/html/2412.04414v2).

Theorem 2 proves syndrome-conditioned logical unitarity for a one-logical-qubit CSS code with even-weight stabilizers and odd X/Z distances, when the physical unitary commutes with the specified time-reversal operation. This includes all products of single-qubit unitaries. Section VI.2.3 discusses other LDPC codes and multiple logical qubits.

The support criterion gives necessity and sufficiency for arbitrary qubit stabilizer codes and any number of logical qubits within the product-unitary class; its five-qubit example is non-CSS. The physical classes differ: the time-reversal theorem also covers certain entangling unitaries, while the classification here is specifically for products.

## Huang–Doherty–Flammia: coherent-error scaling

E. Huang, A. C. Doherty, and S. T. Flammia, *Performance of quantum error correction with coherent errors*, Physical Review A **99**, 022313 (2019). [Article](https://doi.org/10.1103/PhysRevA.99.022313) · [arXiv v2](https://arxiv.org/html/1805.08227v2).

Theorem 1 proves a code-dependent bound of order the physical rotation strength raised to distance $d$, using minimum-weight decoding for independent unitary errors. Section III.2, especially Eqs. (31)–(33), derives the even-repetition dephasing probability, including the half-weight contribution of tied error patterns, and interprets syndrome outcomes as weak logical measurements.

The repetition tail is inherited. Its use here is an embedding in every code with finite $r_S$: a minimum even logical circuit produces the complementary error pairs of an even repetition code. The converse optimizes over every subsequent CPTP recovery, and the matching upper bound identifies $r_S$ as the exact worst-case exponent. For a four-site witness, the exact error is $p_4(q)=3q^2-2q^3$, where $q=\sin^2(\theta/2)$ and $0\le q\le1/2$. Rotation conventions differ by a factor of two: here $U_j=\exp(-i\theta E_j/2)$.

## Liu–Chen: syndrome instruments and logical projections

H. Liu and X. Chen, *Coherent error induced phase transition*, [arXiv:2506.00650v3](https://arxiv.org/html/2506.00650v3) (2026 revision); [accepted Physical Review B article](https://journals.aps.org/prb/accepted/10.1103/jly5-q9mx).

Section II.3 distinguishes state-specific return probabilities from channel-level recoverability. Section II.4, Proposition 1, equates logical branch unitarity, scalar branch effects, equal nonzero singular values, and input-independent syndrome probabilities. Its Clifford analysis, Eqs. (57)–(67), gives the logical-projector form and the number of logical Pauli observables learned by the syndrome.

These are precedents for the instrument interpretation and product-Clifford witness. The support theorem classifies which stabilizer codes satisfy scalar branch effects for all product rotations; the invariant determines the optimal weak-angle loss. The $r$ denoting measurement rank in Liu–Chen is different from the physical Pauli weight $r_S$ here.

## Coherent-noise diagnostics, thresholds, and parity

S. J. Beale, J. J. Wallman, M. Gutiérrez, K. R. Brown, and R. Laflamme, *Quantum Error Correction Decoheres Noise*, Physical Review Letters **121**, 190501 (2018). [Article](https://doi.org/10.1103/PhysRevLett.121.190501) · [arXiv v2](https://arxiv.org/html/1805.08802v2).

Theorem 2 bounds syndrome-averaged off-diagonal logical-noise entries under a product channel in terms of physical infidelity and distance. Suppression of logical coherence is different from the information no recovery can restore after a specified coherent realization: exact reversibility may coexist with nontrivial syndrome-conditioned logical unitaries.

Z. Han, Y. Zhao, Y. Xu, Y. Wang, and Z.-W. Liu, *Coherent error threshold for quantum LDPC codes*, [arXiv:2609.20537v1](https://arxiv.org/html/2609.20537v1) (2026).

Theorems 1–2 and the appendices establish code-capacity threshold bounds for qLDPC families with distance growing at least logarithmically, using cluster resummation for local coherent noise and more general local channels. They treat optimal recovery and minimum-weight decoding. The theorem here is an exact fixed-code weak-angle exponent; its constants may grow with the code, and the concatenated examples do not supply a qLDPC threshold.

E. Kubischta and I. Teixeira, *Time-Reversal Selection Rules for Quantum Error Correction*, [arXiv:2608.06304v1](https://arxiv.org/html/2608.06304v1) (2026).

This work derives scalar compression of even-weight Paulis for a time-reversal-invariant logical qubit on an odd number of spins, including codes beyond stabilizers, and relates time reversal to shadow enumerators. Parity is therefore established background. The support criterion treats arbitrary stabilizer codes without requiring global even/odd weight separation, and its circuit form retains the local Pauli labels needed for the exponent.

C. Cao, M. J. Gullans, B. Lackey, and Z. Wang, *Quantum Lego Expansion Pack: Enumerators from Tensor Networks*, PRX Quantum **5**, 030313 (2024). [Article](https://doi.org/10.1103/PRXQuantum.5.030313) · [arXiv v2](https://arxiv.org/html/2308.05152v2).

Sections III.2–III.3, III.6–III.8, and Appendix D develop generalized and coset enumerators for coherent channels, logical error probabilities, effective-distance diagnostics, and decoding. Refining ordinary distance and evaluating coherent logical channels algebraically are established ideas. The specific operational conclusion here identifies minimum even logical-circuit weight with the optimal worst-case weak-angle exponent and classifies its infinite case by exact support. This is not the effective-distance definition of their Section III.3.

## Five-qubit code and stabilizer foundations

The five-qubit code is established work. R. Laflamme, C. Miquel, J. P. Paz, and W. H. Zurek introduced a perfect five-qubit code in *Perfect Quantum Error Correcting Code*, Physical Review Letters **77**, 198–201 (1996): [article](https://doi.org/10.1103/PhysRevLett.77.198), [arXiv](https://arxiv.org/abs/quant-ph/9602019). A five-qubit construction also appears in C. H. Bennett, D. P. DiVincenzo, J. A. Smolin, and W. K. Wootters, *Mixed-state entanglement and quantum error correction*, Physical Review A **54**, 3824–3851 (1996): [article](https://doi.org/10.1103/PhysRevA.54.3824), [arXiv v2](https://arxiv.org/abs/quant-ph/9604024v2).

For the explicit stabilizer convention, see Table 4 and the binary symplectic representation in Section 3.1 of D. Gottesman, *An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation*, Proceedings of Symposia in Applied Mathematics **68**, 13–58 (2010), [arXiv:0904.2557v1](https://arxiv.org/abs/0904.2557v1). Sections 2.2–2.3 develop error-correction conditions, stabilizers, syndrome measurement, and distance.

The [family proof](../science/code_families.md) establishes the required concatenated parity and distance properties, then analyzes the change caused by one physical CNOT. Neither the underlying code nor concatenation is claimed as new. For unchanged product blocks, equality of a stabilizer support and a nontrivial logical support is excluded on an active block. Global weight parity alone would not exclude logical strings active on an even number of blocks.

## Coexistence with protected logical generators

The companion [Error-Transparent Locality account](https://github.com/GoGoKo699/Error-Transparent-Locality) studies static Hamiltonians with equal logical rates on every Pauli-error image through a specified weight. Its native resources are maximum physical Pauli-term weight $\ell$ and the sum $J$ of absolute nonidentity coefficients. It uses the inherited five-qubit code, canonical protected base Hamiltonians, and stabilizer dressing. The joint [gate-readout theorem](../science/gate_readout_tradeoff.md) combines its interaction-weight moments with the circuit theorem here: at odd $d$, attaining the general gate minimum $\ell=d+t$ forces $r_S=d+1$, so $r_S>d+1$ requires $\ell\ge d+t+1$. Fixed ancillary qubits do not evade the theorem.

The following coexistence constructions retain the same code and original syndrome, with $r_S=\infty$. Put $Q_m=[[5^m,1,3^m]]$ and $d=3^m$.

| Protected generator on unchanged codes | Derived gate statement |
|---|---|
| [One-error gates](../science/minimum_readout_safe_gates.md) on $Q_m$ | Each logical Pauli axis attains $\ell=d+2$ and $J/\omega=d+1$. |
| [Two-error gates](../science/growing_two_error_gates.md), $m\ge2$ | Correlated inactive-frontier padding and block moments give $\ell=d+4$ and $J/\omega=(4d^2+21d+8)/8$. The [four-class twenty-five-qubit case](../science/two_error_gates.md) has order thirteen and budget $521/8$. |
| [Three-error gates](../science/growing_three_error_gates.md), $m\ge2$ | A strict total-degree increment lifting and cubic block moments give $\ell=d+6$ and $J/\omega=J_3(d)=(4d^3+51d^2+170d+24)/24$. The [seven-class twenty-five-qubit case](../science/three_error_gates.md) has order fifteen and budget $2867/8$. |
| [Pure logical Pauli strings](../science/protected_entangling_gates.md) on $Q_m^{\otimes M}$, $m\ge2$ | A target of logical weight $p$ has coset minimum $D=pd$, distinct from ordinary distance $d$, and three-error optima $\ell=D+6$, $J=|\Omega|J_3(D)$ for $h_L=\Omega P+cI$. |

These are minimum orders and minimum budgets at those orders for the specified families and targets. The three-error proof includes a full error-image zero-factor argument; the pure-string lower bound projects moments onto the target logical-Pauli coefficient rather than restricting competing Hamiltonians to a tensor ansatz. The finite protected entangler and Pauli Lie-algebra criterion give universal encoded control through complete protected Hamiltonians. They do not give a common weighted-sum optimum for arbitrary mixed logical Hamiltonians at one locality cap or a physical synthesis of many-body interactions.

Description size is treated separately. Standard nonnegative support reduction, discussed in [Hank Yang's background notes](https://hankyang.seas.harvard.edu/Semidefinite/background.html), gives $3n+1$ single-error terms in the [sparse account](../science/sparse_protected_gates.md). A constructive cyclic-orbit recursion gives at most $4n$ with $O(\log n)$-bit coefficients, including a 75-term twenty-five-qubit witness. The [two-error sparse compiler](../science/sparse_two_error_synthesis.md) combines positive component certificates, pairwise-independent coupling, and rational torus polygon areas, including an exact 125-qubit witness. Numerical support selection is followed by independent rational and integer verification. Three-error polynomial-size existence does not yet provide a practical sparse compiler.

The gate and readout experiments keep their separate assumptions. Bounded-weight Pauli faults during protected continuous evolution and a known product unitary before ideal syndrome readout are different guarantees. Their coexistence and encoded universality do not establish a joint noisy threshold, hardware feasibility, or exhaustive literature priority. This comparison positions the results against the specified sources; it does not assert that every related formulation has been located.
