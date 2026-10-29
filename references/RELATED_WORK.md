# Related work and attribution

The question here is whether a prescribed stabilizer syndrome measurement preserves the information carried by a code after a known product unitary. The unitary acts first, the original full syndrome is measured next, and an arbitrary recovery may act afterward. The quantifiers and channel model are fixed in [the model](../science/model.md).

The contribution consists of the [exact support criterion](../science/support_criterion.md), the [sharp optimal-recovery exponent](../science/error_exponent.md), and the [equal-parameter code families](../science/code_families.md). Several ingredients in their proofs have established precedents. The comparisons below identify those ingredients and the additional statements proved here. Bibliographic records are in [references.bib](references.bib).

## Gregoratti–Werner: record-assisted recovery

M. Gregoratti and R. F. Werner, *Quantum lost and found*, Journal of Modern Optics **50**, 915–933 (2003). [Article](https://doi.org/10.1080/09500340308234541) · [arXiv version](https://arxiv.org/abs/quant-ph/0209025).

For a fixed observed Kraus decomposition, Proposition 2 gives exact recovery precisely when each Kraus effect is a scalar identity. Proposition 8, in Section V, gives the maximum entanglement fidelity and a recovery from polar decomposition, including completion on the unused output subspace. Its Eqs. (16)–(20) specialize directly to the fidelity formula and positive-Kraus channel used in [the model](../science/model.md).

These recovery facts are inherited. Our classification asks when their hypothesis holds for **every** product unitary on a stabilizer code, and then determines the sharp weak-rotation exponent when it fails. Polar recovery is used as an entanglement-fidelity optimizer and a constructive diamond-error upper bound; it is not asserted to minimize diamond error for every instrument.

## Cheng et al.: symmetry-protected logical unitarity

Z. Cheng, E. Huang, V. Khemani, M. J. Gullans, and M. Ippoliti, *Emergent unitary designs for encoded qubits from coherent errors and syndrome measurements*, PRX Quantum **6**, 030333 (2025). [Article](https://doi.org/10.1103/bnld-2chd) · [arXiv v2](https://arxiv.org/html/2412.04414v2).

Theorem 2 proves syndrome-conditioned logical unitarity for a one-logical-qubit CSS code with even-weight stabilizers and odd X/Z distances, when the physical unitary commutes with the specified time-reversal operation. This includes all products of single-qubit unitaries. Section VI.2.3 discusses extensions to other LDPC codes and multiple logical qubits.

The support criterion applies to arbitrary qubit stabilizer codes and any number of logical qubits. It supplies necessity as well as sufficiency within the product-unitary class. The physical noise classes should be compared carefully: the time-reversal theorem also covers certain entangling unitaries, while our classification is specifically for products. The five-qubit example is non-CSS.

## Huang–Doherty–Flammia: coherent-error scaling

E. Huang, A. C. Doherty, and S. T. Flammia, *Performance of quantum error correction with coherent errors*, Physical Review A **99**, 022313 (2019). [Article](https://doi.org/10.1103/PhysRevA.99.022313) · [arXiv v2](https://arxiv.org/html/1805.08227v2).

Theorem 1 proves a code-dependent bound of order the physical rotation strength raised to distance $d$, using minimum-weight decoding for independent unitary errors. Section III.2, especially Eqs. (31)–(33), already derives the even-repetition dephasing probability, including the half-weight contribution from tied error patterns. It explains the syndrome outcomes as weak logical measurements.

Our repetition tail is this established formula. The new use is its embedding in every code with a finite support-collision invariant: a minimum even logical circuit produces exactly the complementary error pairs of an even repetition code. The converse optimizes over all subsequent CPTP recoveries, and the matching upper bound identifies $r$, rather than a general $d$-based estimate, as the exact worst-case exponent. Rotation conventions differ by a factor of two: here $U_j=\exp(-i\theta E_j/2)$.

## Liu–Chen: syndrome instruments and logical projections

H. Liu and X. Chen, *Coherent error induced phase transition*, [arXiv:2506.00650v3](https://arxiv.org/html/2506.00650v3) (2026 revision); [accepted Physical Review B article](https://journals.aps.org/prb/accepted/10.1103/jly5-q9mx).

Section II.3 distinguishes state-specific return probabilities from channel-level recoverability. Section II.4, Proposition 1, equates logical branch unitarity, a scalar branch effect, equal nonzero singular values, and input-independent syndrome probabilities. Its Clifford analysis, Eqs. (57)–(67), gives the logical-projector form and the number of logical Pauli observables learned by the syndrome.

These are direct precedents for the instrument interpretation and the product-Clifford witness used here. The support theorem classifies which stabilizer codes satisfy the scalar-effect condition for **all** product rotations; the invariant then resolves the optimal small-angle loss. The $r$ used for a measurement rank in Liu–Chen is a different quantity from this repository's $r$, which is a physical Pauli weight.

## Beale et al.: suppression of logical coherence

S. J. Beale, J. J. Wallman, M. Gutiérrez, K. R. Brown, and R. Laflamme, *Quantum Error Correction Decoheres Noise*, Physical Review Letters **121**, 190501 (2018). [Article](https://doi.org/10.1103/PhysRevLett.121.190501) · [arXiv v2](https://arxiv.org/html/1805.08802v2).

Theorem 2 bounds the syndrome-averaged off-diagonal entries of the logical noise under a product channel in terms of physical infidelity and code distance. This establishes an important role of syndrome measurement in suppressing logical coherence. The present quantity instead measures information that no recovery can restore after a specified coherent realization. Exact reversibility may coexist with nontrivial syndrome-conditioned logical unitaries, so suppression of coherent channel entries and irreversible loss are distinct diagnostics.

## Han et al.: uniform thresholds for local coherent noise

Z. Han, Y. Zhao, Y. Xu, Y. Wang, and Z.-W. Liu, *Coherent error threshold for quantum LDPC codes*, [arXiv:2609.20537v1](https://arxiv.org/html/2609.20537v1) (2026).

Theorems 1–2 and the appendices establish code-capacity threshold bounds for qLDPC families with distance growing at least logarithmically, using cluster resummation for local coherent noise and more general local channels. Both optimal recovery and minimum-weight decoding are treated.

Our fixed-code theorem gives an exact exponent as the rotation angle tends to zero. Its constants can grow with the code; the concatenated families do not supply a qLDPC threshold result. The comparison therefore separates exact finite-code structure from uniform control across code families.

## Kubischta–Teixeira: time-reversal selection rules

E. Kubischta and I. Teixeira, *Time-Reversal Selection Rules for Quantum Error Correction*, [arXiv:2608.06304v1](https://arxiv.org/html/2608.06304v1) (2026).

This work derives scalar compression of even-weight Paulis for a time-reversal-invariant logical qubit on an odd number of spins, including codes beyond the stabilizer setting. It also relates time reversal to shadow enumerators. The parity mechanism is therefore established background. Our support criterion resolves arbitrary stabilizer codes without assuming that global even/odd weight separation holds, and its circuit form retains the local Pauli labels needed to determine the exponent.

## Cao et al.: enumerators and noise-dependent diagnostics

C. Cao, M. J. Gullans, B. Lackey, and Z. Wang, *Quantum Lego Expansion Pack: Enumerators from Tensor Networks*, PRX Quantum **5**, 030313 (2024). [Article](https://doi.org/10.1103/PRXQuantum.5.030313) · [arXiv v2](https://arxiv.org/html/2308.05152v2).

Sections III.2–III.3, III.6–III.8, and Appendix D develop generalized and coset enumerators for coherent channels, logical error probabilities, effective-distance diagnostics, and decoding. Thus refinements of ordinary distance and algebraic evaluation of coherent logical channels are existing ideas. The present $r$ has a particular operational characterization: it is the minimum even logical-circuit weight and equals the optimal worst-case weak-angle exponent. The exact-support theorem identifies its infinite case. These statements are not the effective-distance definitions in Section III.3.

## Five-qubit code and stabilizer foundations

The five-qubit code is established work. R. Laflamme, C. Miquel, J. P. Paz, and W. H. Zurek introduced a perfect five-qubit code in *Perfect Quantum Error Correcting Code*, Physical Review Letters **77**, 198–201 (1996): [article](https://doi.org/10.1103/PhysRevLett.77.198), [arXiv](https://arxiv.org/abs/quant-ph/9602019). A five-qubit construction also appears in C. H. Bennett, D. P. DiVincenzo, J. A. Smolin, and W. K. Wootters, *Mixed-state entanglement and quantum error correction*, Physical Review A **54**, 3824–3851 (1996): [article](https://doi.org/10.1103/PhysRevA.54.3824), [arXiv v2](https://arxiv.org/abs/quant-ph/9604024v2).

For the explicit stabilizer convention used here, see Table 4 and the binary symplectic representation in Section 3.1 of D. Gottesman, *An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation*, Proceedings of Symposia in Applied Mathematics **68**, 13–58 (2010), [arXiv:0904.2557v1](https://arxiv.org/abs/0904.2557v1). Sections 2.2–2.3 develop the error-correction condition, stabilizers, syndrome measurement, and distance.

The family result in [code families](../science/code_families.md) proves the required parity and distance properties, and then analyzes the change caused by one physical CNOT. It does not attribute novelty to the underlying five-qubit code or to concatenation.

## Connection to protected logical generators

The [concatenated equality construction](../science/minimum_readout_safe_gates.md) uses the existing five-qubit code and canonical error-transparent base Hamiltonians together with this project's support criterion. At unbounded distance $d=3^m$, its recursive signed Pauli classes attain single-error-transparent interaction order $d+2$ and native coupling budget $(d+1)\omega$ for each logical Pauli axis while preserving $r_S=\infty$. The code and general stabilizer-dressing framework are inherited; the supplied gate construction and resource proof are recorded in the companion [Error-Transparent Locality account](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/minimum_readout_safe_gates.md). This is a mathematical coexistence statement for the declared gate and readout models, rather than a combined noise threshold or an exhaustive literature-priority claim.

## Sparse protected-gate descriptions

The [sparse gate account](../science/sparse_protected_gates.md) preserves the same infinite syndrome reversibility while replacing the dense protected-gate expansion with a linear number of native Pauli terms. Standard nonnegative support reduction gives at most $3n+1$ terms; a constructive cyclic-orbit recursion gives at most $4n$ with $O(\log n)$-bit coefficients. The linear-programming support principle is inherited (see [Hank Yang's background notes](https://hankyang.seas.harvard.edu/Semidefinite/background.html)); the exact count-vector reduction, recursive construction, and resource consequences are proved in the companion account. The code's original known-product-unitary readout assumptions remain part of the result.

## Two-error protected gates on the same reversible-readout code

The [two-error account](../science/two_error_gates.md) gives exact protected gates for all Pauli errors through weight two on the $[[25,1,9]]$ code while keeping $r_S=\infty$. Its order-thirteen interaction and native budget $521\omega/8$ attain code-specific lower bounds. The code and protected-Hamiltonian framework are inherited; the four-class equality construction and inner-block moment refinement are the derived statements supplied here. The sparse rational witness is checked independently of the numerical search that selected its support.

## Growing two-error gates with unchanged readout reversibility

The [all-depth two-error account](../science/growing_two_error_gates.md) extends the protected-gate equality to every concatenation depth $m\ge2$ without changing the code's infinite syndrome reversibility. The new construction uses correlated inactive-frontier padding; its resource lower bound uses squared bottom-block interaction weights. The concatenated code and the original support criterion remain inherited ingredients, and the gate and readout noise assumptions remain distinct.

## Reading the comparison

The cited results support the recovery framework, parity examples, coherent-error estimates, and code constructions. The additional conclusions are stated and proved in the topic pages linked at the top. A bibliography comparison establishes the relationship to these specified sources; it is not a guarantee that every related formulation in the literature has been located.
