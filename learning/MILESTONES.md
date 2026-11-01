# Thirty milestones for Syndrome Reversibility Distance

A prescribed error-syndrome measurement can reveal logical information after an otherwise reversible physical rotation. This map follows the objects and proof mechanisms that decide when that happens, quantify the unavoidable loss, and distinguish it from ordinary code distance. Its final part connects the readout invariant to protected logical Hamiltonians on the same code.

The map has **12 background/preliminaries milestones** and **18 narrative/innovation milestones**. The background joins physical Paulis, stabilizer measurements, logical channels and recovery into the complete experiment. The narrative joins exact-support orbits, logical circuits, recovery bounds and code constructions into the main readout results, then derives their consequences for protected control. The allocation gives most of the second phase to the readout classification and sharp exponent. The gate results enter as consequences and coexistence constructions, with their separate operational assumptions intact. The milestones supply global understanding, while the [scientific pages](../science/model.md) retain the complete specialist arguments.

The one common educational anchor for this project and [Error-Transparent Locality](https://github.com/GoGoKo699/Error-Transparent-Locality) is **John Preskill, *Quantum Information*, Chapter 7, “Quantum Error Correction,” updated March 2026**: [official chapter](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf), linked from the [author's course index](https://preskill.caltech.edu/ph229/index.html). Below, **P** refers to that chapter. The [Preskill bridge](PRESKILL_BRIDGE.md) translates its conventions and supplies the additional branch-recovery, channel-metric, Clifford/CNOT and protected-control steps locally. Research citations retain their attribution role and do not add educational anchors.

| Verified anchor location | Foundation used here |
|---|---|
| P §7.2, “Conditions for quantum error correction,” including Eq. (7.17) | Correctable error spans, recovery and degeneracy |
| P §§7.3.1–7.3.3 | Distinguishability, logical-information leakage and reference decoupling |
| P §7.4, especially §§7.4.1 and 7.4.3 | Pauli expansions, weight, distance and detection |
| P §§7.6–7.8 | Binary linear codes, CSS checks and encoded superpositions |
| P §§7.10.1–7.10.2 | Stabilizers, the real-Pauli convention and binary symplectic representation |
| P §7.10.4, “Encoded qubits” | Logical Pauli labels and the quotient by stabilizers |
| P §7.11, “The 5-Qubit Code,” including Eq. (7.150) | The seed code and syndrome measurement |
| P §7.17, “Concatenated codes” | Recursive encoding; the project proves its exact family parity and distance properties |

Each milestone gives a central question, an understanding target, earlier dependencies and its role in the story, a conceptual connection, and source pointers. Research-level recovery results and prior coherent-error mechanisms retain their separate attribution in [related work](../references/RELATED_WORK.md).

## Phase I: Background and preliminaries

### M01: Physical and logical quantum information

**Central question.** What is being preserved when several physical qubits encode fewer logical qubits?

**Understanding target.** Distinguish the physical tensor-product Hilbert space from the logical space and the codespace inside it. Explain the encoding isometry $V$, the projector $P=VV^\dagger$, and why restoration must work for arbitrary encoded superpositions, including states entangled with a reference.

**Dependencies and story role.** None. This fixes the input and output objects before introducing errors or measurements.

**Conceptual connection.** Encoding changes where information is represented. It does not change the requirement to preserve its quantum relationships, which will later make a classical syndrome record insufficient to undo a logical measurement.

**Source pointers.** P §§7.2–7.4 for code subspaces, correctability and reference systems; [Model §1](../science/model.md#1-code-phases-and-syndrome) and [Model §3](../science/model.md#3-performance-measures) for this project's spaces and reference convention.

### M02: Pauli labels, physical signs and support

**Central question.** Which features of a Pauli operator survive when its overall phase is ignored?

**Understanding target.** Identify support and weight, determine commutation from local anticommutations, and separate phase-free Pauli labels from signed Hermitian physical operators. Translate Preskill's real-Pauli convention into the Hermitian tensors used here. Explain why support calculations may discard global phases while syndrome projectors must retain the stabilizer signs.

**Dependencies and story role.** M01. This provides the local algebra used throughout both the criterion and its physical implementation.

**Conceptual connection.** The theorem compares exact sets of active qubits. Its recovery interpretation still depends on genuine observables, so combinatorial labels and physical operators must remain coordinated.

**Source pointers.** P §§7.4 and 7.10.1–7.10.2 for Pauli weight, the real-Pauli convention and symplectic labels; [Bridge §3](PRESKILL_BRIDGE.md#3-translate-the-notation-before-comparing-statements) translates Hermitian phases; [Model §1](../science/model.md#1-code-phases-and-syndrome) fixes physical signs.

### M03: A code as a common eigenspace

**Central question.** How do commuting Pauli checks specify an encoded quantum system?

**Understanding target.** Describe a stabilizer group that excludes $-I$, its independent generators, and its common positive-eigenvalue codespace. Relate $n-k$ independent checks to $2^k$ logical dimensions. Distinguish a generator list from the entire group it generates.

**Dependencies and story role.** M01, M02. The individual observables now define the global subspace on which logical information lives.

**Conceptual connection.** Products of generators can carry supports absent from the displayed generator list. That distinction becomes essential when the support criterion quantifies over the full stabilizer group.

**Source pointers.** P §7.10.1; [Model §1](../science/model.md#1-code-phases-and-syndrome) and the full-group qualification in [Support criterion §1](../science/support_criterion.md#1-the-support-criterion).

### M04: Syndrome sectors and error labels

**Central question.** What does the full syndrome measure, and what quantum state remains afterward?

**Understanding target.** Interpret syndrome bits as eigenvalue signs of the independent checks. Explain the orthogonal projectors $\Pi_s$, the equal dimension of their sectors, and the Pauli identification of each sector with the logical space. Keep the recorded outcome separate from the surviving quantum state.

**Dependencies and story role.** M02, M03. This makes the measurement operation concrete before considering its effect on logical information.

**Conceptual connection.** A Pauli error moves a code state to a labeled sector. A coherent superposition of errors can interfere within a sector, so the syndrome is not always a label for a single classical event.

**Source pointers.** P §§7.10.1–7.10.2 and §7.11 for checks and syndrome measurement; [Model §§1–2](../science/model.md) supplies the exact sector identifications and recorded instrument.

### M05: Logical operators and ordinary distance

**Central question.** Which undetected Pauli actions change the encoded information?

**Understanding target.** Distinguish detectable Paulis, stabilizers, and nontrivial logical Paulis. Identify the centralizer $N$, the logical quotient $N/S$, and the ordinary distance $d$ as the minimum weight of a nontrivial logical label. Explain why changing only a global phase does not produce a logical operation.

**Dependencies and story role.** M03, M04. This separates harmless and harmful zero-syndrome actions and establishes the conventional parameter to which $r$ will be compared.

**Conceptual connection.** Ordinary distance concerns the existence of a logical Pauli of a given weight. The later invariant will add a condition about how the fixed measurement can expose that operator.

**Source pointers.** P §§7.4.1, 7.4.3 and 7.10.4; [Model §1](../science/model.md#1-code-phases-and-syndrome) and [Support criterion §5](../science/support_criterion.md#5-relation-to-ordinary-distance).

### M06: Coherent rotations as Pauli amplitudes

**Central question.** How does a reversible local rotation produce several possible syndrome paths?

**Understanding target.** Expand a single-qubit Pauli-axis rotation into identity and Pauli terms, then form the product expansion across sites. Track complex amplitudes before measurement and probabilities afterward. Distinguish one specified coherent unitary from an unknown random mixture of unitaries.

**Dependencies and story role.** M02, M04. This joins the physical noise model to the syndrome structure.

**Conceptual connection.** Two Pauli terms with the same syndrome remain coherent inside a branch. Their relative phase will determine whether that branch acts as a reversible logical rotation or a weak logical measurement.

**Source pointers.** P §§7.2 and 7.4 for error-operator expansions and coherent errors; [Model §2](../science/model.md#2-the-prescribed-order-of-operations), [Bridge §§3–5](PRESKILL_BRIDGE.md#3-translate-the-notation-before-comparing-statements) and [Error exponent §3](../science/error_exponent.md#3-an-embedded-even-repetition-witness) develop the specified rotations and their interference.

### M07: The recorded measurement instrument

**Central question.** How is a full syndrome measurement represented as a channel on the logical input?

**Understanding target.** Construct the channel that retains a classical syndrome and the corresponding physical output. Interpret its logical branch map $K_s$ and effect $M_s=K_s^\dagger K_s$. Explain why the probability of outcome $s$ is $\mathrm{Tr}(M_s\rho)$ and why the effects sum to the identity.

**Dependencies and story role.** M01, M04, M06. This organizes the physical process into the objects optimized in the theorems.

**Conceptual connection.** A physical readout induces a measurement on the encoded input. The effect operators reveal exactly which input information the outcome probabilities can depend on.

**Source pointers.** P §§7.2–7.3 and 7.10.1 provide operations, information-leakage and syndrome foundations; the explicit instrument is developed in [Model §2](../science/model.md#2-the-prescribed-order-of-operations) and [Bridge §4](PRESKILL_BRIDGE.md#4-the-syndrome-induces-a-measurement-on-the-logical-state).

### M08: Scalar effects and reversible branches

**Central question.** When can knowing a measurement outcome allow perfect restoration?

**Understanding target.** Explain why $M_s=p_sI$ means that the outcome reveals no logical-input information. For this one-Kraus-per-outcome instrument, connect scalar effects to branches proportional to isometries and to exact recovery. Recognize why the same scalar test would be insufficient for a general branch with unobserved Kraus multiplicity.

**Dependencies and story role.** M05, M07. This gives the operational recovery condition that the new code criterion will classify.

**Conceptual connection.** Preserving all logical superpositions is equivalent here to preventing the classical readout from resolving their amplitudes. The condition is an inherited recovery fact, with its source attribution retained.

**Source pointers.** P §§7.2–7.3 for correction conditions and absence of logical-information leakage; [Model §5](../science/model.md#5-exact-reversibility-and-logical-information-in-the-syndrome) gives the precise one-Kraus branch theorem and its research attribution.

### M09: Polar recovery and the fidelity optimum

**Central question.** Which part of a known branch map can recovery undo?

**Understanding target.** Separate a branch into its polar unitary and positive factor. Explain why undoing the unitary leaves $\sqrt{M_s}$, and how its singular values determine the optimal entanglement fidelity. Include singular and zero branches without assuming an inverse of $K_s$ exists.

**Dependencies and story role.** M07, M08. Exact reversibility becomes a quantitative benchmark for imperfect recovery.

**Conceptual connection.** Rotations of a branch can be corrected with full knowledge of the unitary. Unequal singular values encode the residual distortion that determines the recovery optimum.

**Source pointers.** P §7.2 supplies correction context. The polar optimization is an additional local derivation of the attributed Gregoratti–Werner result, supplied in [Model §4](../science/model.md#4-exact-optimal-entanglement-fidelity) and [Bridge §6](PRESKILL_BRIDGE.md#6-recovery-and-error-measures-supplied-by-the-project); it is not assigned as a theorem in Preskill.

### M10: Errors that include a reference system

**Central question.** How should recovery error be measured when the logical system may be entangled with something else?

**Understanding target.** Distinguish entanglement infidelity from half-diamond distance. Explain the optimization over recovery for each metric and the inequality $e_{\mathrm{opt}}\leq\varepsilon_{\mathrm{opt}}$. Understand that an entanglement-fidelity optimizer need not optimize diamond distance for a general unitary.

**Dependencies and story role.** M07, M09. This specifies the strength of the eventual upper and lower bounds.

**Conceptual connection.** A lower bound derived from a maximally entangled input already constrains every recovery in the stronger channel metric. An explicit recovered channel will later attain that bound.

**Source pointers.** P §§7.2 and 7.3.3 supply correctability and reference-system context. The optimized entanglement-fidelity and diamond-distance comparison is developed locally in [Model §3](../science/model.md#3-performance-measures) and [Bridge §6](PRESKILL_BRIDGE.md#6-recovery-and-error-measures-supplied-by-the-project), beyond the chapter's treatment.

### M11: Clifford changes and physical locality

**Central question.** Which changes of Pauli basis preserve the physical support structure?

**Understanding target.** Follow Pauli operators under Clifford conjugation. Distinguish local Clifford axis changes and qubit permutations from an entangling CNOT. Explain why an entangling transformation of a code can change its relation to product-unitary noise even though the transformed subspace is unitarily equivalent.

**Dependencies and story role.** M02, M05, M06. This prepares both the local-Clifford witness and the later comparison of two code families.

**Conceptual connection.** The noise class refers to the physical division into qubits. Support-sensitive properties therefore contain information that abstract equivalence of code spaces alone does not preserve.

**Source pointers.** P §§7.10–7.11 provide stabilizer and encoding context. [Bridge §3](PRESKILL_BRIDGE.md#local-derivation-of-the-clifford-and-cnot-rules) derives the required local-Clifford and CNOT rules; [Support criterion §§2 and 5](../science/support_criterion.md) and [Interpretation: the CNOT comparison](../science/interpretation.md#why-the-cnot-comparison-is-consequential) give their project roles.

### M12: The complete operational question

**Central question.** What is optimized, what is fixed, and what must the code preserve?

**Understanding target.** State the full experiment: encode, apply a known product unitary, measure the original complete syndrome ideally, then recover. Explain the quantifier “for every specified unitary, there exists a recovery,” the retained classical record and quantum output, and the unrestricted resources allowed to recovery.

**Dependencies and story role.** M03–M11. This completes the first phase by assembling its local ingredients into the precise global task.

**Conceptual connection.** A prescribed readout can destroy information that a known unitary alone would preserve. The project measures the cost of this fixed order rather than assuming that the inverse rotation can be applied before measurement.

**Source pointers.** P §§7.2–7.4 and 7.10 for the QEC framework; [Model §§2 and 6](../science/model.md) and [Interpretation: operational boundaries](../science/interpretation.md#operational-boundaries) fix the project's exact ordering and access assumptions.

## Phase II: Narrative and innovation

### M13: From the recorded syndrome to exact-support orbits

**Central question.** How can a condition on every measurement branch become a condition on physical Pauli supports?

**Understanding target.** Use Fourier inversion between the full syndrome projectors and the full stabilizer group to replace scalarity of every $M_s$ by scalarity of every $V^\dagger U^\dagger gUV$. Explain why product conjugation spans precisely the tensor Paulis with the same exact support as $g$.

**Dependencies and story role.** M02–M08, M11, M12. This first narrative step connects the operational recovery test to a local algebraic mechanism.

**Conceptual connection.** The complete measurement contains information from every stabilizer product. Local rotations can vary the axes on active sites while preserving the set of active sites.

**Source pointers.** P §§7.10.1–7.10.2 supplies stabilizer and symplectic foundations; [Bridge §3](PRESKILL_BRIDGE.md#local-derivation-of-the-clifford-and-cnot-rules) gives the additional conjugation rules. [Support criterion §1](../science/support_criterion.md#1-the-support-criterion) and [Bridge §7](PRESKILL_BRIDGE.md#7-why-exact-physical-support-classifies-reversibility) prove the orbit argument.

### M14: The exact classification and its measurement witness

**Central question.** Which codes admit exact restoration after this readout for every specified product unitary?

**Understanding target.** Combine the support orbit with the three compression cases: a detectable Pauli gives zero, a stabilizer gives a scalar, and a nontrivial logical Pauli gives a nonscalar logical action. State the if-and-only-if absence of equal exact support between the full stabilizer group and any nontrivial logical Pauli. Explain the product-Clifford witness and its bounds $F_{\mathrm{opt}}\le1/2$ and $\varepsilon_{\mathrm{opt}}\ge1/2$ when a collision exists.

**Dependencies and story role.** M05, M08–M10, M13. The local mechanism now yields the complete exact-reversibility theorem.

**Conceptual connection.** A matching support lets local basis changes turn one measured stabilizer into a logical observable. Recording its eigenvalue destroys coherence between the logical eigenspaces even when the rotation is known.

**Source pointers.** P §§7.2–7.3 and 7.10.4 supplies correction and encoded-operator context; [Support criterion §§1–2](../science/support_criterion.md) proves the classification and finite-angle witness. [Related work](../references/RELATED_WORK.md) separates the classification from inherited branch-recovery and logical-projection facts.

### M15: The distance that counts changed local axes

**Central question.** What quantitative invariant is hidden in an exact-support collision?

**Understanding target.** Define $r(S)=\min\mathrm{wt}(gL)$ over pairs $g\in S$, $L\in N\setminus S$ with equal exact support, with an empty minimum equal to infinity. Translate it into a logical Pauli $E$ that anticommutes locally with one stabilizer at every active site. Derive that finite $r$ is even and at least the ordinary distance $d$.

**Dependencies and story role.** M02, M05, M14. This turns an existence classification into a candidate measure of weak-rotation sensitivity.

**Conceptual connection.** The weight being minimized counts differing local axes within a matching-support pair. It is neither their common support size nor the minimum weight of an arbitrary logical operator.

**Source pointers.** P §§7.4 and 7.10.1–7.10.4 supplies weight, Pauli commutation and logical classes; [Support criterion §3](../science/support_criterion.md#3-syndrome-reversibility-distance) and [Bridge §8](PRESKILL_BRIDGE.md#8-from-a-support-collision-to-a-logical-circuit) define and interpret $r$.

### M16: Minimal syndrome dependencies and even logical circuits

**Central question.** How does a minimal dependence among local syndrome columns characterize $r$?

**Understanding target.** Define a logical circuit using a fixed Pauli word and the absence of a proper nonempty normalizer subword. Explain why a minimum local-anticommutation witness is a circuit and why an even circuit supplies a stabilizer witness by binary row-space/kernel duality. Conclude that $r$ is the minimum even logical-circuit weight and that even ordinary distance forces $r=d$.

**Dependencies and story role.** M04, M05, M15. This identifies the minimal interference structure that the quantitative proof will use.

**Conceptual connection.** A circuit means a minimal binary dependence here, rather than a sequence of quantum gates. Global even weight alone does not guarantee the required minimality or a readout-loss witness.

**Source pointers.** P §7.10.2 supplies binary symplectic foundations; [Support criterion §§4–5](../science/support_criterion.md#4-equivalent-characterization-by-even-logical-circuits) proves the circuit theorem and distance comparison; [Bridge §8](PRESKILL_BRIDGE.md#8-from-a-support-collision-to-a-logical-circuit) develops the additional row-space/kernel connection.

### M17: A uniform upper bound without discarding rare branches

**Central question.** Why does an outcome with very small probability not spoil the weak-angle exponent?

**Understanding target.** Form the collision polynomial $B_S(x)=\sum_{(g,L)}x^{\mathrm{wt}(gL)}$. For bounded local generators, relate the sum of nonscalar branch effects to $B_S(|\theta|)$. Explain how polar recovery, the double-commutator channel identity and a square-root norm estimate bound the error without division by a branch probability. Obtain $e_{\mathrm{opt}}\le\varepsilon_{\mathrm{opt}}\le\min\{1,B_S(|\theta|)\}$ for $|\theta|\le1$.

**Dependencies and story role.** M06, M09, M10, M13, M15. This establishes the uniform half of the quantitative theorem before choosing any adverse rotation.

**Conceptual connection.** A small effect may be highly state dependent after normalization. Bounding the unnormalized recovered channel keeps rare outcomes at their correct total weight and preserves the exponent.

**Source pointers.** P §§7.2 and 7.4 motivates coherent Pauli expansions. The recovery and channel-norm proof is supplied locally in [Error exponent §2](../science/error_exponent.md#2-a-finite-code-collision-polynomial) and [Bridge §9](PRESKILL_BRIDGE.md#9-why-rare-syndrome-outcomes-do-not-weaken-the-upper-exponent), beyond the anchor.

### M18: Complementary paths give an all-recovery converse

**Central question.** Which physical rotation forces every recovery to lose the same leading-order information?

**Understanding target.** Rotate about the fixed labels of a minimum even circuit. Show that only complementary Pauli subwords share a syndrome, and reduce each branch to $aI\pm b\overline E$. Use the two singular values and the exact fidelity optimum to identify the smaller complementary probability as irrecoverable. Explain why this constrains arbitrary quantum recovery rather than only Pauli decoding.

**Dependencies and story role.** M06, M09, M10, M16. This supplies the lower-bound mechanism matching the uniform upper bound.

**Conceptual connection.** Even circuit length makes complementary amplitudes interfere with a real relative phase. The resulting singular-value imbalance is a logical measurement, and cannot be removed by reversing a branch unitary.

**Source pointers.** P §§7.2–7.3 supplies coherent-error and no-information-leakage context; [Error exponent §3](../science/error_exponent.md#3-an-embedded-even-repetition-witness), [Bridge §10](PRESKILL_BRIDGE.md#10-why-the-exponent-is-attained-against-every-recovery) and [related work](../references/RELATED_WORK.md) develop and attribute the repetition mechanism and code-wide converse.

### M19: The exact tail and sharp fixed-code error law

**Central question.** How do the upper and lower mechanisms combine into the operational meaning of $r$?

**Understanding target.** Interpret the exact even-repetition tail $p_r(q)$, including the half-weighted tie at $r/2$, with $q=\sin^2(\theta/2)\le1/2$. Explain why a recovered dephasing channel makes both optimized errors equal to that tail for the witness. State the worst-case laws $\mathcal E_S(\theta),\mathcal D_S(\theta)=\Theta(|\theta|^r)$ for finite $r$ and exact zero for $r=\infty$, with the code fixed and recovery optimized separately for each known rotation.

**Dependencies and story role.** M10, M12, M17, M18. This assembles the quantitative result as a theorem with its complete quantifiers.

**Conceptual connection.** The metrics coincide for the witness while remaining different optimization problems in general. Code-dependent constants and a fixed-code limit do not establish a uniform threshold for a growing family.

**Source pointers.** P §§7.2–7.3 supplies the arbitrary-input and reference-preserving correction setting; [Error exponent §§1, 3–5](../science/error_exponent.md) and [Bridge §§9–10](PRESKILL_BRIDGE.md#9-why-rare-syndrome-outcomes-do-not-weaken-the-upper-exponent) develop the metric argument, exact tail and fixed-code quantifiers.

### M20: The five-qubit seed as a complete parity statement

**Central question.** Which properties of the standard five-qubit code will protect it against this measurement loss?

**Understanding target.** Relate the four displayed stabilizer generators to the complete stabilizer group, single-qubit syndromes and minimum logical representatives. Establish that every stabilizer has even weight and every representative of every nontrivial logical Pauli has odd weight. Keep this full-coset statement separate from a check of selected generators.

**Dependencies and story role.** M03–M05, M14. This begins the construction that separates the new invariant from ordinary code parameters.

**Conceptual connection.** A parity difference excludes exact-support equality. The construction uses the familiar five-qubit code as a verified ingredient in a new support-sensitive comparison.

**Source pointers.** P §§7.10.2 and 7.11, especially Eq. (7.150), supplies the seed and binary language; [Code families §2](../science/code_families.md#2-the-five-qubit-building-block) proves the exact distance and full-coset parity facts; [Bridge §11](PRESKILL_BRIDGE.md#11-from-the-five-qubit-code-to-the-one-cnot-separation) joins them to the support theorem.

### M21: Concatenation and blockwise exact reversibility

**Central question.** Why do the seed's distance and support separation persist at arbitrary depth and on two blocks?

**Understanding target.** Explain the recursive family $Q_m=[[5^m,1,3^m]]$, the matching distance lower and upper bounds, and the parity induction for full stabilizer and logical cosets. On $Q_m\otimes Q_m$, use parity inside each encoded block to exclude every stabilizer/logical support collision, including logical operators of even total weight. Conclude $r=\infty$.

**Dependencies and story role.** M14, M20. A local seed now gives an unbounded-distance family with exact readout recovery.

**Conceptual connection.** The decisive partition is into physical encoded blocks. A logical acting on both blocks may be globally even, while its odd restriction in a nontrivial block still forbids a matching stabilizer support.

**Source pointers.** P §7.17 supplies recursive-code context; [Code families §§3–4](../science/code_families.md#3-concatenation-preserves-the-required-properties) proves this family's exact distance and blockwise readout property, beyond the chapter's general concatenation bound.

### M22: One physical CNOT changes the code without changing distance

**Central question.** How can an entangling Clifford alter readout reversibility while leaving $[[n,k,d]]$ unchanged?

**Understanding target.** Conjugate both the product codespace and its full stabilizer group by one physical CNOT between blocks. Follow a minimum logical that avoids the control to preserve an upper bound on distance. Use inverse conjugation and the minimum-weight product-code structure to exclude a shorter logical. Explain why each resulting code is measured with its own original checks.

**Dependencies and story role.** M11, M20, M21. This introduces the physical transformation used to compare two codes with equal conventional parameters.

**Conceptual connection.** Abstract unitary equivalence does not preserve the tensor-product structure of the allowed rotations. The CNOT belongs to the definition of the second code; the readout experiment still uses product-unitary noise.

**Source pointers.** P §§7.10–7.11 supplies stabilizer and encoding context; [Bridge §3](PRESKILL_BRIDGE.md#local-derivation-of-the-clifford-and-cnot-rules) derives CNOT conjugation, and [Code families §5](../science/code_families.md#5-one-cnot-preserves-the-distance) and [Bridge §11](PRESKILL_BRIDGE.md#11-from-the-five-qubit-code-to-the-one-cnot-separation) prove distance preservation.

### M23: The new even circuit and the same-parameter separation

**Central question.** What readout-loss structure does that single physical CNOT introduce?

**Understanding target.** Construct the transformed logical of weight $3^m+1$ and explain why no proper nonempty fixed-label subword normalizes the transformed stabilizer. Combine its even circuit with parity and distance to obtain $r=3^m+1$. Assemble the pair of $[[2\cdot5^m,2,3^m]]$ families with $r=\infty$ and $r=3^m+1$.

**Dependencies and story role.** M15, M16, M19, M21, M22. This completes the growing-distance separation and attaches its sharp weak-angle error law.

**Conceptual connection.** Ordinary distance fixes the smallest undetected Pauli weight. It does not fix whether a prescribed syndrome measurement can expose a coherent logical observable after local rotations.

**Source pointers.** P §§7.10.4, 7.11 and 7.17 supplies logical and concatenation foundations, with CNOT rules derived in [Bridge §3](PRESKILL_BRIDGE.md#local-derivation-of-the-clifford-and-cnot-rules). [Code families §§1 and 6](../science/code_families.md) and [Error exponent §1](../science/error_exponent.md#1-statement-and-quantifiers) give the project theorem and its quantitative consequence.

### M24: A physical ten-qubit channel checks the construction

**Central question.** How does the smallest pair realize the abstract circuit mechanism in the full physical channel?

**Understanding target.** At $m=1$, connect the weight-four circuit to the exact error $p_4(q)=3q^2-2q^3$. Explain what constructing all signed syndrome sectors and the recovered channel checks beyond a support enumeration. Distinguish this finite certificate from the analytical all-depth proof and from a check of arbitrary noisy extraction circuits.

**Dependencies and story role.** M07, M09, M18, M19, M23. This makes the core theorem concrete before connecting it to a second physical task.

**Conceptual connection.** A support witness, a physical projective measurement and a logical recovered channel are three representations of the same experiment. Their agreement checks signs and phases that a phase-free support calculation omits.

**Source pointers.** P §§7.10–7.11 supplies physical syndrome foundations; [Code families §7](../science/code_families.md#7-the-ten-qubit-example-and-its-physical-channel) and [Verification](../verification/README.md) specify the finite certificate, while [Error exponent §3](../science/error_exponent.md#3-an-embedded-even-repetition-witness) proves its all-recovery meaning.

### M25: Two operational experiments on one stabilizer code

**Central question.** What additional assumptions are needed to relate the readout invariant to protected logical gates?

**Understanding target.** State the error-transparent Hamiltonian identity $HEV=EV(h_L+\lambda_E I)$ for every Pauli of weight at most $t$, allowing consistent scalar error-sector shifts. Distinguish it from strict transparency with zero shifts. Define the maximum native Pauli interaction order $\ell$ and coefficient budget $J$. Keep faults during protected evolution and ideal final recovery separate from the known-product-unitary, prescribed-readout experiment.

**Dependencies and story role.** M01, M05, M12, M19, M23. This supplies the operator framework needed for consequences involving the same code, without changing the main readout model.

**Conceptual connection.** A code can be evaluated for two tasks. An operator identity protecting an error image during a gate and an information-preserving measurement condition after a rotation are compatible properties with distinct hypotheses.

**Source pointers.** P §§7.2–7.4 and 7.10.4 supplies correction and encoded-operator foundations; [Gate/readout theorem §1](../science/gate_readout_tradeoff.md#1-notation-and-hypotheses) and [Bridge §12](PRESKILL_BRIDGE.md#12-two-tasks-on-the-same-code) develop the additional Hamiltonian model and native resources locally.

### M26: Minimum gate locality forces a readout witness

**Central question.** What does attaining the general protected-gate interaction bound imply for readout reversibility?

**Understanding target.** Under odd $d\ge3$, $t\ge1$ and $2t+1\le d$, explain the implication $\ell=d+t\Rightarrow r_S=d+1$. Follow the weight-one-stabilizer reduction, operator weight moments, forced logical term of weight $(d+1)$, and even-circuit argument. Deduce $\ell\ge d+t+1$ whenever $r_S>d+1$, including infinity, without treating this necessary condition as an existence result.

**Dependencies and story role.** M15, M16, M25. The readout characterization now constrains an independently defined gate task.

**Conceptual connection.** Saturating the gate's moment constraints forces an even logical representative that distance makes minimal after fixed qubits are removed. That representative is precisely a circuit capable of exposing logical information to readout.

**Source pointers.** P §§7.10.1–7.10.4 supplies stabilizer and logical foundations. The operator moments and fixed-qubit proof are supplied in [Gate/readout theorem §2](../science/gate_readout_tradeoff.md#2-joint-theorem) and [Bridge §13](PRESKILL_BRIDGE.md#13-why-minimum-gate-order-forces-readout-loss), beyond the educational anchor.

### M27: Optimal one-, two- and three-error gates coexist with $r_S=\infty$

**Central question.** How can the readout-safe concatenated family support exact protected gates at its own minimum interaction order?

**Understanding target.** Explain why its odd logical-coset weights strengthen the gate bound to $\ell\ge d+2t$. Identify the attained one-axis orders and budgets: $d+2,d+1$ for $t=1$; $d+4,J_2(d)$ for $t=2$; and $d+6,J_3(d)$ for $t=3$, with $m\ge1$ for $t=1$ and $m\ge2$ for $t=2,3$. Connect strict lower-order increments and their vanishing action on lighter error images to recursive lifting. Interpret coupling optimality at the stated minimum order.

**Dependencies and story role.** M20, M21, M25, M26. This gives constructive coexistence results on the readout-safe family rather than only an obstruction to equality in the general bound.

**Conceptual connection.** Parity both protects readout and removes candidate interaction weights. Exact gate constructions attain the stronger family-specific bound, while native term count and coefficient size remain separate resources.

**Source pointers.** P §§7.10–7.11 and 7.17 supplies code and concatenation context; [One-error gates](../science/minimum_readout_safe_gates.md), [Two-error gates](../science/growing_two_error_gates.md), [Three-error hierarchy](../science/growing_three_error_gates.md) and [Bridge §14](PRESKILL_BRIDGE.md#14-protected-control-with-the-original-readout-unchanged) develop the additional constructions, exact $J_2,J_3$ formulas and lower bounds.

### M28: Protected entanglers preserve the product code and its checks

**Central question.** Can encoded entanglement coexist with exact product-rotation readout recovery on the unchanged code?

**Understanding target.** On unchanged $Q_m$ product blocks with $m\ge2$, take a nonidentity pure logical Pauli string of logical weight $p\ge1$. Distinguish its minimum physical coset weight $D=pd$ from the product code's ordinary distance $d$. Explain the total-degree-three tensor-increment construction and its optimum $\ell=D+6$, $J=|\Omega|J_3(D)$. Relate the fixed-code $ZZ$ entangler and local logical axes to universal two-qubit control while retaining $r_S=\infty$ block by block.

**Dependencies and story role.** M11, M21, M25, M27. The coexistence statement extends from separate encoded axes to entangling encoded dynamics.

**Conceptual connection.** A protected logical entangler preserves the codespace and stabilizer measurement. The earlier physical CNOT conjugated the code and checks. Their different actions on the fixed physical code explain their different readout consequences.

**Source pointers.** P §§7.10.4 and 7.17 supplies encoded Pauli and product-encoding foundations. [Bridge §3](PRESKILL_BRIDGE.md#local-derivation-of-the-clifford-and-cnot-rules) gives basic gate rules; [Protected entangling gates §§1–7](../science/protected_entangling_gates.md) and [Bridge §14](PRESKILL_BRIDGE.md#14-protected-control-with-the-original-readout-unchanged) prove the additional tensor construction, target-coset lower bound and control argument.

### M29: What exact certificates and resource optima establish

**Central question.** Which parts of the combined account are analytical theorems, finite certificates, and native-model resource statements?

**Understanding target.** Distinguish the general support/exponent proofs, all-depth code and gate recursions, finite physical-channel checks, exact rational class certificates and profile-based error coverage. Explain why compact class descriptions do not make the dense native expansion small. Restrict the pure-string budget optimum to its specified target and minimum order; do not infer an optimum for an arbitrary mixture of logical Pauli terms at one common cap.

**Dependencies and story role.** M17–M19, M23, M24, M27, M28. This fixes the evidential and resource meaning of every ingredient before the final synthesis.

**Conceptual connection.** Verification can check all error profiles in a finite construction without enumerating its native interaction list. Neither exact arithmetic nor Hamiltonian locality establishes hardware realizability or a joint noisy gate/readout threshold.

**Source pointers.** P §§7.2–7.4 and 7.17 supplies correction and concatenation context, while the current resource and verification claims have their own hypotheses. [Verification](../verification/README.md), [Three-error examples](../science/growing_three_error_gates.md#5-examples-description-size-and-verification), [Entangler examples and scope](../science/protected_entangling_gates.md#5-the-first-protected-entangler) and [Bridge §§14–15](PRESKILL_BRIDGE.md#14-protected-control-with-the-original-readout-unchanged) give those exact local statements.

### M30: The complete scientific account

**Central question.** What new information about a stabilizer code is obtained by combining the readout theorem, sharp error law, same-parameter family and gate consequences?

**Understanding target.** Connect logical effects, exact-support separation, even circuits, all-recovery weak-angle loss and the physical-CNOT family into the main readout result. Then explain the joint locality obstruction and the protected one-, two-, three-error and entangling coexistence constructions on unchanged readout-safe codes. State which ingredients are inherited, which are proved here, and which assumptions belong to each operational task.

**Dependencies and story role.** M12–M29. This completes the narrative by integrating its local mechanisms into their global meaning, significance and limits.

**Conceptual connection.** Conventional code parameters leave out how the original syndrome measurement responds to coherent local rotations. The new invariant captures that response, constrains equality in a gate-resource theorem, and still permits exact protected encoded control under its separate model.

**Source pointers.** P §§7.2–7.4, 7.10–7.11 and 7.17 supplies the common correction, stabilizer and concatenation foundation. [Interpretation](../science/interpretation.md), [Gate/readout theorem](../science/gate_readout_tradeoff.md), [Protected entanglers](../science/protected_entangling_gates.md), [Related work](../references/RELATED_WORK.md) and [Bridge §15](PRESKILL_BRIDGE.md#15-assemble-the-complete-statement) give the authoritative synthesis, attribution and separate operational boundaries.
