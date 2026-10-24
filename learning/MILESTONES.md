# Thirty milestones for Syndrome Reversibility Distance

This map follows the scientific structure of the project. It is a common route for any reader: individual Pauli operators and measurements lead to the precise recovery question, then local proof mechanisms join into the classification, sharp error law, and code-family separation. The complete proofs remain in the authoritative [scientific pages](../README.md#read-the-science).

The map has **12 background/preliminaries milestones** and **18 narrative/innovation milestones**. The first phase builds the operational framework from its component objects. The second phase starts with the support mechanism, establishes its quantitative consequence, and ends with the growing-distance construction and the meaning of the combined result. Inherited ingredients, such as polar recovery and the five-qubit code, appear where the argument needs them; placement in the second phase does not claim that those ingredients are new.

The single educational anchor is **Daniel Gottesman, *An Introduction to Quantum Error Correction and Fault-Tolerant Quantum Computation***, [arXiv:0904.2557v1](https://arxiv.org/abs/0904.2557v1), the 2009 preprint of the chapter published in 2010. Below, **G** denotes this version. Use the [Gottesman bridge](GOTTESMAN_BRIDGE.md) for notation, explanations and the steps supplied by this project. The bridge develops the channel metrics and optimized branch recovery beyond the anchor's treatment; it does not attribute those results to Gottesman.

| Anchor location | Role in this route |
|---|---|
| G §2.1 | Pauli operators, coherent error expansions, and linearity of correction |
| G §2.2 | General error-correction conditions and code distance |
| G §2.3 | Stabilizer checks, syndrome measurement, and the codespace |
| G §3.1, including Table 4 and Eq. (56) | Five-qubit construction and binary symplectic representation |
| G §3.4 | Logical Pauli operators and the quotient by stabilizers |
| G §3.5 | Clifford conjugation, including local changes of axis and CNOT |
| G §5.2, concatenation discussion | Recursive encoding; the project's family proof supplies its specific distance and parity claims |

Each milestone states a central question, an understanding target, earlier dependencies and its role in the story, a conceptual connection, and source pointers. These targets cover the connected argument without replacing every specialist calculation.

## Phase I — Background and preliminaries

### M01 — Physical and logical quantum information

**Central question.** What is being preserved when several physical qubits encode fewer logical qubits?

**Understanding target.** Distinguish the physical tensor-product Hilbert space from the logical space and the codespace inside it. Explain the encoding isometry $V$, the projector $P=VV^\dagger$, and why restoration must work for arbitrary encoded superpositions, including states entangled with a reference.

**Dependencies and story role.** None. This fixes the input and output objects before introducing errors or measurements.

**Conceptual connection.** Encoding changes where information is represented. It does not change the requirement to preserve its quantum relationships, which will later make a classical syndrome record insufficient to undo a logical measurement.

**Source pointers.** G §§2.1–2.2 for code subspaces and encoded states; [Model §1](../science/model.md#1-code-phases-and-syndrome) and [Model §3](../science/model.md#3-performance-measures) for this project's spaces and reference convention.

### M02 — Pauli labels, physical signs and support

**Central question.** Which features of a Pauli operator survive when its overall phase is ignored?

**Understanding target.** Identify support and weight, determine commutation from local anticommutations, and separate phase-free Pauli labels from signed Hermitian physical operators. Explain why support calculations may discard global phases while syndrome projectors must retain the stabilizer signs.

**Dependencies and story role.** M01. This provides the local algebra used throughout both the criterion and its physical implementation.

**Conceptual connection.** The theorem compares exact sets of active qubits. Its recovery interpretation still depends on genuine observables, so combinatorial labels and physical operators must remain coordinated.

**Source pointers.** G §2.1, Table 1, and §3.1 for the binary representation; [Model §1](../science/model.md#1-code-phases-and-syndrome).

### M03 — A code as a common eigenspace

**Central question.** How do commuting Pauli checks specify an encoded quantum system?

**Understanding target.** Describe a stabilizer group that excludes $-I$, its independent generators, and its common positive-eigenvalue codespace. Relate $n-k$ independent checks to $2^k$ logical dimensions. Distinguish a generator list from the entire group it generates.

**Dependencies and story role.** M01, M02. The individual observables now define the global subspace on which logical information lives.

**Conceptual connection.** Products of generators can carry supports absent from the displayed generator list. That distinction becomes essential when the support criterion quantifies over the full stabilizer group.

**Source pointers.** G §2.3; [Model §1](../science/model.md#1-code-phases-and-syndrome) and the full-group qualification in [Support criterion §1](../science/support_criterion.md#1-the-support-criterion).

### M04 — Syndrome sectors and error labels

**Central question.** What does the full syndrome measure, and what quantum state remains afterward?

**Understanding target.** Interpret syndrome bits as eigenvalue signs of the independent checks. Explain the orthogonal projectors $\Pi_s$, the equal dimension of their sectors, and the Pauli identification of each sector with the logical space. Keep the recorded outcome separate from the surviving quantum state.

**Dependencies and story role.** M02, M03. This makes the measurement operation concrete before considering its effect on logical information.

**Conceptual connection.** A Pauli error moves a code state to a labeled sector. A coherent superposition of errors can interfere within a sector, so the syndrome is not always a label for a single classical event.

**Source pointers.** G §2.3 and §3.1, Lemma 1; [Model §§1–2](../science/model.md).

### M05 — Logical operators and ordinary distance

**Central question.** Which undetected Pauli actions change the encoded information?

**Understanding target.** Distinguish detectable Paulis, stabilizers, and nontrivial logical Paulis. Identify the centralizer $N$, the logical quotient $N/S$, and the ordinary distance $d$ as the minimum weight of a nontrivial logical label. Explain why changing only a global phase does not produce a logical operation.

**Dependencies and story role.** M03, M04. This separates harmless and harmful zero-syndrome actions and establishes the conventional parameter to which $r$ will be compared.

**Conceptual connection.** Ordinary distance concerns the existence of a logical Pauli of a given weight. The later invariant will add a condition about how the fixed measurement can expose that operator.

**Source pointers.** G §§2.2–2.3 and §3.4; [Model §1](../science/model.md#1-code-phases-and-syndrome) and [Support criterion §5](../science/support_criterion.md#5-relation-to-ordinary-distance).

### M06 — Coherent rotations as Pauli amplitudes

**Central question.** How does a reversible local rotation produce several possible syndrome paths?

**Understanding target.** Expand a single-qubit Pauli-axis rotation into identity and Pauli terms, then form the product expansion across sites. Track complex amplitudes before measurement and probabilities afterward. Distinguish one specified coherent unitary from an unknown random mixture of unitaries.

**Dependencies and story role.** M02, M04. This joins the physical noise model to the syndrome structure.

**Conceptual connection.** Two Pauli terms with the same syndrome remain coherent inside a branch. Their relative phase will determine whether that branch acts as a reversible logical rotation or a weak logical measurement.

**Source pointers.** G §2.1, especially the rotation and error-expansion discussion; [Model §2](../science/model.md#2-the-prescribed-order-of-operations) and [Error exponent §3](../science/error_exponent.md#3-an-embedded-even-repetition-witness).

### M07 — The recorded measurement instrument

**Central question.** How is a full syndrome measurement represented as a channel on the logical input?

**Understanding target.** Construct the channel that retains a classical syndrome and the corresponding physical output. Interpret its logical branch map $K_s$ and effect $M_s=K_s^\dagger K_s$. Explain why the probability of outcome $s$ is $\mathrm{Tr}(M_s\rho)$ and why the effects sum to the identity.

**Dependencies and story role.** M01, M04, M06. This organizes the physical process into the objects optimized in the theorems.

**Conceptual connection.** A physical readout induces a measurement on the encoded input. The effect operators reveal exactly which input information the outcome probabilities can depend on.

**Source pointers.** G §§2.1–2.3 for operations and syndrome measurement; the explicit instrument is developed in [Model §2](../science/model.md#2-the-prescribed-order-of-operations) and the [Gottesman bridge](GOTTESMAN_BRIDGE.md).

### M08 — Scalar effects and reversible branches

**Central question.** When can knowing a measurement outcome allow perfect restoration?

**Understanding target.** Explain why $M_s=p_sI$ means that the outcome reveals no logical-input information. For this one-Kraus-per-outcome instrument, connect scalar effects to branches proportional to isometries and to exact recovery. Recognize why the same scalar test would be insufficient for a general branch with unobserved Kraus multiplicity.

**Dependencies and story role.** M05, M07. This gives the operational recovery condition that the new code criterion will classify.

**Conceptual connection.** Preserving all logical superpositions is equivalent here to preventing the classical readout from resolving their amplitudes. The condition is an inherited recovery fact, with its source attribution retained.

**Source pointers.** G §2.2 for the no-logical-information principle behind correction; [Model §5](../science/model.md#5-exact-reversibility-and-logical-information-in-the-syndrome) for the precise branch theorem and its attribution.

### M09 — Polar recovery and the fidelity optimum

**Central question.** Which part of a known branch map can recovery undo?

**Understanding target.** Separate a branch into its polar unitary and positive factor. Explain why undoing the unitary leaves $\sqrt{M_s}$, and how its singular values determine the optimal entanglement fidelity. Include singular and zero branches without assuming an inverse of $K_s$ exists.

**Dependencies and story role.** M07, M08. Exact reversibility becomes a quantitative benchmark for imperfect recovery.

**Conceptual connection.** Rotations of a branch can be corrected with full knowledge of the unitary. Unequal singular values encode the residual distortion that determines the recovery optimum.

**Source pointers.** G §2.2 supplies correction context. The polar optimization is an additional project derivation of the attributed Gregoratti–Werner result, supplied in [Model §4](../science/model.md#4-exact-optimal-entanglement-fidelity) and the [bridge](GOTTESMAN_BRIDGE.md); it is not assigned as a theorem in Gottesman.

### M10 — Errors that include a reference system

**Central question.** How should recovery error be measured when the logical system may be entangled with something else?

**Understanding target.** Distinguish entanglement infidelity from half-diamond distance. Explain the optimization over recovery for each metric and the inequality $e_{\mathrm{opt}}\leq\varepsilon_{\mathrm{opt}}$. Understand that an entanglement-fidelity optimizer need not optimize diamond distance for a general unitary.

**Dependencies and story role.** M07, M09. This specifies the strength of the eventual upper and lower bounds.

**Conceptual connection.** A lower bound derived from a maximally entangled input already constrains every recovery in the stronger channel metric. An explicit recovered channel will later attain that bound.

**Source pointers.** G §2.2 gives the correction setting. The reference-system metrics and their comparison are project-supplied background in [Model §3](../science/model.md#3-performance-measures) and the [bridge](GOTTESMAN_BRIDGE.md), beyond the anchor's treatment.

### M11 — Clifford changes and physical locality

**Central question.** Which changes of Pauli basis preserve the physical support structure?

**Understanding target.** Follow Pauli operators under Clifford conjugation. Distinguish local Clifford axis changes and qubit permutations from an entangling CNOT. Explain why an entangling transformation of a code can change its relation to product-unitary noise even though the transformed subspace is unitarily equivalent.

**Dependencies and story role.** M02, M05, M06. This prepares both the local-Clifford witness and the later comparison of two code families.

**Conceptual connection.** The noise class refers to the physical division into qubits. Support-sensitive properties therefore contain information that abstract equivalence of code spaces alone does not preserve.

**Source pointers.** G §3.5; [Support criterion §§2 and 5](../science/support_criterion.md) and [Interpretation: the CNOT comparison](../science/interpretation.md#why-the-cnot-comparison-is-consequential).

### M12 — The complete operational question

**Central question.** What is optimized, what is fixed, and what must the code preserve?

**Understanding target.** State the full experiment: encode, apply a known product unitary, measure the original complete syndrome ideally, then recover. Explain the quantifier “for every specified unitary, there exists a recovery,” the retained classical record and quantum output, and the unrestricted resources allowed to recovery.

**Dependencies and story role.** M03–M11. This completes the first phase by assembling its local ingredients into the precise global task.

**Conceptual connection.** A prescribed readout can destroy information that a known unitary alone would preserve. The project measures the cost of this fixed order rather than assuming that the inverse rotation can be applied before measurement.

**Source pointers.** G §§2.1–2.3 for the QEC framework; [Model §§2 and 6](../science/model.md) and [Interpretation: operational boundaries](../science/interpretation.md#operational-boundaries) for the project's exact assumptions.

## Phase II — Narrative and innovation

### M13 — From syndrome effects to the full stabilizer group

**Central question.** Why can all syndrome effects be analyzed through conjugated stabilizers?

**Understanding target.** Explain the binary Fourier relation between the syndrome projectors and the full stabilizer group. Use its inverse to see that scalarity of every $M_s$ is equivalent to scalar compression of every conjugated stabilizer. Identify why checking only independent generators would leave the argument incomplete.

**Dependencies and story role.** M04, M12. This is the first reduction of the global recovery question to algebraic components.

**Conceptual connection.** The project uses the completeness of the measured record to translate an operational condition into a condition on a basis of commuting observables.

**Source pointers.** G §2.3 provides the stabilizer measurement basis; the Fourier reduction is proved in [Support criterion §1](../science/support_criterion.md#1-the-support-criterion).

### M14 — The orbit of one physical support

**Central question.** Which Pauli directions can a product unitary generate from one stabilizer?

**Understanding target.** Explain why conjugation rotates a nonidentity local axis within the span of $X,Y,Z$ and leaves an identity factor unchanged. Show that the linear span of the product-unitary orbit is exactly the Pauli space with the original support, with local Cliffords reaching its Pauli basis elements.

**Dependencies and story role.** M02, M06, M11. This identifies the local geometric restriction that makes the classification possible.

**Conceptual connection.** Product noise explores all axis choices on an existing support. The theorem can therefore use exact support equality, rather than only operator weight or support containment.

**Source pointers.** G §§2.1 and 3.5 supply Pauli and Clifford foundations; [Support criterion §1](../science/support_criterion.md#1-the-support-criterion) proves the orbit statement.

### M15 — Compression into the logical space

**Central question.** What can a Pauli operator become when restricted to the codespace?

**Understanding target.** Derive the three possibilities for $V^\dagger QV$: zero for a detectable Pauli, a signed identity for a stabilizer label, and a traceless nontrivial logical Pauli otherwise. Explain the role of anticommutation in the zero case and in proving that a nontrivial logical action cannot be scalar.

**Dependencies and story role.** M03, M05, M14. This converts the support orbit into a statement about logical information.

**Conceptual connection.** The first two compression outcomes cannot reveal a logical observable. The third is the precise obstruction to scalar syndrome effects.

**Source pointers.** G §§2.3 and 3.4 for detectable and logical actions; [Support criterion §1](../science/support_criterion.md#1-the-support-criterion), Eq. (4), for the trichotomy used here.

### M16 — The exact support classification

**Central question.** Which stabilizer codes remain reversible after every product rotation and the prescribed syndrome readout?

**Understanding target.** Combine the previous reductions to prove both directions: exact recovery for every product unitary is equivalent to disjoint support sets of stabilizers and nontrivial logical Paulis. Keep the universal unitary quantifier, exact support equality, and full stabilizer group visible in the statement.

**Dependencies and story role.** M08, M13, M14, M15. The local reductions now form the first complete global theorem.

**Conceptual connection.** A recovery condition expressed in channels becomes a finite property of the code's Pauli support structure, valid at all rotation angles.

**Source pointers.** G §§2.3, 3.4 and 3.5 are prerequisite background; [Support criterion §1, Theorem 1](../science/support_criterion.md#1-the-support-criterion) is the authoritative new classification.

### M17 — A support collision becomes a logical measurement

**Central question.** What physical loss witnesses failure of the support condition?

**Understanding target.** Map a colliding stabilizer to a nontrivial logical Pauli using local Cliffords. Explain why each reached syndrome branch then has logical rank at most $D/2$, so optimal entanglement fidelity is at most one half and optimal half-diamond error is at least one half.

**Dependencies and story role.** M09, M10, M11, M16. This turns necessity of the criterion into an operationally substantial finite-angle witness.

**Conceptual connection.** The syndrome learns a logical eigenvalue. Keeping that classical answer cannot reconstruct the coherence lost between the eigenvalue sectors.

**Source pointers.** G §§2.3 and 3.5 for measurement and conjugation; [Support criterion §2](../science/support_criterion.md#2-a-collision-exposes-a-logical-observable), with the logical-projection attribution in [Related work](../references/RELATED_WORK.md).

### M18 — The reversibility distance and local anticommutation

**Central question.** How many local axis changes are needed before a support collision can expose a logical operator?

**Understanding target.** Define $r$ as the minimum weight of $gL$ over equal-support stabilizer/logical pairs, with an empty minimum equal to infinity. Translate that definition into a logical Pauli $E$ and a stabilizer that anticommutes with every active local factor of $E$. Explain why finite $r$ is even and at least $d$.

**Dependencies and story role.** M05, M16. The exact yes-or-no classification acquires a quantity that can control weak-noise sensitivity.

**Conceptual connection.** Global commutation forces an even number of local anticommutations. The physical arrangement of those local factors adds a restriction absent from ordinary distance.

**Source pointers.** G §§2.3, 3.1 and 3.4 for commutation and logical classes; [Support criterion §3](../science/support_criterion.md#3-syndrome-reversibility-distance) defines and analyzes $r$.

### M19 — Syndrome columns and logical circuits

**Central question.** When does a fixed Pauli word contain a smaller undetectable subword?

**Understanding target.** Associate a binary syndrome column to each active local Pauli factor. Identify subwords with subsets of columns and zero syndrome with a binary dependence. Define a logical circuit by the absence of proper nonempty normalizer subwords, retaining the original local Pauli labels.

**Dependencies and story role.** M04, M05, M18. This gives the combinatorial language that identifies irreducible witnesses.

**Conceptual connection.** A circuit is a minimal dependence in a binary matrix, rather than a sequence of quantum gates. Its minimality controls exactly which coherent error paths can share a syndrome.

**Source pointers.** G §3.1, Eq. (56) and the syndrome linear-algebra discussion; [Support criterion §4](../science/support_criterion.md#4-equivalent-characterization-by-even-logical-circuits).

### M20 — Why even circuits determine the invariant

**Central question.** Why are the support-collision minimum and the even-circuit minimum identical?

**Understanding target.** Explain both directions: a minimum local-anticommutation witness cannot contain a proper normalizer subword, and an even logical circuit supplies the required stabilizer through binary row-space orthogonality. Deduce $r=d$ for even-distance codes while preserving the distinction between a circuit and an arbitrary even-weight logical Pauli.

**Dependencies and story role.** M18, M19. This completes the combinatorial characterization needed for a sharp lower bound.

**Conceptual connection.** Minimum distance already gives a logical circuit. Its parity decides whether it directly provides the stabilizer witness for measurement-induced irreversibility.

**Source pointers.** G §§3.1 and 3.4 supply binary and logical foundations; [Support criterion §§4–5](../science/support_criterion.md) proves the equivalence and distance corollary.

### M21 — Counting the powers of a small rotation

**Central question.** Why can no nonscalar syndrome contribution appear below order $r$?

**Understanding target.** Bound each changed local Pauli-axis coefficient by the rotation scale. Apply the support and compression structure to bound the nonscalar part $M_s-p_sI$. Interpret the collision polynomial $B_S(x)$ as an explicit finite-code upper certificate and identify its lowest nonzero power as $r$.

**Dependencies and story role.** M06, M13, M15, M18. This connects the combinatorial invariant to perturbation order without choosing a decoder.

**Conceptual connection.** A logical contribution needs a permitted pattern of local axis changes. The polynomial records those patterns and makes the code dependence of the prefactor visible.

**Source pointers.** G §2.1 for local error expansions; the quantitative coefficient and polynomial bounds are in [Error exponent §2](../science/error_exponent.md#2-a-finite-code-collision-polynomial).

### M22 — A recovery bound that includes rare outcomes

**Central question.** How can square-root branch maps retain the exponent when syndrome probabilities are small?

**Understanding target.** Use the positive Kraus operators of polar recovery and the double-commutator identity for the corrected channel. Explain why subtracting scalar square roots and then squaring the operator-norm deviation yields a bound in terms of $\sum_s\Vert M_s-p_sI\Vert $, without dividing by $p_s$.

**Dependencies and story role.** M09, M10, M21. This supplies the uniform-in-axis recovery upper bound and closes the rare-syndrome issue.

**Conceptual connection.** The cancellation imposed by trace preservation matters at the channel level. Estimating each square root in isolation would obscure the full error exponent.

**Source pointers.** G §2.2 is the correction context. The additional operator-norm argument is developed in [Error exponent §2.3](../science/error_exponent.md#23-polar-recovery-without-loss-of-exponent) and the [bridge](GOTTESMAN_BRIDGE.md).

### M23 — Complementary coherent error paths

**Central question.** Why does a minimum circuit turn a general code into a repetition-type instrument?

**Understanding target.** Rotate about the local factors of a minimum even logical circuit. Show that two error subsets have the same syndrome exactly when they coincide or complement one another. Derive the resulting two-term logical branch and explain why even circuit weight makes its relative amplitude real.

**Dependencies and story role.** M06, M19, M20. This constructs an explicit physical witness for the recovery lower bound.

**Conceptual connection.** Circuit minimality fixes the entire interference pattern, so the witness can be analyzed exactly without simulating every code state.

**Source pointers.** G §2.1 for coherent Pauli expansion and §3.1 for syndrome addition; [Error exponent §§3–3.1](../science/error_exponent.md#3-an-embedded-even-repetition-witness).

### M24 — The exact repetition tail and recovery converse

**Central question.** What loss remains even after the best possible recovery of the circuit witness?

**Understanding target.** Compute the optimal fidelity contribution of a two-term branch, count complementary pairs including equal-weight ties, and obtain $p_r(q)$. Explain how majority Pauli recovery realizes a dephasing channel and why the fidelity lower bound proves that arbitrary CPTP recovery cannot achieve a smaller half-diamond error.

**Dependencies and story role.** M09, M10, M23. The witness provides an exact optimized loss, not merely a failure of one selected decoder.

**Conceptual connection.** The established even-repetition mechanism becomes a universal embedded witness for every code with finite $r$. Equality of the two optimized error metrics is proved for this witness; it is not assumed for all instruments.

**Source pointers.** G §§2.1–2.2 provide the correction background; [Error exponent §3.2](../science/error_exponent.md#32-exact-optimal-error) gives the calculation and all-recovery converse. [Related work](../references/RELATED_WORK.md#huangdohertyflammia-coherent-error-scaling) credits the repetition mechanism.

### M25 — The sharp error law as one theorem

**Central question.** How do the upper certificate and exact witness establish the worst-case exponent?

**Understanding target.** Combine the uniform upper bound with the positive leading coefficient of the repetition tail to derive $\Theta(|\theta|^r)$ in both optimized metrics. State the supremum over bounded local Hamiltonians and the separate optimization of recovery for each unitary. Distinguish the fixed-code limit from uniform large-code estimates.

**Dependencies and story role.** M20, M21, M22, M23, M24. This joins the quantitative mechanisms into the second central theorem.

**Conceptual connection.** The support invariant simultaneously classifies exact reversibility through $r=\infty$ and measures the first unavoidable order of loss when $r$ is finite.

**Source pointers.** G §2.1 supplies the local-rotation context; [Error exponent §§1 and 4–5](../science/error_exponent.md) states the sharp law, its limits and the finite-angle consequence.

### M26 — The five-qubit seed and useful logical representatives

**Central question.** Which concrete properties of the five-qubit code make the family construction possible?

**Understanding target.** Use its fifteen distinct nonzero single-site syndromes to exclude normalizers of weight one or two. Identify minimum representatives of all three logical classes. Prove that stabilizers have even weight and every nontrivial logical representative has odd weight, and retain the chosen logical supports needed at the CNOT endpoints.

**Dependencies and story role.** M02, M05, M11, M16. This introduces the established code as the concrete ingredient of the project's separation.

**Conceptual connection.** Blockwise parity excludes matching supports, while carefully chosen minimum representatives make a later entangling deformation fully explicit.

**Source pointers.** G §3.1, Table 4, for the five-qubit stabilizer; [Code families §2](../science/code_families.md#2-the-five-qubit-building-block) for the required parity, distance and endpoint arguments. The underlying code is attributed in [Related work](../references/RELATED_WORK.md#five-qubit-code-and-stabilizer-foundations).

### M27 — Concatenation and blockwise reversibility

**Central question.** How do the seed properties survive while the ordinary distance grows?

**Understanding target.** Lift outer logical operators through inner blocks to prove the parameters $[[5^t,1,3^t]]$ and preserve odd logical/even stabilizer parity. Explain why two such blocks have $r=\infty$ by comparing supports within a block on which a logical operator is nontrivial. Recognize why total parity alone is insufficient for the two-block code.

**Dependencies and story role.** M20, M26. This turns the finite seed into the exactly reversible side of an unbounded-distance comparison.

**Conceptual connection.** Recursive encoding transfers local properties to a family, while tensor-product structure preserves a stronger blockwise obstruction to support collisions.

**Source pointers.** G §5.2, its concatenation discussion, and §3.4 for logical classes; [Code families §§3–4](../science/code_families.md) supplies the family's specific proofs.

### M28 — Preserving distance with one physical CNOT

**Central question.** Why does this particular entangling deformation leave ordinary distance unchanged?

**Understanding target.** Follow a hypothetical shorter logical operator backward through the CNOT. Use its at-most-one weight change and the product-code distance to obtain the contradiction that a minimum logical must both touch and avoid the other block. Exhibit the untouched minimum logical for the matching upper bound.

**Dependencies and story role.** M11, M27. This holds the conventional parameters fixed before comparing the new invariant.

**Conceptual connection.** General Clifford equivalence alone does not preserve physical weight. The construction needs its explicit endpoint choices and a dedicated distance argument.

**Source pointers.** G §3.5 for CNOT conjugation; [Code families §5](../science/code_families.md#5-one-cnot-preserves-the-distance) for the exact lower and upper distance bounds.

### M29 — The new even circuit and the family separation

**Central question.** How does the same CNOT create the smallest finite reversibility distance allowed by odd $d$?

**Understanding target.** Show that a minimum odd target-block logical becomes an even logical of weight $3^t+1$. Analyze each possible subword under inverse CNOT to prove circuit minimality. Combine that witness with parity and the preserved distance to obtain $r=3^t+1$, and apply the exact repetition formula at every level.

**Dependencies and story role.** M20, M25, M27, M28. This completes the same-parameter separation and its operational error consequence.

**Conceptual connection.** A single physical entangling change can create a new way for the prescribed readout to learn logical information, while $n$, $k$ and $d$ remain identical.

**Source pointers.** G §§3.4–3.5 provide logical and conjugation background; [Code families §§6–7](../science/code_families.md) and [Error exponent §3](../science/error_exponent.md#3-an-embedded-even-repetition-witness) establish the circuit and its exact error law.

### M30 — The complete scientific conclusion

**Central question.** What does syndrome reversibility distance add to the understanding of a quantum code?

**Understanding target.** Connect the information leaked by syndrome readout, the exact support criterion, the even-circuit invariant, the sharp optimal error exponent, and the growing-distance CNOT pair. Separate inherited recovery and repetition results from the additional code-level conclusions. State the assumptions about known noise, fixed measurement order, ideal readout, unrestricted recovery and code-dependent constants.

**Dependencies and story role.** M12, M17, M25, M29. This completes the second phase by integrating the local mechanisms into the scope and significance of the whole project.

**Conceptual connection.** The result identifies an information-preservation property of a prescribed measurement relative to physical locality. Finite calculations support the specified examples; the universal classification and asymptotic conclusions rest on their proofs. The shared conventional code parameters do not determine this property.

**Source pointers.** G §§2.2–2.3 for the foundational information-preservation question, with §5.2 as context for the distinct fault-tolerance setting; [Interpretation](../science/interpretation.md), [Related work](../references/RELATED_WORK.md), and [Verification](../verification/README.md) provide the final synthesis, attribution and evidence boundaries.

