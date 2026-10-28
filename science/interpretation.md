# Syndrome measurement and logical information

The [model](model.md) fixes a physical tensor-product structure, a stabilizer code, and an order of operations. The [support theorem](support_criterion.md) classifies exact recoverability in that model. The [exponent theorem](error_exponent.md) quantifies the first unavoidable loss, and the [code families](code_families.md) show why ordinary code parameters cannot determine it.

## A syndrome is also a measurement on the encoded input

After a known product unitary $U$, the syndrome probability for an encoded logical state $\rho$ is

```math
\Pr(s\mid\rho)=\mathrm{Tr}(M_s\rho),
\qquad M_s=V^\dagger U^\dagger\Pi_sUV.
```

Thus the prescribed physical readout induces a POVM on the logical input. If $M_s=p_s I$ for every outcome, the syndrome contains no information about that input. Each conditional quantum transformation is then a scaled isometry, so a recovery that knows $U$ and $s$ can undo it. If an effect is nonscalar, the outcome probabilities depend on the logical input. Its conditional transformation also changes relative singular values, and deterministic recovery cannot restore arbitrary logical states exactly.

This equivalence between outcome information and reversible branches is inherited from outcome-assisted quantum recovery; see Gregoratti–Werner and the QEC specialization in Liu–Chen in [Related work](../references/RELATED_WORK.md). The additional step in this project is to determine exactly when scalarity holds for **every** physical product unitary, directly from the code's stabilizer and logical supports.

The distinction between classical information and coherent quantum information matters. The syndrome is retained as a classical record, but the measurement has removed coherence between different outcomes. Retaining the outcome allows conditional correction; it does not restore those off-diagonal terms. If the syndrome has learned a logical Pauli eigenvalue, states superposed across its eigenspaces have lost phase information.

## Why support enters

A single-qubit rotation changes the axis of a nonidentity Pauli while leaving the identity unchanged. A product rotation therefore moves a stabilizer inside the linear space of Pauli operators with exactly the same support. Compressing that space into the code yields three possibilities: zero for detectable Paulis, scalar operators for stabilizers, and nontrivial encoded observables for logical Paulis.

The support condition says that the third possibility never occurs. Conversely, a matching-support logical Pauli can be reached from the stabilizer by local Clifford changes of axis. The syndrome then reads a logical observable. This establishes an exact boundary, including arbitrarily large product rotations; it is not merely a perturbative sufficient condition.

## What the exponent measures

The ordinary distance $d$ counts the smallest number of physical sites supporting any nontrivial logical Pauli. The syndrome reversibility distance $r$ asks which of those logical operators can be exposed through this particular measurement mechanism. A qualifying operator must have a stabilizer that anticommutes locally with every one of its active factors. The total number of such local anticommutations is even because the logical operator commutes with the stabilizer globally.

A minimum qualifying operator is a circuit of the binary syndrome map: no proper nonempty subword has zero syndrome. Local rotations about its factors consequently produce complementary pairs of physical error paths with the same syndrome. Their real relative amplitudes make the measurement weakly resolve one encoded Pauli. The established even-repetition mechanism then gives an exact error probability, while the general circuit theorem proves that every code with finite $r$ contains such a witness.

This mechanism explains both sides of the exponent theorem. Changing enough local axes to expose logical information costs at least $r$ powers of the small rotation scale. The circuit witness realizes loss at exactly that order. Polar recovery controls the upper bound without dividing by small syndrome probabilities, so the conclusion includes rare branches.

Finite $r$ is even and at least $d$. For even $d$, every minimum logical operator is already an even circuit, giving $r=d$. Odd distance allows further suppression or exact reversibility. The exponent is defined after optimization over quantum recovery; a fixed Pauli decoder can leave a coherent logical rotation that the allowed, fully informed recovery would remove.

## Why the CNOT comparison is consequential

The two families encode two qubits with the same growing ordinary distance. Their difference is one physical CNOT in the definition of the code. The product family has no matching stabilizer/logical supports. Conjugating it by that CNOT turns a minimum odd logical operator into an even circuit of weight $d+1$, creating the smallest finite exponent consistent with odd distance.

An entangling Clifford transformation preserves the stabilizer formalism but can change how operator support meets the physical division into qubits. Product-unitary noise is defined using that division. Consequently, the invariant is preserved by physical permutations and local Clifford changes of basis, while an entangling Clifford can change it. There is no contradiction with unitary equivalence of the two code spaces: conjugating the **entire** experiment also changes the physical noise class.

The construction is an unbounded-distance information-preservation separation. It keeps $k=2$ as the block size grows. Its proof does not require a large simulation, an assumed decoder, or an asymptotic threshold estimate.

## A constraint on protected-gate interaction order

The [joint gate/readout theorem](gate_readout_tradeoff.md) relates this invariant to exact error-transparent logical Hamiltonians. For odd $d$, attaining the interaction-order lower bound $\ell=d+t$ forces $r=d+1$, even for degenerate stabilizer codes with fixed ancillary qubits. Therefore $r>d+1$, including infinity, requires $\ell\ge d+t+1$. The nine-qubit four-body gate supplies an explicit $r=4$ example and a logical-$Z$ readout-loss witness. The five-qubit code attains the stronger five-body bound while keeping $r=\infty$. Higher-error attainability and a joint noise experiment require additional results.

The [concatenated equality construction](minimum_readout_safe_gates.md) attains this stronger bound for one protected error at every distance $d=3^m$. Its codes retain $r=\infty$, while exact logical $X$, $Y$, and $Z$ generators use interaction order $d+2$ and coupling budget $d+1$ per unit logical half-splitting. The readout guarantee still permits recovery to know the product unitary. The [sparse construction](sparse_protected_gates.md) retains these equalities with at most $4n$ native Pauli terms and polynomially bounded rational denominators; the 25-qubit example uses 75 terms.

The [two-error gate construction](two_error_gates.md) strengthens this coexistence on the twenty-five-qubit code: all physical errors of weight at most two are transparent at interaction order thirteen, while $r=\infty$ remains unchanged. An inner-block moment proves the exact coupling optimum $521/8$ per unit logical half-splitting. The stored sparse form has 2,575 native terms per axis. This gate guarantee still shares the code with the prescribed-readout statement under their distinct operational assumptions.

## Operational boundaries

**Measurement order.** The question fixes recovery after syndrome readout. If the known inverse $U^\dagger$ could be applied beforehand, it would restore the code state directly. The present model isolates the information cost of committing to the original syndrome measurement before that reversal. One possible operational setting is a fixed extraction cycle with calibration information available only at recovery time; that motivation is an interpretation of the model, not a demonstrated hardware protocol.

**Knowledge available to recovery.** The quantifier is “for each $U$, there exists a recovery.” It supplies an information-theoretic benchmark when the unitary is known. A single recovery that works without knowing $U$ is a different requirement. The zero-error statement does not estimate the resources needed to compute or implement the conditional correction.

**Measurement and noise assumptions.** The measurement consists of all independent original stabilizer checks and is noiseless. Products of arbitrary one-qubit unitaries are allowed. Correlated multiqubit noise, imperfect extraction, and access to a coherent unmeasured syndrome register change the channel being analyzed.

**Asymptotics.** The exponent theorem takes $\theta\to0$ for a fixed finite code. Its explicit collision-polynomial bound identifies the code dependence of the prefactor. The growing family establishes exact values of $r$ and an exact witness channel at every size; it does not interchange the small-angle and large-code limits. Uniform thresholds require a separate control of those limits.

## Attribution and verification

The branch recovery framework, its optimal entanglement fidelity, and the even-repetition mechanism retain their original attribution. The present results combine the physical-support classification, circuit equivalence, optimized exponent, and same-parameter construction into one code-level account. The [source comparison](../references/RELATED_WORK.md) states the operational differences from prior symmetry criteria, coherent-error analyses, enumerators, and threshold theorems.

The [verification suite](../verification/README.md) checks concrete Pauli identities, distances and circuit witnesses, and the full physical measurement calculation for the ten-qubit instance. These computations support their finite assertions. The general conclusions rest on the proofs in the linked scientific pages.

