# Learning the syndrome reversibility problem

A known local rotation is reversible before measurement. Measuring a stabilizer code's original syndrome can nevertheless expose a logical observable and destroy information that later recovery cannot restore. The learning route explains the exact support criterion for avoiding that loss, the sharp weak-rotation exponent when it occurs, and the same-parameter code families that separate it from ordinary distance.

This project and [Error-Transparent Locality](https://github.com/GoGoKo699/Error-Transparent-Locality) share one educational anchor: **John Preskill, *Quantum Information*, Chapter 7, “Quantum Error Correction,” updated March 2026**. The [author's course index](https://preskill.caltech.edu/ph229/index.html) links the [official chapter](https://www.preskill.caltech.edu/ph219/chap7_26_6A.pdf).

The common foundation is §§7.2–7.4 for exact correction, information leakage, reference systems and distance; §§7.10–7.11 for stabilizers, encoded Paulis and the five-qubit code; and §7.17 for concatenation. Sections 7.6–7.8 supply classical-code and CSS context. The readout-first bridge translates Preskill's real-Pauli convention into the project's Hermitian tensors, derives the needed local Clifford and CNOT rules, and develops the recovery optimization, channel metrics and project-specific proofs locally. The [companion gate bridge](https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/learning/PRESKILL_BRIDGE.md) follows the same foundations into protected dynamics.

| Route | What it supplies |
|---|---|
| [Thirty milestones](MILESTONES.md) | Exactly 12 background/preliminaries and 18 narrative/innovation targets, progressing from individual objects to the complete scientific account |
| [Preskill bridge](PRESKILL_BRIDGE.md) | Verified anchor sections, notation translation, a complete two-qubit example, recovery and metric derivations, and the support/circuit mechanism |
| [Measured channel and recovery](../science/model.md) | Exact ordering, informed-recovery quantifiers, the recorded instrument and optimal entanglement fidelity |
| [Support criterion](../science/support_criterion.md) | The if-and-only-if classification and minimum even logical-circuit characterization |
| [Sharp error exponent](../science/error_exponent.md) | The rare-branch upper bound, all-recovery repetition witness and fixed-code asymptotic law |
| [Same-parameter code families](../science/code_families.md) | Concatenated product codes with infinite reversibility distance and their one-physical-CNOT conjugates with finite distance |
| [Gate/readout connection](../science/gate_readout_tradeoff.md) | The consequence of minimum protected-gate interaction order for the readout invariant |
| [Protected entanglers](../science/protected_entangling_gates.md) | Exact protected encoded control on an unchanged product code, with its original readout retained |
| [Related work](../references/RELATED_WORK.md) and [verification](../verification/README.md) | Attribution, exact finite certificates and the boundaries of computational checks |

Milestones M01–M12 assemble the complete readout experiment. M13–M24 establish and test its classification, sharp error law and code-family separation. M25–M30 connect the readout invariant to protected gates, optimal one-, two- and three-error families, and fixed-code logical entanglers. These final connections compare two tasks on the same code; they retain the different assumptions of gate faults and the known-unitary readout experiment.

For the complete scientific statements and proofs, begin with the [measured-channel model](../science/model.md). For physical meaning and operational boundaries, use [Syndrome measurement and logical information](../science/interpretation.md).
