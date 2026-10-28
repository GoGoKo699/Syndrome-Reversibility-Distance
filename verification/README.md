# Reproducing the finite calculations

The suite checks the finite Pauli algebra, explicit members of the concatenated family, and a physical realization of the ten-qubit measurement instrument. The universal support criterion and error exponent are proved in the [scientific account](../science/support_criterion.md); finite calculations are used to check the stated examples and their implementation.

## Run

The recorded environment is Python **3.12.14** with NumPy **2.3.5**. NumPy is the only external dependency and is pinned in [requirements.txt](requirements.txt).

From the repository root:

```bash
python3 -m pip install -r verification/requirements.txt
python3 verification/run_all.py
```

The runner executes all four programs, checks their mathematical assertions, and compares their full results with the committed records in [expected/](expected/). The algebraic outputs are compared exactly. Numerical outputs use absolute tolerance $2\times10^{-12}$; runtime-version metadata is excluded from comparison. Success prints `"expected_records_matched": true` along with the results. Run Python normally so that assertions remain enabled.

To save a fresh set of output records separately:

```bash
python3 verification/run_all.py --write-results /tmp/syndrome-reversibility-results
```

The default command does not alter committed outputs. A failing assertion or an output mismatch produces a nonzero exit status. Each program can also be executed individually with `python3 verification/<program>.py`.

## What each calculation checks

| Program | Construction and assertions | Saved output |
|---|---|---|
| [verify_code_parameters.py](verify_code_parameters.py) | Enumerates the complete stabilizer and centralizer of the five-qubit seed and both ten-qubit codes; calculates distance and the support-collision definition of $r$ directly; checks the explicit witness and its proper subwords | [Exact parameters](expected/verify_code_parameters.json) |
| [verify_concatenated_family.py](verify_concatenated_family.py) | Builds recursive generators for levels $t=1,2,3$; checks rank, commutation, logical cosets, minimum-weight representatives, the unchanged logical, and the transformed circuit rank | [Family structure](expected/verify_concatenated_family.json) |
| [verify_physical_channel.py](verify_physical_channel.py) | Builds physical code states, applies single-qubit rotations, projects all 256 signed syndrome sectors, computes branch singular values, and applies physical Pauli recovery | [Physical channel](expected/verify_physical_channel.json) |
| [verify_gate_readout_tradeoff.py](verify_gate_readout_tradeoff.py) | Independently enumerates the nine-qubit protected-gate code, five-qubit code, and fixed-ancilla extension; verifies support collisions, logical-$Z$ circuit and all proper subwords, exact Hamiltonian coefficient identities, and the $p_4$ polynomial | [Gate/readout comparison](expected/verify_gate_readout_tradeoff.json) |

The helper [pauli_tools.py](pauli_tools.py) uses exact binary arithmetic and discards Pauli phases for support and commutation calculations. It represents physical qubit 0 by the leftmost printed character and the least significant computational-basis bit. The physical-channel program independently retains signs by applying the original generator matrices and conjugating their action with the actual CNOT permutation. It never constructs a physical projector from a signless Pauli label after entangling conjugation.

## Recorded results

Complete enumeration finds the following values directly from their definitions:

| Code | Physical qubits | Encoded qubits | Distance | $r$ | Support-collision pairs |
|---|---:|---:|---:|---:|---:|
| Five-qubit seed | 5 | 1 | 3 | $\infty$ | 0 |
| Product of two seeds | 10 | 2 | 3 | $\infty$ | 0 |
| Product conjugated by CNOT $0\to5$ | 10 | 2 | 3 | 4 | 4,608 |
| Nine-qubit four-body protected-gate code | 9 | 1 | 3 | 4 | 2,848 |
| Five-qubit seed with a fixed $\lvert0\rangle$ ancilla | 6 | 1 | 3 | $\infty$ | 0 |

For the transformed family at $t=1,2,3$, the constructed even logical circuits have weights $4,10,28$ and syndrome-column ranks $3,9,27$, respectively. Their full subword kernel is therefore one dimensional. The family script does not exhaustively search distances at 50 or 250 physical qubits; the lower bounds for all levels are established in [the family proof](../science/code_families.md).

The physical-channel calculation uses a $1024\times4$ code isometry. For rotations at $\theta=0.1,0.2,0.4,0.8,\pi/2$, it finds exactly eight active syndrome sectors and verifies

```math
e_{\mathrm{opt}}(U_\theta)
=p_4(q)=3q^2-2q^3,
\qquad q=\sin^2(\theta/2).
```

In the committed run, the largest absolute infidelity discrepancy is below $5\times10^{-16}$. The normalized Choi matrix of the physically recovered channel agrees with

```math
\Lambda_\theta(\rho)
=(1-p_4(q))\rho+p_4(q)\overline E\rho\overline E
```

to a Frobenius residual below $1.3\times10^{-15}$. A generic set of nonidentical rotations on the original product code yields a maximum nonscalar-effect residual below $4.8\times10^{-16}$. Slight numerical excursions above unit fidelity in these records are floating-point roundoff and are retained without clipping.

The [joint gate/readout theorem](../science/gate_readout_tradeoff.md) uses the last two entries to check its construction and its fixed-qubit caveat. The standalone script is identical in both companion repositories and requires only the standard library. For the nine-qubit witness $Z_{1368}$ (one-based indexing), it verifies all fourteen proper subwords have nonzero syndrome and the sixteen rotation subwords form eight complementary syndrome pairs. Its exact $p_4$ polynomial calculation checks the circuit formula algebraically; it is separate from the ten-qubit physical matrix calculation above.

## What is analytical and what is numerical

The numerical program constructs the physical instrument without assuming the repetition formula. It obtains optimal entanglement fidelity from the branch singular values using the [proved and attributed recovery formula](../science/model.md#4-exact-optimal-entanglement-fidelity). It then checks the entire recovered channel against the predicted dephasing channel.

The program does not run a numerical optimization over all quantum recoveries. The all-recovery fidelity optimum, its lower bound on half-diamond error, and the exact diamond distance of a Pauli dephasing channel are analytical results in [the exponent proof](../science/error_exponent.md). A Choi-state trace distance alone would not establish a generic diamond distance; here the verified dephasing form is what allows the analytic equality to apply.

Likewise, enumeration verifies the listed finite codes. The exact-support classification, circuit equivalence, uniform-in-axis fixed-code upper bound, and unbounded-distance construction depend on their written proofs. The suite provides no separate claim about efficient decoding, noisy extraction, or uniform threshold constants.

