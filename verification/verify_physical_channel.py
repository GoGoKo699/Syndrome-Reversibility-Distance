"""Direct 1024-dimensional check of the ten-qubit irreversible channel.

The channel is built from physical code states, local rotations, and all
256 signed syndrome projectors. Neither a repetition-code channel nor
the claimed logical error formula is used to construct those objects.
The projected states then determine optimal entanglement fidelity and
the channel produced by minimum-weight Pauli recovery.
"""

from __future__ import annotations

import json
import platform

import numpy as np

from pauli_tools import (
    FIVE_GENERATORS, cnot, from_string, local_factors,
    product_five_generators, span, syndrome, weight,
)


TOLERANCE = 2e-12


def apply_pauli(pauli: int, vectors: np.ndarray, n: int) -> np.ndarray:
    mask = (1 << n) - 1
    x, z = pauli & mask, pauli >> n
    basis = np.arange(1 << n)
    phases = (1j) ** (x & z).bit_count() * np.array(
        [(-1.0) ** (int(b) & z).bit_count() for b in basis]
    )
    result = np.zeros_like(vectors, dtype=complex)
    result[basis ^ x] = phases[:, None] * vectors
    return result


def apply_cnot(vectors: np.ndarray, control: int, target: int) -> np.ndarray:
    basis = np.arange(len(vectors))
    permutation = basis ^ (((basis >> control) & 1) << target)
    return vectors[permutation]


def five_qubit_isometry() -> np.ndarray:
    zero = np.zeros((32, 1), dtype=complex)
    zero[0, 0] = 1
    for word in FIVE_GENERATORS:
        zero = (zero + apply_pauli(from_string(word), zero, 5)) / 2
    zero /= np.linalg.norm(zero)
    one = apply_pauli(from_string("XXXXX"), zero, 5)
    result = np.concatenate((zero, one), axis=1)
    assert np.linalg.norm(result.conj().T @ result - np.eye(2)) < TOLERANCE
    return result


def local_rotations(vectors: np.ndarray, rotations: list[tuple[int, float, np.ndarray]]) -> np.ndarray:
    n = int(np.log2(len(vectors)))
    result = vectors.copy()
    for qubit, angle, axis in rotations:
        axis = np.asarray(axis, dtype=float)
        axis = axis / np.linalg.norm(axis)
        h_result = np.zeros_like(result)
        for coefficient, letter in zip(axis, "XYZ"):
            if coefficient:
                word = "I" * qubit + letter + "I" * (n - qubit - 1)
                h_result += coefficient * apply_pauli(from_string(word), result, n)
        result = np.cos(angle / 2) * result - 1j * np.sin(angle / 2) * h_result
    return result


def project_syndromes(vectors: np.ndarray, generators: list[int], transformed: bool) -> dict[int, np.ndarray]:
    """Build every projector physically; preserve the signs under CNOT."""
    states = {0: vectors}
    for j, generator in enumerate(generators):
        following = {}
        for label, state in states.items():
            if transformed:
                acted = apply_cnot(apply_pauli(generator, apply_cnot(state, 0, 5), 10), 0, 5)
            else:
                acted = apply_pauli(generator, state, 10)
            following[label] = (state + acted) / 2
            following[label | (1 << j)] = (state - acted) / 2
        states = following
    return states


def effects_and_fidelity(states: dict[int, np.ndarray]) -> tuple[list[np.ndarray], float]:
    effects = [state.conj().T @ state for state in states.values()]
    assert np.linalg.norm(sum(effects) - np.eye(4)) < TOLERANCE
    # Singular values of the physical branch map are square roots of its
    # effect eigenvalues, without unstable square roots of numerical zeros.
    fidelity = sum(np.linalg.svd(state, compute_uv=False).sum() ** 2
                   for state in states.values()) / 16
    return effects, float(fidelity)


def choi(kraus: list[np.ndarray]) -> np.ndarray:
    return sum(np.outer(k.reshape(-1, order="F"), k.reshape(-1, order="F").conj())
               for k in kraus) / 4


def main() -> dict:
    v5 = five_qubit_isometry()
    original_v = np.kron(v5, v5)
    transformed_v = apply_cnot(original_v, 0, 5)
    generators = product_five_generators()
    transformed_generators = [cnot(g, 0, 5, 10) for g in generators]
    assert np.linalg.norm(transformed_v.conj().T @ transformed_v - np.eye(4)) < TOLERANCE
    for generator in generators:
        acted = apply_cnot(apply_pauli(generator, apply_cnot(transformed_v, 0, 5), 10), 0, 5)
        assert np.linalg.norm(acted - transformed_v) < TOLERANCE

    witness = from_string("ZIIIIZYIIY")
    logical = transformed_v.conj().T @ apply_pauli(witness, transformed_v, 10)
    assert abs(np.trace(logical)) < TOLERANCE
    assert np.linalg.norm(logical @ logical - np.eye(4)) < TOLERANCE
    representatives: dict[int, int] = {}
    for error in sorted(span(local_factors(witness, 10)), key=lambda p: (weight(p, 10), p)):
        representatives.setdefault(syndrome(error, transformed_generators, 10), error)
    assert len(representatives) == 8

    records = []
    identity_choi = choi([np.eye(4)])
    for angle in (0.1, 0.2, 0.4, 0.8, np.pi / 2):
        rotations = [(0, angle, np.array([0, 0, 1])),
                     (5, angle, np.array([0, 0, 1])),
                     (6, angle, np.array([0, 1, 0])),
                     (9, angle, np.array([0, 1, 0]))]
        states = project_syndromes(local_rotations(transformed_v, rotations), generators, True)
        effects, fidelity = effects_and_fidelity(states)
        active = [label for label, state in states.items() if np.linalg.norm(state) ** 2 > 1e-20]
        assert set(active) == set(representatives)
        recovered_kraus = []
        for label in active:
            returned = apply_pauli(representatives[label], states[label], 10)
            logical_kraus = transformed_v.conj().T @ returned
            assert np.linalg.norm(returned - transformed_v @ logical_kraus) < TOLERANCE
            recovered_kraus.append(logical_kraus)
        assert np.linalg.norm(sum(k.conj().T @ k for k in recovered_kraus) - np.eye(4)) < TOLERANCE
        q = float(np.sin(angle / 2) ** 2)
        predicted = 3 * q ** 2 - 2 * q ** 3
        recovered_choi = choi(recovered_kraus)
        dephasing_choi = (1 - predicted) * identity_choi + predicted * choi([logical])
        channel_residual = float(np.linalg.norm(recovered_choi - dephasing_choi))
        choi_distance = float(np.abs(np.linalg.eigvalsh(recovered_choi - identity_choi)).sum() / 2)
        assert abs((1 - fidelity) - predicted) < TOLERANCE
        assert channel_residual < TOLERANCE
        assert abs(choi_distance - predicted) < TOLERANCE
        records.append({"theta": float(angle), "q": q,
                        "optimal_entanglement_infidelity": float(1 - fidelity),
                        "predicted_error": predicted, "active_syndromes": len(active),
                        "recovered_choi_distance_from_identity": choi_distance,
                        "dephasing_channel_residual_frobenius": channel_residual})

    generic_rotations = [(j, 0.173 * (j + 1), np.array([1 + j, 2 - j, 3 + (j % 3)]))
                         for j in range(10)]
    states = project_syndromes(local_rotations(original_v, generic_rotations), generators, False)
    effects, fidelity = effects_and_fidelity(states)
    nonscalar = float(max(np.linalg.norm(m - np.trace(m) * np.eye(4) / 4) for m in effects))
    assert nonscalar < TOLERANCE and abs(1 - fidelity) < TOLERANCE
    return {
        "check": "full_physical_channel",
        "physical_dimension": 1024, "logical_dimension": 4,
        "syndrome_sectors": len(states), "tolerance": TOLERANCE,
        "witness_rotations": records,
        "product_code_generic_rotation": {"maximum_nonscalar_effect_frobenius": nonscalar,
                                           "optimal_entanglement_fidelity": fidelity},
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "checks_passed": True,
    }


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True))
