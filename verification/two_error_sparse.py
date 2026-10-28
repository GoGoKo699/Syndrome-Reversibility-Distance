"""Exact sparse certificate for the twenty-five-qubit two-error gate.

A stored list of cyclic-orbit seeds is expanded and checked using integer
arithmetic. No optimizer or numerical package is needed. A numerical search
selected the support; all coefficients and physical intertwiners below are
independently certified over the rationals.
"""

from fractions import Fraction
from itertools import combinations, product
from math import lcm
from pathlib import Path
import json

import verify_readout_safe_family as family


HERE = Path(__file__).resolve().parent
N = 25
MASK = (1 << N) - 1
EXPECTED_MASSES = {
    "9": Fraction(143, 8),
    "11": -Fraction(117, 4),
    "13A": Fraction(243, 16),
    "13B": -Fraction(45, 16),
}


def _word(label):
    return "".join(family.axis_at(label, site, N) for site in range(N))


def _translate(word, outer_shift, inner_shift):
    """Literal permutation, shifting both base-five coordinates independently."""
    result = ["I"] * N
    for site, local_axis in enumerate(word):
        block, inner = divmod(site, 5)
        target = 5 * ((block + outer_shift) % 5) + (inner + inner_shift) % 5
        result[target] = local_axis
    return "".join(result)


def _cycle_axes(word, repetitions):
    cycle = {"I": "I", "X": "Y", "Y": "Z", "Z": "X"}
    for _ in range(repetitions):
        word = "".join(cycle[axis] for axis in word)
    return word


def _classify(word):
    """Classify from actual block weights, independently of the seed metadata."""
    block_weights = tuple(sum(axis != "I" for axis in word[5 * block:5 * block + 5])
                          for block in range(5))
    nonzero = sorted(w for w in block_weights if w)
    classification = {
        (3, 3, 3): "9",
        (3, 3, 5): "11",
        (3, 5, 5): "13A",
        (3, 3, 3, 4): "13B",
    }
    assert tuple(nonzero) in classification, block_weights
    return classification[tuple(nonzero)]


def _errors():
    answer = [0]
    for site in range(N):
        for axis in "XYZ":
            answer.append(family.pauli("I" * site + axis + "I" * (N - site - 1)))
    for left, right in combinations(range(N), 2):
        for a, b in product("XYZ", repeat=2):
            word = ["I"] * N
            word[left], word[right] = a, b
            answer.append(family.pauli("".join(word)))
    assert len(answer) == 1 + 3 * N + 9 * N * (N - 1) // 2 == 2776
    assert len(set(answer)) == len(answer)
    return answer


def check_sparse_construction():
    """Verify all three sparse logical axes and return a deterministic record.

    Every term's signed logical-coset membership is checked by phase-aware
    elimination of the actual signed stabilizer generators. Since each term
    preserves every Pauli error image, integer conjugation sums certify the
    full physical-space identities, rather than only projected energies.
    """
    specification = json.loads((HERE / "data" / "two_error_sparse_orbits.json").read_text())
    assert specification["n"] == N and specification["axis"] == "Z"
    assert Fraction(specification["J"]) == Fraction(521, 8)
    assert len(specification["orbits"]) == 103
    generators = family.code_generators(2)
    basis = family.signed_basis(generators, N)
    assert len(basis) == len(generators) == N - 1
    errors = _errors()
    records = []
    for repetitions, axis in enumerate("ZXY"):
        terms = {}
        classes = {}
        for entry in specification["orbits"]:
            seed = _cycle_axes(entry["seed"], repetitions)
            assert len(seed) == N
            classification = _classify(seed)
            assert classification == entry["class"]
            coefficient = Fraction(entry["orbit_mass"]) / 25
            assert (coefficient > 0) == (classification in ("9", "13A"))
            claimed_sign = entry["logical_coset_sign"]
            assert claimed_sign in (-1, 1)
            orbit = {family.pauli(_translate(seed, outer, inner))
                     for outer, inner in product(range(5), repeat=2)}
            assert len(orbit) == 25
            for label in orbit:
                assert label not in terms, "Stored seed orbits overlap."
                assert _classify(_word(label)) == classification
                terms[label] = coefficient, claimed_sign
                classes[label] = classification
        assert len(terms) == 2575
        logical = family.pauli(axis * N)
        for label, (coefficient, claimed_sign) in terms.items():
            difference, phase = family.multiply(label, logical, N)
            phase = (phase + family.stabilizer_phase(difference, basis, N)) % 4
            assert phase % 2 == 0 and (-1) ** (phase // 2) == claimed_sign
            assert family.weight(label, N) in (9, 11, 13)
        masses = {classification: sum(coefficient for label, (coefficient, _) in terms.items()
                                      if classes[label] == classification)
                  for classification in EXPECTED_MASSES}
        assert masses == EXPECTED_MASSES
        denominator = lcm(*(coefficient.denominator for coefficient, _ in terms.values()))
        labels = tuple(terms)
        x = tuple(label & MASK for label in labels)
        z = tuple(label >> N for label in labels)
        numerators = tuple(int(terms[label][0] * denominator) for label in labels)
        assert sum(numerators) == denominator
        assert sum(abs(value) for value in numerators) == Fraction(521, 8) * denominator
        for error in errors:
            ex, ez = error & MASK, error >> N
            actual = sum(-coefficient if (((qx & ez) ^ (qz & ex)).bit_count() & 1)
                         else coefficient for qx, qz, coefficient in zip(x, z, numerators))
            assert actual == denominator, (axis, error, Fraction(actual, denominator))
        records.append({
            "logical_axis": axis,
            "seed_orbits": len(specification["orbits"]),
            "orbit_size": 25,
            "native_pauli_terms": len(terms),
            "interaction_order": max(family.weight(label, N) for label in terms),
            "coupling_budget_per_logical_half_splitting": "521/8",
            "signed_coset_memberships_checked": len(terms),
            "full_error_image_intertwiners_checked": len(errors),
            "class_masses": {name: str(mass) for name, mass in masses.items()},
            "common_coefficient_denominator": str(denominator),
            "coefficient_denominator_bit_length": denominator.bit_length(),
        })
    return {"checks_passed": True, "axes": records}


if __name__ == "__main__":
    import sys
    if sys.flags.optimize:
        raise SystemExit("Assertions are part of the certificates; do not use -O.")
    print(json.dumps(check_sparse_construction(), indent=2, sort_keys=True))
