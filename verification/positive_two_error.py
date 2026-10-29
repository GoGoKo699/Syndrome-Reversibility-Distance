"""Exact positive depth-two classes for response-preserving sparse recursion.

The committed data stores positive orbit masses for logical Z. Coordinatewise
C5 x C5 translation averaging and a global cyclic axis relabeling supply each
of X,Y,Z. The loader expands native words; no unrecorded numerical solver is
used. The checker independently verifies signs, class patterns and complete
weight-at-most-two response profiles using integer orbit correlations.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
import json
from math import lcm
from pathlib import Path
import sys

from sparse_pauli import translation_orbit
from verify_readout_safe_family import (
    axis_at, code_generators, multiply, pauli, signed_basis, stabilizer_phase, weight,
)

CLASS_NAMES = ("A0", "A1", "A2", "F")
RESPONSES = {
    "identity": (Fraction(1),) * 4,
    "single": (Fraction(13, 25), Fraction(31, 75), Fraction(23, 75), Fraction(23, 75)),
    "same_block": (Fraction(9, 25), Fraction(89, 225), Fraction(97, 225), Fraction(11, 75)),
    "different_blocks": (Fraction(29, 125), Fraction(13, 125), Fraction(1, 225), Fraction(23, 375)),
}


def _cycle_axes(label, n):
    """Physical X->Y->Z->X, with no Hermitian sign changes."""
    mask = (1 << n) - 1
    x, z = label & mask, label >> n
    return (x ^ z) | (x << n)


@lru_cache(maxsize=None)
def seed_classes(axis="Z"):
    """Return class -> immutable (packed Pauli, coset sign, orbit mass) atoms."""
    assert axis in "XYZ"
    data = json.loads((Path(__file__).resolve().parent / "data" /
                       "two_error_positive_classes.json").read_text())
    assert data["logical_axis"] == "Z" and data["features"] == 112
    cycles = {"Z": 0, "X": 1, "Y": 2}[axis]
    answer = {}
    for name in CLASS_NAMES:
        atoms = []
        for record in data["classes"][name]:
            label = pauli(record["pauli"])
            for _ in range(cycles):
                label = _cycle_axes(label, 25)
            sign, mass = record["coset_sign"], Fraction(record["mass"])
            assert sign in (-1, 1) and mass > 0
            assert weight(label, 25) == record["weight"]
            atoms.append((label, sign, mass))
        assert len(atoms) <= 112 and sum(mass for _, _, mass in atoms) == 1
        answer[name] = tuple(atoms)
    return answer


@lru_cache(maxsize=None)
def _expanded(axis):
    result = {}
    for name, atoms in seed_classes(axis).items():
        words = {}
        for label, sign, mass in atoms:
            for translated in translation_orbit(label, 2):
                if translated in words:
                    old_sign, old_mass = words[translated]
                    assert old_sign == sign
                    words[translated] = sign, old_mass + mass / 25
                else:
                    words[translated] = sign, mass / 25
        assert sum(mass for _, mass in words.values()) == 1
        result[name] = tuple((label, sign, mass) for label, (sign, mass)
                             in sorted(words.items()))
    return result


def load_classes(level=2, axis="Z"):
    """Return positive normalized native lists, class -> [(label,sign,mass)].

    This loader supplies the depth-two base of the general sparse recursion.
    All atoms are actual physical Pauli words, not average operators.
    """
    assert level == 2 and axis in "XYZ"
    return {name: list(atoms) for name, atoms in _expanded(axis).items()}


def _negative(delta):
    return ((-(delta % 5)) % 5) + 5 * ((-(delta // 5)) % 5)


DELTAS = tuple(delta for delta in range(1, 25) if delta < _negative(delta))
FEATURES = (("identity",),) + tuple(("single", a) for a in range(3)) + tuple(
    ("pair", delta, a, b) for delta in DELTAS for a in range(3) for b in range(3))
assert len(FEATURES) == 112
FEATURE_INDEX = {feature: index for index, feature in enumerate(FEATURES)}


def _integer_orbit_column(label):
    """Sum all 25 translated conjugation characters without FFT/arithmetic tolerances."""
    axes = [axis_at(label, site, 25) for site in range(25)]
    character = [[1 if local in ("I", error_axis) else -1 for local in axes]
                 for error_axis in "XYZ"]
    answer = [25] + [sum(row) for row in character]
    for delta in DELTAS:
        shifted = [((site % 5 + delta % 5) % 5) +
                   5 * ((site // 5 + delta // 5) % 5) for site in range(25)]
        for a in range(3):
            for b in range(3):
                answer.append(sum(character[a][site] * character[b][shifted[site]]
                                  for site in range(25)))
    return tuple(answer)


def _pattern(label):
    block_weights = [sum(axis_at(label, 5 * block + j, 25) != "I" for j in range(5))
                     for block in range(5)]
    return tuple(sorted(block_weights))


PATTERNS = {"A0": (0, 0, 3, 3, 3), "A1": (0, 0, 3, 3, 5),
            "A2": (0, 0, 3, 5, 5), "F": (0, 3, 3, 3, 4)}


def _response_index(first, a, second, b):
    delta = ((second % 5 - first % 5) % 5) + 5 * ((second // 5 - first // 5) % 5)
    if delta > _negative(delta):
        delta, a, b = _negative(delta), b, a
    return FEATURE_INDEX["pair", delta, a, b]


def check_positive_classes():
    generators = code_generators(2)
    basis = signed_basis(generators, 25)
    reports = []
    for axis in "XYZ":
        logical = pauli(axis * 25)
        native = load_classes(axis=axis)
        classes = {}
        for number, name in enumerate(CLASS_NAMES):
            atoms = seed_classes(axis)[name]
            denominator = lcm(*(mass.denominator for _, _, mass in atoms))
            numerators = [int(mass * denominator) for _, _, mass in atoms]
            totals = [0] * 112
            for (label, _, _), numerator in zip(atoms, numerators):
                assert _pattern(label) == PATTERNS[name]
                column = _integer_orbit_column(label)
                totals = [old + numerator * value for old, value in zip(totals, column)]
            response = tuple(Fraction(total, 25 * denominator) for total in totals)
            # Independently verify every expanded native term, including actual
            # Hermitian stabilizer phases rather than only label membership.
            for label, sign, mass in native[name]:
                assert mass > 0 and _pattern(label) == PATTERNS[name]
                difference, phase = multiply(label, logical, 25)
                eigenphase = (phase + stabilizer_phase(difference, basis, 25)) % 4
                assert eigenphase % 2 == 0 and (-1) ** (eigenphase // 2) == sign
            counts = {key: 0 for key in RESPONSES}
            assert response[0] == RESPONSES["identity"][number]
            counts["identity"] += 1
            for site in range(25):
                for a in range(3):
                    assert response[FEATURE_INDEX["single", a]] == RESPONSES["single"][number]
                    counts["single"] += 1
            for first, second in combinations(range(25), 2):
                kind = "same_block" if first // 5 == second // 5 else "different_blocks"
                for a, b in product(range(3), repeat=2):
                    assert response[_response_index(first, a, second, b)] == RESPONSES[kind][number]
                    counts[kind] += 1
            assert counts == {"identity": 1, "single": 75,
                              "same_block": 450, "different_blocks": 2250}
            classes[name] = {
                "orbit_seeds": len(atoms), "native_terms": len(native[name]),
                "orbit_coefficient_denominator_bits": denominator.bit_length(),
                "full_response_equations_checked": sum(counts.values()),
                "response_by_error_pattern": {key: str(values[number])
                                               for key, values in RESPONSES.items()},
                "signed_native_memberships_checked": len(native[name]),
            }
        reports.append({"logical_axis": axis, "classes": classes})
    return {"positive_class_profiles_verified": True, "n": 25,
            "orbit_response_coordinates": len(FEATURES), "axes": reports}


if __name__ == "__main__":
    if sys.flags.optimize:
        raise SystemExit("Assertions are scientific checks; do not use -O.")
    print(json.dumps(check_positive_classes(), indent=2, sort_keys=True))
