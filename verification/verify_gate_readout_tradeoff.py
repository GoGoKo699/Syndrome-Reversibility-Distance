#!/usr/bin/env python3
"""Exact finite certificates for the gate-locality/readout tradeoff.

The same standalone, standard-library program is kept in both companion
repositories. It enumerates the full phase-free Pauli centralizer and computes
syndrome reversibility distance directly from equal-support collisions, rather
than assuming the joint theorem. Strings list physical sites from left to right.

Code and Hamiltonian inputs are the explicit constructions in:
https://github.com/GoGoKo699/Error-Transparent-Locality/blob/main/science/four_body_construction.md
https://github.com/GoGoKo699/Syndrome-Reversibility-Distance/blob/main/science/code_families.md

This is a finite-example certificate, not a computational proof of the general
odd-distance theorem. Run without Python's -O option: assertions are checks.
"""

from fractions import Fraction
from math import comb
import json
import sys


NINE_GENERATORS = (
    "IIIIXIIIX", "XXIXIXIII", "XIXXIIXII", "IXXXXIIXI",
    "ZIIZIIIZI", "IZIIIZIZI", "IIZIIIZZI", "IZZZZIIIZ",
)
FIVE_GENERATORS = ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
NINE_Z = "IIZZIZIII"
NINE_TERMS = (
    ("IIZZIZIII", "1"), ("ZIIIIZZII", "1"),
    ("IIIIZIIZZ", "3/2"), ("ZZZIIIIII", "1/2"),
    ("IZIIZZIIZ", "-1/2"), ("ZIZIIZIZI", "-1"),
    ("IIIZIZZZI", "-1/2"), ("ZIIZZIIIZ", "-1/2"),
    ("IIZIZIZIZ", "-1/2"),
)


def pauli(word):
    n = len(word)
    x = sum((axis in "XY") << site for site, axis in enumerate(word))
    z = sum((axis in "ZY") << site for site, axis in enumerate(word))
    assert set(word) <= set("IXYZ")
    return x | (z << n)


def support(label, n):
    return (label & ((1 << n) - 1)) | (label >> n)


def weight(label, n):
    return support(label, n).bit_count()


def anticommutes(left, right, n):
    mask = (1 << n) - 1
    return (((left & mask) & (right >> n)) ^
            ((left >> n) & (right & mask))).bit_count() % 2


def span(generators):
    values = {0}
    for generator in generators:
        values |= {value ^ generator for value in values}
    return values


def syndrome(label, generators, n):
    return sum(anticommutes(label, generator, n) << index
               for index, generator in enumerate(generators))


def analyze_code(words):
    """Exhaustively compute distance and the defining collision minimum."""
    n = len(words[0])
    assert all(len(word) == n for word in words)
    generators = tuple(map(pauli, words))
    assert all(not anticommutes(left, right, n)
               for left in generators for right in generators)
    stabilizers = span(generators)
    assert len(stabilizers) == 1 << len(generators)
    normalizer = {label for label in range(1 << (2 * n))
                  if not syndrome(label, generators, n)}
    assert len(normalizer) == 1 << (2 * n - len(generators))
    logicals = normalizer - stabilizers
    distance = min(weight(label, n) for label in logicals)
    by_support = {}
    for logical in logicals:
        by_support.setdefault(support(logical, n), []).append(logical)
    collision_weights = []
    for stabilizer in stabilizers:
        for logical in by_support.get(support(stabilizer, n), ()):
            collision_weights.append(weight(stabilizer ^ logical, n))
    assert all(value >= distance and value % 2 == 0
               for value in collision_weights)
    result = {
        "n": n,
        "k": n - len(generators),
        "distance": distance,
        "stabilizer_count": len(stabilizers),
        "centralizer_count": len(normalizer),
        "support_collision_count": len(collision_weights),
        "syndrome_reversibility_distance": min(collision_weights)
        if collision_weights else "infinity",
    }
    return result, stabilizers, normalizer, generators


def subwords(label, n):
    active = [site for site in range(n) if support(label, n) >> site & 1]
    words = []
    for selected in range(1 << len(active)):
        sites = sum(1 << site for index, site in enumerate(active)
                    if selected >> index & 1)
        words.append(label & (sites | (sites << n)))
    return words


def literal_product(left, right):
    """Return phase i**power and literal Hermitian word for a Pauli product."""
    forward = {("X", "Y"): "Z", ("Y", "Z"): "X", ("Z", "X"): "Y"}
    power = 0
    result = []
    for a, b in zip(left, right):
        if a == "I":
            result.append(b)
        elif b == "I":
            result.append(a)
        elif a == b:
            result.append("I")
        elif (a, b) in forward:
            power += 1
            result.append(forward[a, b])
        else:
            power -= 1
            result.append(forward[b, a])
    return power % 4, "".join(result)


def check_nine_witness(stabilizers, normalizer, generators):
    n = 9
    e_word, g_word, l_word = "ZIZIIZIZI", "XIXIXXIXI", "YIYIXYIYI"
    e, g, logical = map(pauli, (e_word, g_word, l_word))
    assert g in stabilizers and logical in normalizer - stabilizers
    assert e in normalizer - stabilizers
    assert e ^ pauli(NINE_Z) in stabilizers
    assert g ^ logical == e and weight(e, n) == 4
    assert support(g, n) == support(logical, n)
    assert literal_product(g_word, e_word) == (0, l_word)
    pieces = subwords(e, n)
    proper = [piece for piece in pieces if piece not in (0, e)]
    assert len(proper) == 14
    assert all(piece not in normalizer for piece in proper)
    grouped = {}
    for piece in pieces:
        grouped.setdefault(syndrome(piece, generators, n), []).append(piece)
    assert len(grouped) == 8
    assert all(len(pair) == 2 and pair[0] ^ pair[1] == e
               for pair in grouped.values())
    return {
        "logical_circuit": e_word,
        "stabilizer_witness": g_word,
        "equal_support_logical": l_word,
        "literal_product_phase": "+1",
        "logical_class": "Z",
        "proper_subwords_with_nonzero_syndrome": len(proper),
        "complementary_syndrome_pairs": len(grouped),
    }


def check_nine_hamiltonian(stabilizers):
    n = 9
    logical_z = pauli(NINE_Z)
    positive_z_stabilizers = span(tuple(map(pauli, NINE_GENERATORS[4:])))
    assert positive_z_stabilizers <= stabilizers
    terms = [(pauli(word), Fraction(coefficient)) for word, coefficient in NINE_TERMS]
    # The differences are pure-Z stabilizer products, so no hidden signs occur.
    assert all(label ^ logical_z in positive_z_stabilizers for label, _ in terms)
    assert sum(coefficient for _, coefficient in terms) == 1
    for site in range(n):
        assert sum(coefficient for label, coefficient in terms
                   if support(label, n) >> site & 1) == 0
        for axis in "XYZ":
            error = pauli("I" * site + axis + "I" * (n - site - 1))
            assert sum((-1) ** anticommutes(label, error, n) * coefficient
                       for label, coefficient in terms) == 1
    masses = {w: sum(coefficient for label, coefficient in terms
                     if weight(label, n) == w) for w in (3, 4)}
    assert masses == {3: 4, 4: -3}
    budget = sum(abs(coefficient) for _, coefficient in terms)
    assert budget == 7
    return {
        "maximum_pauli_weight": max(weight(label, n) for label, _ in terms),
        "logical_weight_masses": {str(w): str(value) for w, value in masses.items()},
        "coupling_budget_per_logical_half_splitting": str(budget),
        "single_pauli_error_identities": 27,
    }


def check_repetition_polynomial():
    # Complementary error subsets produce half of the weight-two probability
    # plus the weight-three and weight-four binomial tail for 0 <= q <= 1/2.
    coefficients = [Fraction(0) for _ in range(5)]
    for weight_, factor in ((2, Fraction(1, 2)), (3, Fraction(1)), (4, Fraction(1))):
        for exponent in range(5 - weight_):
            coefficients[weight_ + exponent] += (
                factor * comb(4, weight_) * comb(4 - weight_, exponent) * (-1) ** exponent)
    assert coefficients == [0, 0, 3, -2, 0]
    return {"q_polynomial_coefficients": [int(value) for value in coefficients],
            "theta_fourth_order_coefficient": str(Fraction(3, 16))}


def main():
    if sys.version_info < (3, 10):
        raise SystemExit("Python 3.10 or later is required.")
    if sys.flags.optimize:
        raise SystemExit("Assertions are part of the certificates; do not use -O.")
    nine, nine_s, nine_n, nine_g = analyze_code(NINE_GENERATORS)
    five, _, _, _ = analyze_code(FIVE_GENERATORS)
    padded_words = tuple(word + "I" for word in FIVE_GENERATORS) + ("IIIIIZ",)
    padded, padded_s, padded_n, _ = analyze_code(padded_words)
    assert (nine["distance"], nine["syndrome_reversibility_distance"]) == (3, 4)
    assert (five["distance"], five["syndrome_reversibility_distance"]) == (3, "infinity")
    assert (padded["distance"], padded["syndrome_reversibility_distance"]) == (3, "infinity")
    padded_logical = "IYYIXZ"
    label = pauli(padded_logical)
    fixed = pauli("IIIIIZ")
    assert label in padded_n - padded_s and weight(label, 6) == 4
    assert fixed in padded_s and label ^ fixed in padded_n - padded_s
    assert fixed in subwords(label, 6)
    padded["weight_four_logical"] = padded_logical
    padded["proper_stabilizer_subword"] = "IIIIIZ"
    padded["weight_four_word_is_a_circuit"] = False
    result = {
        "checks_passed": True,
        "nine_qubit": nine,
        "five_qubit": five,
        "five_qubit_with_fixed_ancilla": padded,
        "nine_qubit_support_witness": check_nine_witness(nine_s, nine_n, nine_g),
        "nine_qubit_hamiltonian": check_nine_hamiltonian(nine_s),
        "four_site_rotation_error": check_repetition_polynomial(),
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
