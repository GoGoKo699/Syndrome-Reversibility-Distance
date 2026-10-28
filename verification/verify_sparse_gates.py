#!/usr/bin/env python3
"""Exact sparse-gate certificates, using only Python's standard library.

The closed twenty-five-qubit example is reconstructed from three orbit seeds.
The separate recursive engine is checked at depths one, two, and three for all
three logical axes. Signed normalizer membership, every physical error-image
identity, and the code's coordinate-translation symmetries are checked without
assuming that the sparse engine's coefficient formulas are correct.

These finite checks do not exhaustively enumerate large-code distances. The
distance, infinite syndrome reversibility, and arbitrary-depth resource claims
depend on the written inductive proofs.
"""

from fractions import Fraction
from itertools import product
from math import lcm
from pathlib import Path
import json
import sys

import sparse_pauli as sparse
import verify_readout_safe_family as family


HERE = Path(__file__).resolve().parent


def literal_word(label, n):
    return "".join(family.axis_at(label, site, n) for site in range(n))


def translate(label, level, shifts):
    """Independent literal site permutation by base-five coordinate shifts."""
    n = 5 ** level
    assert len(shifts) == level
    result = ["I"] * n
    for site in range(n):
        target = 0
        quotient = site
        for coordinate, shift in enumerate(shifts):
            quotient, digit = divmod(quotient, 5)
            target += ((digit + shift) % 5) * 5 ** coordinate
        result[target] = family.axis_at(label, site, n)
    return family.pauli("".join(result))


def check_code_symmetry(level):
    n = 5 ** level
    generators = family.code_generators(level)
    assert all(not family.anticommutes(left, right, n)
               for left in generators for right in generators)
    basis = family.signed_basis(generators, n)
    assert len(basis) == len(generators) == n - 1
    assert all(family.weight(generator, n) % 2 == 0 for generator in generators)
    # It suffices to check the m generators of the coordinate-translation
    # group. Products then preserve the same signed stabilizer group.
    for coordinate in range(level):
        shifts = tuple(int(index == coordinate) for index in range(level))
        for generator in generators:
            moved = translate(generator, level, shifts)
            assert family.stabilizer_phase(moved, basis, n) == 0
        for axis in "XYZ":
            logical = family.pauli(axis * n)
            assert translate(logical, level, shifts) == logical
            assert all(not family.anticommutes(logical, g, n) for g in generators)
            assert family.weight(logical, n) % 2 == 1
    return basis


def check_hamiltonian(level, axis, terms, basis):
    """Certify full intertwiners from signed coset membership and conjugation.

    The mapping is physical label -> (positive-logical-coset sign, coefficient
    relative to that signed representative). The bare coefficient is their
    product. Membership ensures each individual term preserves every Pauli
    error image; the exact conjugation sums determine its logical rate.
    """
    n, d = 5 ** level, 3 ** level
    logical = family.pauli(axis * n)
    assert terms and all(coefficient for _, coefficient in terms.values())
    for label, (claimed_sign, coefficient) in terms.items():
        assert claimed_sign in (-1, 1) and isinstance(coefficient, Fraction)
        difference, phase = family.multiply(label, logical, n)
        phase = (phase + family.stabilizer_phase(difference, basis, n)) % 4
        assert phase % 2 == 0 and (-1) ** (phase // 2) == claimed_sign
        w = family.weight(label, n)
        assert w in (d, d + 2)
        assert (coefficient > 0) == (w == d)
    assert sum(coefficient for _, coefficient in terms.values()) == 1
    masses = {w: sum(coefficient for label, (_, coefficient) in terms.items()
                     if family.weight(label, n) == w) for w in (d, d + 2)}
    assert masses == {d: Fraction(d + 2, 2), d + 2: -Fraction(d, 2)}
    # All three local-label incidence sums vanish independently at each site.
    for site in range(n):
        for axis_ in "XYZ":
            assert sum(coefficient for label, (_, coefficient) in terms.items()
                       if family.axis_at(label, site, n) == axis_) == 0
    # Directly evaluate every no-error/single-Pauli logical rate as a separate
    # exact physical conjugation check, rather than using the expected formula.
    for site in range(n):
        for error_axis in "XYZ":
            error = family.pauli("I" * site + error_axis + "I" * (n - site - 1))
            assert sum(coefficient * (-1) ** family.anticommutes(label, error, n)
                       for label, (_, coefficient) in terms.items()) == 1
    budget = sum(abs(coefficient) for _, coefficient in terms.values())
    assert budget == d + 1
    assert budget == Fraction((d + 2) + d, (d + 2) - d)
    denominators = [coefficient.denominator for _, coefficient in terms.values()]
    common_denominator = lcm(*denominators)
    return {
        "axis": axis,
        "native_pauli_terms": len(terms),
        "interaction_order": max(family.weight(label, n) for label in terms),
        "coupling_budget_per_logical_half_splitting": str(budget),
        "logical_weight_masses": {str(w): str(mass) for w, mass in masses.items()},
        "signed_coset_memberships_checked": len(terms),
        "full_error_image_intertwiners_checked": 1 + 3 * n,
        "common_coefficient_denominator": str(common_denominator),
        "largest_denominator_bit_length": max(denominators).bit_length(),
    }


def closed_example():
    specification = json.loads((HERE / "data" / "sparse_25_qubit.json").read_text())
    assert specification["physical_qubits"] == 25
    assert specification["logical_axis"] == "Z"
    assert len(specification["orbits"]) == 3
    cycle = specification["axis_cycle"]
    assert cycle == {"I": "I", "X": "Y", "Y": "Z", "Z": "X"}
    basis = check_code_symmetry(2)
    records = []
    for axis in "ZXY":
        terms = {}
        orbit_records = []
        for entry in specification["orbits"]:
            seed = entry["seed"]
            for _ in range("ZXY".index(axis)):
                seed = "".join(cycle[letter] for letter in seed)
            label = family.pauli(seed)
            assert family.weight(label, 25) == entry["weight"]
            orbit = {translate(label, 2, shifts)
                     for shifts in product(range(5), repeat=2)}
            assert len(orbit) == 25
            assert orbit == set(sparse.translation_orbit(label, 2))
            sign = entry["logical_coset_sign"]
            bare_coefficient = Fraction(entry["bare_coefficient_per_orbit_word"])
            assert bare_coefficient > 0
            for physical_word in orbit:
                assert physical_word not in terms
                terms[physical_word] = sign, sign * bare_coefficient
            orbit_records.append({
                "seed": seed,
                "weight": entry["weight"],
                "orbit_size": len(orbit),
                "bare_coefficient_per_term": str(bare_coefficient),
                "logical_coset_sign": sign,
            })
        assert len(terms) == 75
        assert {sign * coefficient for sign, coefficient in terms.values()} == {
            Fraction(2, 25), Fraction(7, 50), Fraction(9, 50)}
        result = check_hamiltonian(2, axis, terms, basis)
        assert result["common_coefficient_denominator"] == "50"
        result["orbits"] = orbit_records
        records.append(result)
    return records


def recursive_examples():
    cache = {}
    levels = []
    for level in (1, 2, 3):
        n, d = 5 ** level, 3 ** level
        basis = check_code_symmetry(level)
        classes = sparse.sparse_seed_classes(level, cache)
        records = []
        for axis in "XYZ":
            for w, seeds in classes[axis].items():
                assert w in (d, d + 2) and len(seeds) <= 3
                assert sum(mass for _, mass in seeds.values()) == 1
                assert all(mass > 0 for _, mass in seeds.values())
                for label_axis in "XYZ":
                    assert sum(sum(family.axis_at(label, site, n) == label_axis
                                   for site in range(n)) * mass
                               for label, (_, mass) in seeds.items()) == Fraction(w, 3)
            seed_hamiltonian = sparse.sparse_seed_hamiltonian(level, axis, cache)
            assert len(seed_hamiltonian) <= 4
            terms = sparse.sparse_hamiltonian(level, axis, cache)
            assert len(terms) <= 4 * n
            # Verify each engine orbit by the separate site-permutation code.
            for label in seed_hamiltonian:
                assert set(sparse.translation_orbit(label, level)) == {
                    translate(label, level, shifts)
                    for shifts in product(range(5), repeat=level)}
            result = check_hamiltonian(level, axis, terms, basis)
            assert int(result["common_coefficient_denominator"]) <= 16 * n ** 4
            result["seed_orbits"] = len(seed_hamiltonian)
            result["four_n_term_bound"] = 4 * n
            result["coefficient_denominator_bound"] = str(16 * n ** 4)
            records.append(result)
        levels.append({
            "depth": level,
            "physical_qubits": n,
            "distance_from_written_concatenation_proof": d,
            "signed_stabilizer_rank": n - 1,
            "coordinate_shift_generators_checked": level,
            "even_stabilizer_odd_logical_parity_checked": True,
            "axes": records,
        })
    return levels


def main():
    if sys.version_info < (3, 10):
        raise SystemExit("Python 3.10 or later is required.")
    if sys.flags.optimize:
        raise SystemExit("Assertions are part of the certificates; do not use -O.")
    native_compressed = sparse.incidence_compressed_hamiltonian(2, "Z")
    assert len(native_compressed) <= 76
    native_record = check_hamiltonian(2, "Z", native_compressed, check_code_symmetry(2))
    native_record["three_n_plus_one_term_bound"] = 76
    print(json.dumps({
        "checks_passed": True,
        "closed_twenty_five_qubit_construction": closed_example(),
        "recursive_sparse_constructions": recursive_examples(),
        "optional_native_incidence_reduction": native_record,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
