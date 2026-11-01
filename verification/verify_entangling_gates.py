#!/usr/bin/env python3
"""Exact certificates for protected logical Pauli strings and entanglers.

All 540,376 distinct physical Pauli errors through weight three on two
25-qubit blocks are enumerated. Axis-independent child character profiles
then check every pair of prescribed logical Pauli axes. The large native
Hamiltonian is specified by exact tensor factors, not expanded in full.
Only Python's standard library is required; assertions remain enabled.
"""

from collections import Counter
from fractions import Fraction as F
from itertools import combinations, islice, product
from math import comb, lcm
import json
import sys

import three_error_recursion as child
import verify_three_error_gates as base
import verify_readout_safe_family as family


PROFILE_NAMES = tuple(base.PROFILES)
PROFILE_INDEX = {base.PROFILES[name]: i for i, name in enumerate(PROFILE_NAMES)}


def _label(locations, n):
    x, z = 0, 0
    for site, axis in locations:
        if axis in "XY":
            x |= 1 << site
        if axis in "YZ":
            z |= 1 << site
    return x | (z << n)


def _independent_profile(signature, configurations):
    """Average integer inner characters without the base verifier's denominator trick."""
    active = tuple(i for i, values in enumerate(signature) if values != (30,) * 8)
    denominator = 30 ** len(active)
    result = []
    for configs in configurations:
        total = 0
        for config in configs:
            value = 1
            for block in active:
                value *= signature[block][base.OPTION_INDEX[config[block]]]
            total += value
        result.append(F(total, denominator * len(configs)))
    return tuple(result)


def checked_child_lookup():
    """Map every actual child error to a profile, proving independence of logical axis."""
    group, cosets = base.signed_cosets()
    assert len(group) == 16
    table = base.character_tables(cosets)
    configs = {axis: base.configurations(axis, cosets) for axis in "XYZ"}
    lookup, cache = {}, {}
    by_weight = Counter()
    for size in range(4):
        for sites in combinations(range(25), size):
            for axes in product("XYZ", repeat=size):
                locations = tuple(zip(sites, axes))
                word = ["I"] * 25
                for site, axis in locations:
                    word[site] = axis
                signature = tuple(table["".join(word[5 * block:5 * block + 5])]
                                  for block in range(5))
                if signature not in cache:
                    values = tuple(_independent_profile(signature, configs[axis])
                                   for axis in "XYZ")
                    assert values[0] == values[1] == values[2]
                    assert values[0] in PROFILE_INDEX
                    cache[signature] = PROFILE_INDEX[values[0]]
                label = _label(locations, 25)
                assert label not in lookup and family.weight(label, 25) == size
                lookup[label] = cache[signature]
                by_weight[cache[signature], size] += 1
    assert len(lookup) == 64876
    assert sum(by_weight.values()) == 64876
    for index, name in enumerate(PROFILE_NAMES):
        assert sum(count for (profile, _), count in by_weight.items()
                   if profile == index) == base.EXPECTED_PROFILE_COUNTS[name]
    return lookup, {"physical_child_errors_checked": len(lookup),
                    "logical_axes": ["X", "Y", "Z"],
                    "axis_independent_character_signatures": len(cache),
                    "profiles": len(PROFILE_NAMES)}


def _base_increments():
    gates = ((F(1), F(0), F(0), F(0), F(0), F(0), F(0)),
             tuple(map(F, ("11/2", "-9/2", "0", "0", "0", "0", "0"))),
             tuple(map(F, ("143/8", "-117/4", "243/16", "0", "-45/16", "0", "0"))),
             base.MASSES)
    return (gates[0],) + tuple(tuple(a - b for a, b in zip(gates[t], gates[t - 1]))
                              for t in range(1, 4))


def check_all_fifty_qubit_errors():
    lookup, child_report = checked_child_lookup()
    coefficients = _base_increments()
    responses = tuple(tuple(sum(a * b for a, b in zip(row, base.PROFILES[name]))
                            for row in coefficients) for name in PROFILE_NAMES)
    pair_responses = {(left, right): sum(responses[left][a] * responses[right][b]
                                         for a in range(4) for b in range(4) if a + b <= 3)
                      for left in range(7) for right in range(7)}
    seen, profiles, weight_counts = set(), Counter(), Counter()
    side_weights = Counter()
    mask = (1 << 25) - 1
    for size in range(4):
        for sites in combinations(range(50), size):
            for axes in product("XYZ", repeat=size):
                locations = tuple(zip(sites, axes))
                label = _label(locations, 50)
                assert label not in seen and family.weight(label, 50) == size
                seen.add(label)
                x, z = label & ((1 << 50) - 1), label >> 50
                left = (x & mask) | ((z & mask) << 25)
                right = (x >> 25) | ((z >> 25) << 25)
                key = lookup[left], lookup[right]
                assert pair_responses[key] == 1
                profiles[key] += 1
                weight_counts[size] += 1
                side_weights[family.weight(left, 25), family.weight(right, 25)] += 1
    assert len(seen) == sum(3 ** size * comb(50, size) for size in range(4)) == 540376
    assert weight_counts == {0: 1, 1: 150, 2: 11025, 3: 529200}
    assert len(profiles) == 18
    for (left_weight, right_weight), count in side_weights.items():
        assert count == 3 ** (left_weight + right_weight) * comb(25, left_weight) * comb(25, right_weight)
    return {"child_profile_certificate": child_report,
            "actual_distinct_physical_errors_enumerated": len(seen),
            "counts_by_physical_weight": {str(w): weight_counts[w] for w in range(4)},
            "ordered_child_profile_pairs": len(profiles),
            "logical_axis_products_certified": [a + b for a in "XYZ" for b in "XYZ"],
            "physical_errors_per_axis_product": len(seen),
            "axis_independence_reuses_the_complete_child_certificate": True,
            "profile_counts": [{"left": PROFILE_NAMES[left], "right": PROFILE_NAMES[right],
                                "physical_errors": count}
                               for (left, right), count in sorted(profiles.items())],
            "full_native_hamiltonian_expanded": False}


def _padd(*values):
    result = Counter()
    for value in values:
        result.update(value)
    return {key: coefficient for key, coefficient in result.items() if coefficient}


def _pmul(left, right):
    result = Counter()
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            result[i + k, j + ell] += a * b
    return {key: coefficient for key, coefficient in result.items() if coefficient}


def check_budget_semigroup():
    """Bivariate coefficient identity B(x,z)B(y,z)=B(x+y,z) mod z^4."""
    one_variable = ((F(1),), (F(0), F(1)), (F(0), F(13, 8), F(1, 2)),
                    (F(0), F(107, 24), F(13, 8), F(1, 6)))
    xp = tuple({(degree, 0): c for degree, c in enumerate(row) if c}
               for row in one_variable)
    yp = tuple({(0, degree): c for degree, c in enumerate(row) if c}
               for row in one_variable)
    for order in range(4):
        actual = _padd(*(_pmul(xp[a], yp[order - a]) for a in range(order + 1)))
        expected = {}
        for degree, coefficient in enumerate(one_variable[order]):
            for x_degree in range(degree + 1):
                key = x_degree, degree - x_degree
                expected[key] = coefficient * comb(degree, x_degree)
        expected = {key: coefficient for key, coefficient in expected.items() if coefficient}
        assert actual == expected
    examples = []
    for level in range(2, 6):
        d = 3 ** level
        for active in range(1, 6):
            effective = active * d
            increments = (F(1), F(d), F(d * (4 * d + 13), 8),
                          F(d * (4 * d * d + 39 * d + 107), 24))
            coefficients = (F(1), F(0), F(0), F(0))
            for _ in range(active):
                coefficients = tuple(sum(coefficients[a] * increments[order - a]
                                         for a in range(order + 1)) for order in range(4))
            expected = F(4 * effective ** 3 + 51 * effective ** 2 + 170 * effective + 24, 24)
            assert sum(coefficients) == expected
            examples.append({"level": level, "active_logical_blocks": active,
                             "target_coset_distance": effective,
                             "interaction_order": effective + 6,
                             "native_budget_per_target_half_splitting": str(expected)})
    assert examples[1]["native_budget_per_target_half_splitting"] == "1789"
    return {"bivariate_semigroup_checked_coefficientwise": True, "examples": examples}


def check_projected_lower_bound():
    from verify_growing_three_error import symbolic_certificates
    inherited = symbolic_certificates()
    assert inherited["coefficientwise_mass_and_dual_identities"]
    signs = (1, -1, 1, -1, -1, 1, -1)
    examples = []
    for effective in (9, 18, 27, 36, 54, 81, 108, 162):
        masses = (F((effective + 2) * (effective + 4) * (effective + 6), 48),
                  -F(effective * (effective + 4) * (effective + 6), 16),
                  F(effective * (2 * effective + 9) * (effective + 6), 32),
                  -F(effective * (2 * effective + 15) * (effective + 6), 96),
                  -F(5 * effective * (effective + 6), 32),
                  F(5 * effective * (effective + 6), 32), -F(effective, 6))
        columns = tuple((1, effective + a, (effective + a) ** 2, (effective + a) ** 3,
                         3 * effective + b, (effective + a) * (3 * effective + b),
                         9 * effective + c)
                        for a, b, c in ((0, 0, 0), (2, 16, 98), (4, 32, 196),
                                        (6, 48, 294), (4, 16, 64), (6, 32, 162),
                                        (6, 18, 54)))
        budget = F(4 * effective ** 3 + 51 * effective ** 2 + 170 * effective + 24, 24)
        dual = (budget, -F(12 * effective ** 2 + 111 * effective + 404, 24),
                F(effective + 5, 2), -F(1, 6), F(effective + 38, 8), -F(1, 8), -F(1, 2))
        assert tuple(sum(mass * column[row] for mass, column in zip(masses, columns))
                     for row in range(7)) == (1, 0, 0, 0, 0, 0, 0)
        assert tuple(sum(y * x for y, x in zip(dual, column)) for column in columns) == signs
        assert sum(map(abs, masses)) == budget
        if effective == 18:
            assert budget == 1789
        examples.append({"target_coset_distance": effective,
                         "projected_pattern_masses": list(map(str, masses)),
                         "native_budget_lower_bound": str(budget)})
    orthogonality = 0
    for target in (a + b for a in "XYZ" for b in "XYZ"):
        for logical in (a + b for a in "IXYZ" for b in "IXYZ"):
            word, phase = base.multiply(target, logical)
            if target == logical:
                assert word == "II" and phase == 0
            else:
                assert word != "II"
            orthogonality += 1
    return {"target_pauli_coefficient_isolated_from_other_logical_cosets": True,
            "generic_moment_mass_and_dual_identities_checked_coefficientwise": True,
            "two_logical_pauli_projection_orthogonality_checks": orthogonality,
            "ordinary_code_distance_used_as_target_distance": False,
            "moment_rows": ["1", "W", "W^2", "W^3", "R", "W*R", "S"],
            "examples": examples}


def check_target_configurations():
    """Classify every target-coset shell configuration through total weight24."""
    _, cosets = base.signed_cosets()
    configurations = {axis: base.configurations(axis, cosets) for axis in "XYZ"}
    patterns = ((6, 0, 0), (5, 0, 1), (4, 0, 2), (3, 0, 3),
                (6, 1, 0), (5, 1, 1), (8, 0, 0))
    expected = (100, 600, 1500, 2000, 400, 2400, 120)
    for axis in "XYZ":
        assert min(sum(w for _, w in configuration)
                   for group in configurations[axis] for configuration in group) == 9
        assert min(3 * base.weight(outer) for shell in cosets[axis].values()
                   for outer in shell) == 9
    for left_axis, right_axis in product("XYZ", repeat=2):
        actual = Counter()
        for left_group in configurations[left_axis]:
            for right_group in configurations[right_axis]:
                for left, right in product(left_group, right_group):
                    weights = tuple(w for _, w in left + right)
                    if sum(weights) > 24:
                        continue
                    assert 18 <= sum(weights) <= 24 and sum(weights) % 2 == 0
                    counts = Counter(weights)
                    pattern = tuple(counts[w] for w in (3, 4, 5))
                    assert pattern in patterns
                    actual[pattern] += 1
        assert tuple(actual[pattern] for pattern in patterns) == expected
    # Collect the truncated increment expansion into its disjoint base-class pairs.
    increments = _base_increments()
    pairs = {}
    for left, right in product(range(7), repeat=2):
        mass = sum(increments[a][left] * increments[b][right]
                   for a in range(4) for b in range(4) if a + b <= 3)
        if mass:
            pairs[left, right] = mass
    assert len(pairs) == 18
    native_count = sum(base.COUNTS[left] * base.COUNTS[right] for left, right in pairs)
    assert native_count == 41032000000
    signs = (1, -1, 1, -1, -1, 1, -1)
    bottom_patterns = ((3, 0, 0), (2, 0, 1), (1, 0, 2), (0, 0, 3),
                       (3, 1, 0), (2, 1, 1), (5, 0, 0))
    masses = [F(0)] * 7
    for (left, right), mass in pairs.items():
        pattern = tuple(a + b for a, b in zip(bottom_patterns[left], bottom_patterns[right]))
        index = patterns.index(pattern)
        assert mass * signs[index] > 0
        masses[index] += mass
    assert tuple(masses) == tuple(map(F, ("220", "-594", "1215/2", "-459/2",
                                         "-135/2", "135/2", "-3")))
    assert sum(masses) == 1 and sum(map(abs, pairs.values())) == 1789
    native_coefficients = tuple(mass / (base.COUNTS[left] * base.COUNTS[right])
                                for (left, right), mass in pairs.items())
    denominator = lcm(*(value.denominator for value in native_coefficients))
    assert denominator == 32000000000 and denominator.bit_length() == 35
    assert (3200000 ** 2) % denominator == 0
    return {"ordinary_product_code_distance": 9, "target_ZZ_coset_distance": 18,
            "target_logical_axis_products_classified": 9,
            "complete_allowed_shell_configurations_per_target": sum(expected),
            "configuration_counts_by_pattern": dict(zip(base.CLASS_NAMES, expected)),
            "signed_pattern_masses": dict(zip(base.CLASS_NAMES, map(str, masses))),
            "nonzero_disjoint_tensor_class_pairs": len(pairs),
            "dense_native_terms_specified_per_target": native_count,
            "dense_native_expansion_enumerated": False,
            "collected_native_common_denominator": str(denominator),
            "collected_native_common_denominator_bits": denominator.bit_length()}


def _gmul(a, b):
    return a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]


def check_universal_control():
    generators = {"XI", "YI", "ZI", "IX", "IY", "IZ", "ZZ"}
    closure = set(generators)
    while True:
        added = set()
        for left, right in combinations(sorted(closure), 2):
            if base.anticommutes(left, right):
                word, phase = base.multiply(left, right)
                assert phase % 2 == 1 and word != "II"
                added.add(word)
        if added <= closure:
            break
        closure.update(added)
    assert closure == {a + b for a in "IXYZ" for b in "IXYZ"} - {"II"}
    # exp(-i*pi*ZZ/4)|++> has these Gaussian integer numerators over sqrt(8).
    amplitudes = ((1, -1), (1, 1), (1, 1), (1, -1))
    norm = sum(a * a + b * b for a, b in amplitudes)
    assert norm == 8
    reduced = []
    for left in range(2):
        row = []
        for right in range(2):
            numerator = [0, 0]
            for other in range(2):
                a = amplitudes[2 * left + other]
                b = amplitudes[2 * right + other]
                value = _gmul(a, (b[0], -b[1]))
                numerator = [x + y for x, y in zip(numerator, value)]
            row.append((F(numerator[0], norm), F(numerator[1], norm)))
        reduced.append(row)
    assert reduced == [[(F(1, 2), F(0)), (F(0), F(0))],
                       [(F(0), F(0)), (F(1, 2), F(0))]]
    return {"logical_control_generators": sorted(generators),
            "nonidentity_pauli_lie_basis": sorted(closure), "lie_algebra_dimension": 15,
            "ZZ_pi_over_four_maps_plus_plus_to_maximally_entangled_state": True,
            "exact_reduced_density_matrix": [["1/2", "0"], ["0", "1/2"]],
            "noise_threshold_claimed": False}


def check_product_readout():
    group, cosets = base.signed_cosets()
    assert all(base.weight(word) % 2 == 0 for word in group)
    assert all(base.weight(word) % 2 == 1 for axis in "XYZ"
               for words in cosets[axis].values() for word in words)
    normalizer = set(group) | {word for axis in "XYZ" for words in cosets[axis].values() for word in words}
    pair_checks = 0
    for left, right in product(normalizer, repeat=2):
        word, _ = base.multiply(left, right)
        assert base.weight(word) % 2 == ((base.weight(left) + base.weight(right)
                                         + base.anticommutes(left, right)) % 2)
        pair_checks += 1
    for level in range(2, 5):
        n = 5 ** level
        assert all(family.weight(label, n) % 2 == 0 for label in family.code_generators(level))
        assert all(family.weight(family.pauli(axis * n), n) % 2 == 1 for axis in "XYZ")
    return {"seed_stabilizer_even_logical_odd_parity_verified": True,
            "seed_normalizer_parity_product_identities": pair_checks,
            "product_support_criterion": "A nontrivial logical word has an odd-weight restriction on at least one code block; every stabilizer restriction there has even weight, so their exact supports differ.",
            "even_global_ZZ_weight_does_not_give_support_collision": True,
            "syndrome_reversibility_distance": "infinity",
            "known_product_unitary_and_ideal_original_syndrome_assumptions_retained": True}


def check_constructor():
    import logical_pauli_gates as gates
    reports = []
    component_checks = 0
    for level in range(2, 5):
        block_n, d = 5 ** level, 3 ** level
        for word in ("ZZ", "XY", "ZIZ", "IXY", "XYZ", "IYI"):
            active = sum(axis != "I" for axis in word)
            effective = active * d
            expected = F(4 * effective ** 3 + 51 * effective ** 2 + 170 * effective + 24, 24)
            assert gates.target_distance(level, word) == effective
            assert gates.interaction_order(level, word) == effective + 6
            assert gates.budget(level, word) == expected
            assert gates.coefficient_denominator_bound(level, word) == child.coefficient_denominator_bound(level) ** active
            for protection in range(4):
                actual = tuple(gates.increment_components(word, protection))
                expected_indices = {indices for indices in product(range(protection + 1), repeat=active)
                                    if sum(indices) <= protection}
                assert len(actual) == len(set(actual)) == comb(active + protection, protection)
                assert set(actual) == expected_indices
                component_checks += 1
            errors = [(), ((0, "X"),), ((0, "X"), (1, "Y")),
                      ((0, "X"), (1, "Y"), (2, "X")),
                      ((0, "X"), (block_n, "Z")),
                      ((0, "X"), (block_n, "Y"), (block_n + 1, "Z"))]
            if len(word) >= 3:
                errors.append(((block_n, "X"), (block_n + 1, "Y"), (block_n + 2, "Z")))
            for error in errors:
                assert gates.gate_response(level, word, 3, error) == 1
            reports.append({"level": level, "logical_word": word,
                            "target_coset_distance": effective,
                            "interaction_order": effective + 6,
                            "native_budget": str(expected),
                            "actual_recursive_response_checks": len(errors)})
    increments = _base_increments()
    for protection in range(4):
        actual = gates.seed_pair_coefficients(protection)
        expected = {}
        for left, right in product(range(7), repeat=2):
            mass = sum(increments[a][left] * increments[b][right]
                       for a in range(4) for b in range(4) if a + b <= protection)
            if mass:
                expected[base.CLASS_NAMES[left], base.CLASS_NAMES[right]] = mass
        assert actual == expected
    generators = tuple(family.embed(g, 25, 50, block)
                       for block in range(2) for g in family.code_generators(2))
    basis = family.signed_basis(generators, 50)
    memberships = 0
    denominator = child.coefficient_denominator_bound(2) ** 2
    for word in (a + b for a in "XYZ" for b in "XYZ"):
        logical = family.pauli(word[0] * 25 + word[1] * 25)
        contributions = tuple(islice(gates.native_contributions(2, word, 3), 8))
        assert len(contributions) == 8
        for label, coefficient in contributions:
            difference, phase = family.multiply(label, logical, 50)
            phase = (phase + family.stabilizer_phase(difference, basis, 50)) % 4
            assert phase % 2 == 0
            assert coefficient * (-1) ** (phase // 2) > 0
            assert family.weight(label, 50) <= 24
            assert (coefficient * denominator).denominator == 1
            memberships += 1
    return {"recursive_constructor_examples": reports,
            "weak_composition_enumerations_checked": component_checks,
            "collected_seed_pair_tables_checked": 4,
            "native_iterator_prefix_memberships": memberships,
            "prefixes_are_not_full_native_expansions": True}


def main():
    if sys.flags.optimize:
        raise SystemExit("Assertions are scientific checks; do not use -O.")
    result = {"checks_passed": True,
              "fifty_qubit_certificate": check_all_fifty_qubit_errors(),
              "budget_semigroup": check_budget_semigroup(),
              "projected_coset_lower_bound": check_projected_lower_bound(),
              "finite_target_configurations": check_target_configurations(),
              "universal_control": check_universal_control(),
              "product_readout": check_product_readout(),
              "constructor": check_constructor()}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
