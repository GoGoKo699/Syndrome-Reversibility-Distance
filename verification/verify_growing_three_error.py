#!/usr/bin/env python3
"""Exact certificates for recursive three-error gates at growing distance.

Coefficientwise polynomial checks certify the moment masses, dual and cost
recurrences. Independent finite class profiles give complete coverage of all
8,649,376 protected errors on 125 qubits without expanding the large gate.
Only the standard library is required. Assertions must remain enabled.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, islice, permutations, product
from math import comb
import json
import sys

import verify_three_error_gates as base


def poly(*values):
    values = tuple(map(F, values))
    while len(values) > 1 and not values[-1]:
        values = values[:-1]
    return values


def add(*values):
    return poly(*(sum((row[i] if i < len(row) else 0) for row in values)
                  for i in range(max(map(len, values)))))


def scale(values, factor):
    return poly(*(factor * value for value in values))


def mul(left, right):
    out = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return poly(*out)


def power(values, exponent):
    answer = poly(1)
    for _ in range(exponent):
        answer = mul(answer, values)
    return answer


def evaluate(values, value):
    return sum(coefficient * value ** i for i, coefficient in enumerate(values))


def symbolic_certificates():
    """Check every coefficient, rather than fitting values at finite depths."""
    d = poly(0, 1)
    columns = []
    for offset, square, cube in ((0, 0, 0), (2, 16, 98), (4, 32, 196),
                                  (6, 48, 294), (4, 16, 64), (6, 32, 162),
                                  (6, 18, 54)):
        w, r, s = poly(offset, 1), poly(square, 3), poly(cube, 9)
        columns.append((poly(1), w, power(w, 2), power(w, 3), r, mul(w, r), s))
    masses = (
        poly(1, F(11, 12), F(1, 4), F(1, 48)),
        poly(0, -F(3, 2), -F(5, 8), -F(1, 16)),
        poly(0, F(27, 16), F(21, 32), F(1, 16)),
        poly(0, -F(15, 16), -F(9, 32), -F(1, 48)),
        poly(0, -F(15, 16), -F(5, 32)),
        poly(0, F(15, 16), F(5, 32)),
        poly(0, -F(1, 6)),
    )
    polarity = (1, -1, 1, -1, -1, 1, -1)
    budget = poly(1, F(85, 12), F(17, 8), F(1, 6))
    dual = (budget, poly(-F(101, 6), -F(37, 8), -F(1, 2)),
            poly(F(5, 2), F(1, 2)), poly(-F(1, 6)),
            poly(F(19, 4), F(1, 8)), poly(-F(1, 8)), poly(-F(1, 2)))
    for row in range(7):
        assert add(*(mul(mass, column[row]) for mass, column in zip(masses, columns))) == poly(int(row == 0))
    for sign, column in zip(polarity, columns):
        assert add(*(mul(y, entry) for y, entry in zip(dual, column))) == poly(sign)
    assert add(*(scale(mass, sign) for mass, sign in zip(masses, polarity))) == budget
    increments = (poly(1), d, poly(0, F(13, 8), F(1, 2)),
                  poly(0, F(107, 24), F(13, 8), F(1, 6)))
    assert add(*increments) == budget
    for order in range(4):
        cube_coefficient = add(*(mul(mul(increments[a], increments[b]), increments[c])
                                 for a, b, c in product(range(4), repeat=3)
                                 if a + b + c == order))
        substituted = poly(*(coefficient * 3 ** i for i, coefficient in enumerate(increments[order])))
        assert cube_coefficient == substituted
    records = []
    for level in range(2, 7):
        distance = 3 ** level
        actual_masses = tuple(evaluate(row, distance) for row in masses)
        assert all(value * sign > 0 for value, sign in zip(actual_masses, polarity))
        actual_budget = evaluate(budget, distance)
        assert sum(map(abs, actual_masses)) == actual_budget
        records.append({"level": level, "n": 5 ** level, "distance": distance,
                        "interaction_order": distance + 6,
                        "class_masses": list(map(str, actual_masses)),
                        "J_per_logical_half_splitting": str(actual_budget)})
    return {"coefficientwise_mass_and_dual_identities": True,
            "coefficientwise_increment_cost_cube_recurrence": True,
            "mass_polynomials_in_d": [list(map(str, row)) for row in masses],
            "budget_polynomial_in_d": list(map(str, budget)),
            "finite_hierarchy_values": records}


PROFILE_NAMES = ("identity", "single_equivalent", "double_in_one_block_equivalent",
                 "two_distinct_blocks", "normalizer_triple_in_one_block",
                 "two_plus_one_blocks", "three_distinct_blocks")


def checked_child_profiles():
    """Prove profiles are axis independent and retain exact error-weight counts."""
    _, cosets = base.signed_cosets()
    for axis in "XYZ":
        supports = {tuple(i for i, value in enumerate(word) if value != "I")
                    for word in cosets[axis][3]}
        assert supports == set(combinations(range(5), 3))
    table = base.character_tables(cosets)
    records, _ = base.physical_errors(table)
    profile_index = {base.PROFILES[name]: i for i, name in enumerate(PROFILE_NAMES)}
    previous = None
    by_weight = Counter()
    for axis in "XYZ":
        configurations = base.configurations(axis, cosets)
        encoded = [[tuple(base.OPTION_INDEX[item] for item in configuration)
                    for configuration in group] for group in configurations]
        cache = {}
        local_counts = Counter()
        for signature, pattern, _ in records:
            if signature not in cache:
                cache[signature] = base.factorized_profile(signature, encoded)
            profile = cache[signature]
            index = profile_index[profile]
            local_counts[index, sum(pattern)] += 1
        if previous is not None:
            assert cache == previous, "A child profile must not depend on its prescribed logical axis."
        previous = cache
        if by_weight:
            assert by_weight == local_counts
        by_weight = local_counts
    assert sum(by_weight.values()) == 64876
    return by_weight


def smolyak_components():
    zero = (F(0),) * 7
    gates = ((F(1),) + zero[1:],
             tuple(map(F, ("11/2", "-9/2", "0", "0", "0", "0", "0"))),
             tuple(map(F, ("143/8", "-117/4", "243/16", "0", "-45/16", "0", "0"))),
             base.MASSES)
    increments = (gates[0],) + tuple(tuple(a - b for a, b in zip(gates[order], gates[order - 1]))
                                    for order in range(1, 4))
    answer = defaultdict(F)
    for orders in product(range(4), repeat=3):
        if sum(orders) > 3:
            continue
        for labels in product(range(7), repeat=3):
            coefficient = F(1)
            for order, label in zip(orders, labels):
                coefficient *= increments[order][label]
            if coefficient:
                answer[tuple(sorted(labels))] += coefficient
    return {labels: coefficient for labels, coefficient in answer.items() if coefficient}


def component_response(labels, errors):
    """Independent outer-support tensor average, using the checked child rows."""
    total, count = F(0), 0
    for active in combinations(range(5), 3):
        for order in set(permutations(labels)):
            value = F(1)
            for block, profile in enumerate(errors):
                if block in active:
                    value *= base.PROFILES[PROFILE_NAMES[profile]][order[active.index(block)]]
            total += value
            count += 1
    return total / count


def check_depth_three():
    child_counts = checked_child_profiles()
    states = {((), 0): 1}
    for _ in range(5):
        updated = defaultdict(int)
        for (profiles, old_weight), count in states.items():
            for (profile, weight), multiplicity in child_counts.items():
                if old_weight + weight <= 3:
                    new_profiles = profiles if profile == 0 else tuple(sorted(profiles + (profile,)))
                    updated[new_profiles, old_weight + weight] += count * multiplicity
        states = updated
    coverage = Counter()
    for (profiles, _), count in states.items():
        coverage[profiles] += count
    assert len(coverage) == 11 and sum(coverage.values()) == 8649376
    assert sum(coverage.values()) == sum(3 ** weight * comb(125, weight) for weight in range(4))
    coefficients = smolyak_components()
    expected = {(0, 0, 0): F(9889, 16), (0, 0, 1): -F(27621, 16),
                (0, 0, 2): F(24057, 32), (0, 1, 1): F(8019, 8),
                (0, 0, 3): -F(4455, 32), (0, 1, 2): -F(6561, 16),
                (1, 1, 1): -F(729, 8), (0, 0, 4): -F(4455, 32),
                (0, 0, 5): F(2025, 32), (0, 1, 4): F(1215, 16),
                (0, 0, 6): -F(9, 2)}
    assert coefficients == expected
    for errors in coverage:
        assert sum(coefficient * component_response(labels, errors)
                   for labels, coefficient in coefficients.items()) == 1
    assert sum(map(abs, coefficients.values())) == F(40175, 8)
    return {"logical_axis_independent_child_profiles": True,
            "exact_parent_profile_equations": len(coverage),
            "protected_physical_errors_covered": sum(coverage.values()),
            "profile_coverage": [{"nonidentity_child_profiles": [PROFILE_NAMES[i] for i in errors],
                                   "physical_error_labels": count}
                                  for errors, count in sorted(coverage.items())],
            "signed_component_masses": [{"child_classes": [base.CLASS_NAMES[i] for i in labels],
                                         "mass": str(coefficient)}
                                        for labels, coefficient in sorted(coefficients.items())],
            "interaction_order": 33, "J_per_logical_half_splitting": "40175/8"}


def check_truncation():
    counts = []
    for protection in range(4):
        partitions = 0
        omitted = 0
        for error_weights in product(range(4), repeat=3):
            if sum(error_weights) > protection:
                continue
            partitions += 1
            for orders in product(range(4), repeat=3):
                if sum(orders) > protection:
                    assert any(order > weight for order, weight in zip(orders, error_weights))
                    omitted += 1
        counts.append({"protected_error_weight": protection,
                       "active_child_error_count_partitions": partitions,
                       "omitted_tensor_terms_with_zero_factor": omitted})
    return {"truncated_tensor_zero_property_verified": True, "orders": counts}


def check_pattern_polarities():
    ranks = (0, 1, 2, 3, 2, 3, 3)
    signs = (1, -1, 1, -1, -1, 1, -1)
    # Added active bottom blocks, number of weight-five factors, number of
    # weight-four stabilizer factors, relative to the minimum logical tree.
    patterns = ((0, 0, 0), (0, 1, 0), (0, 2, 0), (0, 3, 0),
                (0, 0, 1), (0, 1, 1), (2, 0, 0))
    index = {pattern: i for i, pattern in enumerate(patterns)}
    checks = 0
    for labels in product(range(7), repeat=3):
        if sum(ranks[label] for label in labels) > 3:
            continue
        combined = tuple(sum(patterns[label][coordinate] for label in labels)
                         for coordinate in range(3))
        expected = signs[index[combined]]
        actual = signs[labels[0]] * signs[labels[1]] * signs[labels[2]]
        assert actual == expected
        checks += 1
    return {"recursive_native_pattern_sign_coherence": True,
            "ordered_child_class_compositions_checked": checks}


def check_reference_responses():
    import three_error_recursion as recursion
    reports = []
    for level in range(2, 7):
        n, child_n = 5 ** level, 5 ** (level - 1)
        errors = [(), ((0, "X"),), ((0, "X"), (1, "X")),
                  ((0, "X"), (5, "X")),
                  ((0, "X"), (1, "Y"), (2, "X")),
                  ((0, "X"), (1, "X"), (5, "X")),
                  ((0, "X"), (5, "X"), (10, "X")),
                  ((0, "X"), (child_n, "X")),
                  ((0, "X"), (child_n, "Y"), (2 * child_n, "Z")),
                  ((0, "X"), (1, "X"), (child_n, "X")),
                  ((0, "X"), (5, "X"), (child_n, "X"))]
        errors = list(dict.fromkeys(error for error in errors
                                   if len({site for site, _ in error}) == len(error)))
        d = 3 ** level
        expected_masses = (F((d + 2) * (d + 4) * (d + 6), 48),
                           -F(d * (d + 4) * (d + 6), 16),
                           F(d * (2 * d + 9) * (d + 6), 32),
                           -F(d * (2 * d + 15) * (d + 6), 96),
                           -F(5 * d * (d + 6), 32), F(5 * d * (d + 6), 32), -F(d, 6))
        assert tuple(recursion.class_masses(level).values()) == expected_masses
        expected_increments = (F(1), F(d), F(d * (4 * d + 13), 8),
                               F(d * (4 * d * d + 39 * d + 107), 24))
        assert recursion.increment_budgets(level) == expected_increments
        assert recursion.budget(level) == sum(map(abs, expected_masses)) == sum(expected_increments)
        for coefficients in recursion.increment_coefficients(level):
            assert all(coefficient * sign >= 0 for coefficient, sign in zip(
                coefficients, (1, -1, 1, -1, -1, 1, -1)))
        for column, sign in zip(recursion.moment_columns(level), (1, -1, 1, -1, -1, 1, -1)):
            assert sum(y * x for y, x in zip(recursion.dual(level), column)) == sign
        equations = 0
        for axis in "XYZ":
            for error in errors:
                assert all(site < n for site, _ in error)
                for protection in range(len(error), 4):
                    assert recursion.gate_response(level, protection, axis, error) == 1
                    equations += 1
                for order in range(len(error) + 1, 4):
                    assert recursion.increment_response(level, order, axis, error) == 0
                    equations += 1
        reports.append({"level": level, "n": n, "axes": ["X", "Y", "Z"],
                        "physical_error_witnesses_per_axis": len(errors),
                        "exact_gate_and_zero_increment_equations": equations})
    return {"recursive_reference_responses_verified": True, "levels": reports}


def check_native_specification():
    import three_error_recursion as recursion
    import verify_readout_safe_family as family
    denominator = recursion.coefficient_denominator_bound(2)
    assert denominator == 3200000
    class_coefficients = 0
    for protection in range(4):
        for mass, count in zip(recursion.class_masses(2, protection).values(), base.COUNTS):
            if mass:
                assert (mass * denominator / count).denominator == 1
                class_coefficients += 1
    for coefficients in recursion.increment_coefficients(2):
        for mass, count in zip(coefficients, base.COUNTS):
            if mass:
                assert (mass * denominator / count).denominator == 1
                class_coefficients += 1
    hierarchy = []
    for level in range(2, 7):
        power_three = 3 ** (level - 2)
        expected = 3200000 ** power_three * 10 ** ((power_three - 1) // 2)
        assert recursion.coefficient_denominator_bound(level) == expected
        if level > 2:
            assert expected == 10 * recursion.coefficient_denominator_bound(level - 1) ** 3
        hierarchy.append({"level": level, "safe_common_denominator_bits": expected.bit_length()})
    basis = family.signed_basis(family.code_generators(3), 125)
    memberships = 0
    for axis in "XYZ":
        logical = family.pauli(axis * 125)
        for order in range(4):
            words = tuple(islice(recursion.increment_native_contributions(3, order, axis), 20))
            assert len(words) == 20
            for label, coefficient in words:
                difference, phase = family.multiply(label, logical, 125)
                phase = (phase + family.stabilizer_phase(difference, basis, 125)) % 4
                assert phase % 2 == 0
                relative = coefficient * (-1) ** (phase // 2)
                assert relative > 0  # These deterministic prefixes belong to A0.
                assert family.weight(label, 125) <= 27 + 2 * order
                assert (coefficient * recursion.coefficient_denominator_bound(3)).denominator == 1
                memberships += 1
    return {"base_gate_and_increment_uniform_class_coefficients_checked": class_coefficients,
            "denominator_bounds": hierarchy,
            "deterministic_iterator_prefix_memberships_checked": memberships,
            "iterator_prefixes_are_not_full_native_expansions": True}


def main():
    if sys.flags.optimize:
        raise SystemExit("Assertions are scientific checks; do not use -O.")
    result = {"checks_passed": True, "symbolic_certificates": symbolic_certificates(),
              "depth_three_profile_certificate": check_depth_three(),
              "truncation_certificate": check_truncation(),
              "pattern_polarity_certificate": check_pattern_polarities(),
              "reference_response_checks": check_reference_responses(),
              "native_specification_checks": check_native_specification()}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
