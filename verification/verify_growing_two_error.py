#!/usr/bin/env python3
"""Exact certificates for growing-distance two-error-transparent gates.

The all-depth response cancellation is checked as polynomial identities, not
by fitting numerical data. Independent finite checks include the literal
five-qubit Pauli algebra, all 10,000 depth-three active-support patterns,
frontier-padding probabilities, and classification of all 70,126 protected
Pauli errors on 125 qubits. No large Hilbert-space matrix or dense logical
coset is constructed. Only Python's standard library is required.
"""

from fractions import Fraction as F
from itertools import combinations, product
from math import comb
import json
import sys

import verify_two_error_gates as base


# Polynomials are tuples of rational coefficients in increasing powers of a.
def poly(*values):
    out = tuple(map(F, values))
    while len(out) > 1 and out[-1] == 0:
        out = out[:-1]
    return out


def add(*polynomials):
    out = [F(0)] * max(map(len, polynomials))
    for values in polynomials:
        for index, value in enumerate(values):
            out[index] += value
    return poly(*out)


def scale(values, factor):
    return poly(*(factor * value for value in values))


def mul(left, right):
    out = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return poly(*out)


def dot(left, right):
    return add(*(mul(a, b) for a, b in zip(left, right)))


def assert_zero(values):
    assert values == (F(0),), values


def symbolic_identities():
    """Prove the rational response identities for every integer a >= 3."""
    one, a = poly(1), poly(0, 1)
    d = scale(a, 3)
    masses = (
        scale(mul(add(d, poly(2)), add(d, poly(4))), F(1, 8)),
        scale(mul(d, add(d, poly(4))), -F(1, 4)),
        scale(mul(d, add(scale(d, 2), poly(9))), F(1, 16)),
        scale(d, -F(5, 16)),
    )
    weights = (d, add(d, poly(2)), add(d, poly(4)), add(d, poly(4)))
    block_squares = (scale(d, 3), add(scale(d, 3), poly(16)),
                     add(scale(d, 3), poly(32)), add(scale(d, 3), poly(16)))
    assert_zero(add(dot(masses, (one,) * 4), poly(-1)))
    assert_zero(dot(masses, weights))
    assert_zero(dot(masses, tuple(mul(w, w) for w in weights)))
    assert_zero(dot(masses, block_squares))
    budget = add(masses[0], scale(masses[1], -1), masses[2], scale(masses[3], -1))
    expected_budget = scale(add(scale(mul(d, d), 4), scale(d, 21), poly(8)), F(1, 8))
    assert_zero(add(budget, scale(expected_budget, -1)))

    low, high, padding = F(1, 5), -F(1, 3), -F(1, 15)
    low_pair, high_pair = -F(1, 15), F(1, 9)
    delta = high - low
    # Single-error responses are 1 + these expressions / N.
    single_numerators = tuple(add(scale(a, low - 1), poly(k * delta))
                              for k in range(3)) + (
                                  add(scale(a, low - 1), poly(padding - 1)),)
    assert_zero(dot(masses, single_numerators))
    same_numerators = tuple(add(scale(a, low_pair - 1), poly(k * (high_pair - low_pair)))
                            for k in range(3)) + (
                                add(scale(a, low_pair - 1), poly(padding - 1)),)
    assert_zero(dot(masses, same_numerators))

    # For distinct bottom blocks, the remaining coefficient multiplies rho,
    # their joint active probability. Clear the nonzero denominator a(a-1).
    denominator = mul(a, add(a, poly(-1)))
    pair_numerators = tuple(add(
        scale(denominator, (low - 1) ** 2),
        scale(add(a, poly(-1)), 2 * k * (low - 1) * delta),
        poly(k * (k - 1) * delta ** 2)) for k in range(3)) + (
            add(scale(denominator, (low - 1) ** 2),
                scale(a, 2 * (low - 1) * (padding - 1))),)
    assert_zero(dot(masses, pair_numerators))
    return {
        "polynomial_variable": "a = 3**(m-1), number of active bottom blocks",
        "domain": "m >= 2, hence a >= 3",
        "mass_polynomials_in_a": [[str(x) for x in values] for values in masses],
        "normalization_identity": True,
        "global_first_second_moments": True,
        "bottom_block_squared_weight_moment": True,
        "single_error_inverse_N_coefficient_zero": True,
        "same_bottom_pair_inverse_N_coefficient_zero": True,
        "distinct_bottom_pair_rho_coefficient_zero": True,
        "coupling_budget_polynomial_in_a": [str(x) for x in budget],
    }


def check_base_pauli_algebra():
    group = base.signed_group()
    cosets, normalizer = base.base_cosets(group)
    single, double = base.inner_errors()
    tables = base.inner_response_tables(cosets, group, normalizer, single, double)
    for error in double:
        assert sum(tables[axis, 3, error] for axis in "XYZ") / 3 == -F(1, 15)
        assert sum(tables[axis, 5, error] for axis in "XYZ") / 3 == F(1, 9)
    for error in single + double:
        assert tables["I", 4, error] == -F(1, 15)
    assert all(tables[axis, w, "IIIII"] == 1
               for axis, classes in cosets.items() for w in classes)
    return {
        "physical_five_qubit_labels_enumerated": 1024,
        "signed_stabilizers": len(group),
        "normalizer_labels": len(normalizer),
        "identity_single_double_errors": 106,
        "logical_axes": ["X", "Y", "Z"],
        "all_inner_responses_checked": True,
    }


def joint_probabilities(height, divergence):
    """Independent recursion for active-active and active-frontier-pad pairs."""
    assert 1 <= divergence <= height
    a, inner_a = 3 ** height, 3 ** (height - 1)
    inner_N = 5 ** (height - 1)
    inner_p = F(inner_a, inner_N)
    if divergence == 1:
        rho = F(3, 10) * inner_p ** 2
        # The target root child is inactive (one frontier) or active
        # (inner_a-1 frontiers); its partner's active event is independent.
        eta = F(3, 10) * inner_p * inner_a / ((a - 1) * inner_N)
    else:
        old_rho, old_eta = joint_probabilities(height - 1, divergence - 1)
        rho = F(3, 5) * old_rho
        eta = F(3, 5) * F(inner_a - 1, a - 1) * old_eta
    return rho, eta


def responses(a, N, rho=None):
    """Rows for I, one error, same-bottom pair, or a distinct-bottom pair."""
    p = F(a, N)
    low, high, padding = F(1, 5), -F(1, 3), -F(1, 15)
    low_pair, high_pair = -F(1, 15), F(1, 9)
    single, same, pair = [], [], []
    for k in range(3):
        r = low + F(k, a) * (high - low)
        r_pair = low_pair + F(k, a) * (high_pair - low_pair)
        single.append(1 - p + p * r)
        same.append(1 - p + p * r_pair)
        if rho is not None:
            s = low ** 2 + 2 * F(k, a) * low * (high - low) + (
                F(k * (k - 1), a * (a - 1)) * (high - low) ** 2)
            pair.append(1 + 2 * p * (r - 1) + rho * (1 - 2 * r + s))
    single.append(1 + p * (low - 1) + F(padding - 1, N))
    same.append(1 + p * (low_pair - 1) + F(padding - 1, N))
    if rho is not None:
        pair.append(1 + 2 * p * (low - 1) + 2 * F(padding - 1, N)
                    + rho * ((low - 1) ** 2 + 2 * (low - 1) * (padding - 1) / (a - 1)))
    return (F(1),) * 4, tuple(single), tuple(same), tuple(pair)


def finite_hierarchy_checks():
    records = []
    for depth in range(2, 7):
        height = depth - 1
        n, d, a, N = 5 ** depth, 3 ** depth, 3 ** height, 5 ** height
        matrix = ((1, 1, 1, 1), (d, d + 2, d + 4, d + 4),
                  (d * d, (d + 2) ** 2, (d + 4) ** 2, (d + 4) ** 2),
                  (3 * d, 3 * d + 16, 3 * d + 32, 3 * d + 16))
        masses = base.solve_square(matrix, (1, 0, 0, 0))
        expected = (F((d + 2) * (d + 4), 8), -F(d * (d + 4), 4),
                    F(d * (2 * d + 9), 16), -F(5 * d, 16))
        assert masses == expected
        signs = (1, -1, 1, -1)
        assert tuple(1 if mass > 0 else -1 for mass in masses) == signs
        dual = base.solve_square(tuple(zip(*matrix)), signs)
        budget = sum(map(abs, masses))
        assert dual[0] == budget == F(4 * d * d + 21 * d + 8, 8)
        assert all(sum(dual[i] * matrix[i][j] for i in range(4)) == signs[j]
                   for j in range(4))
        hierarchy = []
        identity, single, same, _ = responses(a, N)
        for row in (identity, single, same):
            assert sum(mass * response for mass, response in zip(masses, row)) == 1
        for divergence in range(1, height + 1):
            rho, eta = joint_probabilities(height, divergence)
            assert rho == F(3, 5) ** (divergence - 1) * F(3, 10) * F(3, 5) ** (2 * (height - divergence))
            assert eta == rho / (a - 1)
            row = responses(a, N, rho)[3]
            assert sum(mass * response for mass, response in zip(masses, row)) == 1
            pairs = 10 * 5 ** (2 * height - divergence - 1)
            hierarchy.append({
                "first_divergence_depth": divergence,
                "bottom_block_pairs": pairs,
                "physical_two_error_labels": 225 * pairs,
                "joint_active_probability": str(rho),
                "active_and_frontier_padding_probability": str(eta),
                "class_responses": [str(value) for value in row],
            })
        assert sum(item["bottom_block_pairs"] for item in hierarchy) == comb(N, 2)
        assert 1 + 3 * n + 90 * N + sum(item["physical_two_error_labels"] for item in hierarchy) == 1 + 3 * n + 9 * comb(n, 2)
        # Uniform inactive frontier choice has N^{-1} marginal at every leaf.
        child_a, child_N = 3 ** (height - 1), 5 ** (height - 1)
        padding_marginal = (F(2, 5) + F(3, 5) * (child_a - 1)) / ((a - 1) * child_N)
        assert padding_marginal == F(1, N)
        frontier_counts = [2 * 3 ** (r - 1) for r in range(1, height + 1)]
        assert sum(frontier_counts) == a - 1
        records.append({
            "depth": depth, "physical_qubits": n, "ordinary_distance_from_written_proof": d,
            "protected_error_weight": 2, "interaction_order": d + 4,
            "class_masses": [str(value) for value in masses],
            "block_moment_dual": [str(value) for value in dual],
            "coupling_budget_per_logical_half_splitting": str(budget),
            "frontier_count": a - 1, "frontiers_by_depth": frontier_counts,
            "single_error_class_responses": [str(value) for value in single],
            "same_bottom_pair_class_responses": [str(value) for value in same],
            "total_identity_single_double_errors": 1 + 3 * n + 9 * comb(n, 2),
            "distinct_bottom_pair_hierarchy": hierarchy,
        })
    # At depth two the frontier is just the two inactive bottom blocks.
    first = records[0]
    assert first["single_error_class_responses"] == [str(x) for x in base.EXPECTED_RESPONSES["single"]]
    assert first["same_bottom_pair_class_responses"] == [str(x) for x in base.EXPECTED_RESPONSES["two_in_one_block"]]
    assert first["distinct_bottom_pair_hierarchy"][0]["class_responses"] == [str(x) for x in base.EXPECTED_RESPONSES["two_in_distinct_blocks"]]
    return records


def brute_depth_three_supports():
    """Enumerate support patterns independently of the probability recurrence."""
    pairs = ((0, 1), (0, 5))
    q_counts = [0, 0]
    eta_units = [0, 0]
    active_count = padding_units = count = 0
    for parents in combinations(range(5), 3):
        for children in product(tuple(combinations(range(5), 3)), repeat=3):
            active = {5 * parent + child for parent, selected in zip(parents, children)
                      for child in selected}
            # There are 8 frontiers: 2 root holes, each holding 5 bottom blocks,
            # and 6 bottom holes. Use the common probability denominator 40.
            units = {}
            for parent in range(5):
                if parent not in parents:
                    units.update({5 * parent + child: 1 for child in range(5)})
            for parent, selected in zip(parents, children):
                units.update({5 * parent + child: 5 for child in range(5) if child not in selected})
            assert len(active) == 9 and len(units) == 16
            assert set(units).isdisjoint(active) and sum(units.values()) == 40
            count += 1
            active_count += 0 in active
            padding_units += units.get(0, 0)
            for index, (left, right) in enumerate(pairs):
                q_counts[index] += left in active and right in active
                eta_units[index] += int(left in active) * units.get(right, 0)
    assert count == 10000
    assert F(active_count, count) == F(9, 25)
    assert F(padding_units, 40 * count) == F(1, 25)
    expected_divergence = (2, 1)
    actual = []
    for index, divergence in enumerate(expected_divergence):
        rho, eta = joint_probabilities(2, divergence)
        assert F(q_counts[index], count) == rho
        assert F(eta_units[index], 40 * count) == eta
        # Negative control: the superficially similar uniform-inactive-leaf
        # rule fails, even though its one-point padding marginal is correct.
        naive_eta = (F(9, 25) - rho) / 16
        naive_row = list(responses(9, 25, rho)[3])
        naive_row[-1] += F(128, 75) * (naive_eta - eta)
        masses = (F(899, 8), -F(837, 4), F(1701, 16), -F(135, 16))
        naive_rate = sum(mass * value for mass, value in zip(masses, naive_row))
        assert naive_rate == {1: F(2419, 2500), 2: F(581, 500)}[divergence]
        assert naive_rate != 1
        actual.append({"bottom_block_pair": list(pairs[index]),
                       "joint_active_probability": str(rho),
                       "active_and_padding_probability": str(eta),
                       "uniform_inactive_padding_rate_negative_control": str(naive_rate)})
    return {"minimum_support_patterns_enumerated": count,
            "frontier_probabilities_checked": True, "pair_cases": actual}


def enumerate_depth_three_errors():
    """Enumerate all125-qubit error labels by their distinct sites and axes."""
    n, height = 125, 2
    counts = {"identity": 1, "single": 0, "same_bottom_pair": 0,
              "distinct_bottom_divergence_1": 0, "distinct_bottom_divergence_2": 0}
    labels = {()}
    for site in range(n):
        for axis in "XYZ":
            label = ((site, axis),)
            assert label not in labels
            labels.add(label)
            counts["single"] += 1
    for left, right in combinations(range(n), 2):
        a, b = left // 5, right // 5
        if a == b:
            key = "same_bottom_pair"
        else:
            digits_a = tuple((a // 5 ** j) % 5 for j in reversed(range(height)))
            digits_b = tuple((b // 5 ** j) % 5 for j in reversed(range(height)))
            divergence = next(j + 1 for j, (x, y) in enumerate(zip(digits_a, digits_b)) if x != y)
            key = f"distinct_bottom_divergence_{divergence}"
        for axes in product("XYZ", repeat=2):
            label = ((left, axes[0]), (right, axes[1]))
            assert label not in labels
            labels.add(label)
            counts[key] += 1
    assert counts == {"identity": 1, "single": 375, "same_bottom_pair": 2250,
                      "distinct_bottom_divergence_1": 56250,
                      "distinct_bottom_divergence_2": 11250}
    assert len(labels) == sum(counts.values()) == 70126
    return {"physical_qubits": n, "identity_single_double_error_labels_enumerated": len(labels),
            "error_category_counts": counts, "classification_complete": True}


def main():
    if sys.version_info < (3, 10):
        raise SystemExit("Python 3.10 or later is required.")
    if sys.flags.optimize:
        raise SystemExit("Assertions are scientific checks; do not use -O.")
    print(json.dumps({"checks_passed": True,
                      "symbolic_all_depth_response_identities": symbolic_identities(),
                      "base_signed_pauli_checks": check_base_pauli_algebra(),
                      "finite_hierarchy_instances": finite_hierarchy_checks(),
                      "independent_depth_three_support_enumeration": brute_depth_three_supports(),
                      "independent_depth_three_error_enumeration": enumerate_depth_three_errors()},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
