#!/usr/bin/env python3
"""Exact two-error-transparent gates on the concatenated five-qubit code.

Only the standard library is required. The program independently constructs
the signed five-qubit Pauli group using literal local multiplication. It then
checks all 2,776 physical Pauli errors of weight at most two, for each logical
axis, by exact tensor-factor response calculations. It does not enumerate a
2**25-dimensional Hilbert space or all 338,800 native terms per axis.

The full-space conclusion uses signed logical-coset membership of every
factor, not compressed energies alone. Universal lower bounds still depend
on the written proof; the finite matrix and its dual certificate are checked
here exactly. Do not disable assertions with Python's -O option.
"""

from fractions import Fraction
from itertools import combinations, product
import json
from math import prod
import sys

from two_error_sparse import check_sparse_construction


GENERATORS = ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
CLASS_NAMES = ("A0", "A1", "A2", "B")
MASSES = (Fraction(143, 8), -Fraction(117, 4),
          Fraction(243, 16), -Fraction(45, 16))
MOMENT_MATRIX = (
    (1, 1, 1, 1),
    (9, 11, 13, 13),
    (81, 121, 169, 169),
    (27, 43, 59, 43),
)
EXPECTED_RESPONSES = {
    "identity": (Fraction(1),) * 4,
    "single": (Fraction(13, 25), Fraction(31, 75),
               Fraction(23, 75), Fraction(23, 75)),
    "two_in_one_block": (Fraction(9, 25), Fraction(89, 225),
                         Fraction(97, 225), Fraction(11, 75)),
    "two_in_distinct_blocks": (Fraction(29, 125), Fraction(13, 125),
                               Fraction(1, 225), Fraction(23, 375)),
}


def multiply_words(left, right):
    """Return (bare word, power) with P(left)P(right)=i**power P(word)."""
    assert len(left) == len(right)
    forward = {("X", "Y"): "Z", ("Y", "Z"): "X", ("Z", "X"): "Y"}
    answer, phase = [], 0
    for a, b in zip(left, right):
        if a == "I":
            answer.append(b)
        elif b == "I":
            answer.append(a)
        elif a == b:
            answer.append("I")
        elif (a, b) in forward:
            answer.append(forward[a, b])
            phase += 1
        else:
            answer.append(forward[b, a])
            phase -= 1
    return "".join(answer), phase % 4


def anticommutes(left, right):
    return sum(a != "I" and b != "I" and a != b
               for a, b in zip(left, right)) % 2


def word_weight(word):
    return len(word) - word.count("I")


def signed_group():
    group = {"IIIII": 1}
    for generator in GENERATORS:
        added = {}
        for word, sign in group.items():
            new_word, phase = multiply_words(word, generator)
            assert phase % 2 == 0
            assert new_word not in group and new_word not in added
            added[new_word] = sign * (-1) ** (phase // 2)
        group.update(added)
    assert len(group) == 16
    assert group["IIIII"] == 1
    assert all(word_weight(word) == 4 for word in group if word != "IIIII")
    for a, sign_a in group.items():
        for b, sign_b in group.items():
            assert not anticommutes(a, b)
            word, phase = multiply_words(a, b)
            assert phase % 2 == 0
            assert group[word] == sign_a * sign_b * (-1) ** (phase // 2)
    return group


def base_cosets(group):
    cosets = {"I": {0: {"IIIII": 1},
                     4: {word: sign for word, sign in group.items() if word != "IIIII"}}}
    labels_seen = set(group)
    for axis in "XYZ":
        classes = {3: {}, 5: {}}
        logical = axis * 5
        assert all(not anticommutes(logical, g) for g in GENERATORS)
        for stabilizer, sign in group.items():
            word, phase = multiply_words(logical, stabilizer)
            assert phase % 2 == 0
            assert word not in labels_seen
            labels_seen.add(word)
            classes[word_weight(word)][word] = sign * (-1) ** (phase // 2)
        assert len(classes[3]) == 10 and len(classes[5]) == 6
        # Every outer support occurs once; each axis occurs twice at each site.
        supports = {tuple(i for i, a in enumerate(word) if a != "I")
                    for word in classes[3]}
        assert supports == set(combinations(range(5), 3))
        assert all(sum(word[site] == a for word in classes[3]) == 2
                   for site in range(5) for a in "XYZ")
        cosets[axis] = classes
    # Exhaust all 4**5 labels to check that the four signed cosets classify the
    # complete normalizer, not only the supplied representatives.
    normalizer = {"".join(letters) for letters in product("IXYZ", repeat=5)
                  if all(not anticommutes("".join(letters), g) for g in GENERATORS)}
    assert normalizer == labels_seen and len(normalizer) == 64
    assert min(word_weight(word) for word in normalizer - set(group)) == 3
    return cosets, normalizer


def inner_errors():
    single = []
    for site in range(5):
        for axis in "XYZ":
            single.append("I" * site + axis + "I" * (4 - site))
    double = []
    for sites in combinations(range(5), 2):
        for axes in product("XYZ", repeat=2):
            word = ["I"] * 5
            for site, axis in zip(sites, axes):
                word[site] = axis
            double.append("".join(word))
    assert len(single) == 15 and len(double) == 90
    return single, double


def syndrome(word):
    return tuple(anticommutes(word, generator) for generator in GENERATORS)


def inner_response_tables(cosets, group, normalizer, single, double):
    correction = {syndrome(error): error for error in single}
    assert len(correction) == 15 and (0,) * 4 not in correction
    logical_labels = normalizer - set(group)
    for error in double:
        assert syndrome(error) in correction
        returned, _ = multiply_words(correction[syndrome(error)], error)
        assert returned in logical_labels and word_weight(returned) == 3
    tables = {}
    for axis, classes in cosets.items():
        for w, terms in classes.items():
            for error in ("IIIII", *single, *double):
                tables[axis, w, error] = Fraction(
                    sum((-1) ** anticommutes(word, error) for word in terms), len(terms))
    assert all(tables["I", 4, error] == -Fraction(1, 15)
               for error in single + double)
    for axis in "XYZ":
        assert all(tables[axis, 3, error] == Fraction(1, 5) and
                   tables[axis, 5, error] == -Fraction(1, 3) for error in single)
    return tables


def class_configurations(axis, cosets):
    configurations = [[] for _ in CLASS_NAMES]
    for outer, outer_sign in cosets[axis][3].items():
        assert outer_sign in (-1, 1)
        active = [i for i, a in enumerate(outer) if a != "I"]
        inactive = [i for i, a in enumerate(outer) if a == "I"]
        for index, promotions in ((0, [()]),
                                  (1, list(combinations(active, 1))),
                                  (2, list(combinations(active, 2)))):
            for promoted in promotions:
                configurations[index].append(tuple(
                    (outer[i], 0 if i in inactive else 5 if i in promoted else 3)
                    for i in range(5)))
        for padded in inactive:
            configurations[3].append(tuple(
                (outer[i], 4 if i == padded else 0 if i in inactive else 3)
                for i in range(5)))
    assert tuple(map(len, configurations)) == (10, 30, 30, 20)
    flattened = [config for configs in configurations for config in configs]
    assert len(flattened) == len(set(flattened)) == 90
    complete_low_weight_signatures = set()
    for outer_class in cosets[axis].values():
        for outer in outer_class:
            options = [(0, 4) if a == "I" else (3, 5) for a in outer]
            for block_weights in product(*options):
                if sum(block_weights) <= 13:
                    complete_low_weight_signatures.add(tuple(zip(outer, block_weights)))
    assert complete_low_weight_signatures == set(flattened)
    # The five inner logical labels and weights identify the class and every
    # configuration. Within one configuration the tensor product is injective.
    # Complete disjoint inner cosets were verified above, so these signatures
    # establish distinct native words without expanding the large product.
    counts = tuple(sum(prod(len(cosets[a][w]) for a, w in config)
                       for config in configs) for configs in configurations)
    assert counts == (10000, 18000, 10800, 300000)
    block_square = tuple({sum(w * w for _, w in config) for config in configs}
                         for configs in configurations)
    assert block_square == ({27}, {43}, {59}, {43})
    assert tuple({sum(w for _, w in config) for config in configs}
                 for configs in configurations) == ({9}, {11}, {13}, {13})
    return configurations, counts


def physical_errors(single, double):
    cases = {"identity": [()], "single": [], "two_in_one_block": [],
             "two_in_distinct_blocks": []}
    for block in range(5):
        cases["single"] += [((block, error),) for error in single]
        cases["two_in_one_block"] += [((block, error),) for error in double]
    for first, second in combinations(range(5), 2):
        cases["two_in_distinct_blocks"] += [((first, a), (second, b))
                                            for a, b in product(single, repeat=2)]
    assert tuple(map(len, cases.values())) == (1, 75, 450, 2250)
    labels = set()
    for errors in cases.values():
        for error in errors:
            full = ["I"] * 25
            for block, word in error:
                full[5 * block:5 * block + 5] = word
            label = "".join(full)
            assert word_weight(label) <= 2 and label not in labels
            labels.add(label)
    assert len(labels) == 2776
    return cases


def response(configurations, error, tables):
    total = Fraction(0)
    for configuration in configurations:
        term = Fraction(1)
        for block, word in error:
            axis, w = configuration[block]
            term *= tables[axis, w, word]
        total += term
    return total / len(configurations)


def solve_square(matrix, rhs):
    n = len(matrix)
    rows = [[Fraction(value) for value in row] + [Fraction(value)]
            for row, value in zip(matrix, rhs)]
    assert all(len(row) == n + 1 for row in rows)
    for column in range(n):
        pivot = next(row for row in range(column, n) if rows[row][column])
        rows[column], rows[pivot] = rows[pivot], rows[column]
        divisor = rows[column][column]
        rows[column] = [value / divisor for value in rows[column]]
        for row in range(n):
            if row != column:
                factor = rows[row][column]
                rows[row] = [a - factor * b for a, b in zip(rows[row], rows[column])]
    return tuple(row[-1] for row in rows)


def moment_certificate():
    masses = solve_square(MOMENT_MATRIX, (1, 0, 0, 0))
    assert masses == MASSES
    signs = tuple(1 if mass > 0 else -1 for mass in masses)
    dual = solve_square(tuple(zip(*MOMENT_MATRIX)), signs)
    assert all(sum(dual[i] * MOMENT_MATRIX[i][j] for i in range(4)) == signs[j]
               for j in range(4))
    budget = sum(map(abs, masses))
    assert budget == dual[0] == Fraction(521, 8)
    scalar_bound = 1 + Fraction(2 * 9 * 13, 4)
    assert scalar_bound == Fraction(119, 2)
    assert budget - scalar_bound == Fraction(45, 8)
    table_masses = solve_square(tuple(EXPECTED_RESPONSES.values()), (1, 1, 1, 1))
    assert table_masses == masses
    return {
        "class_names": list(CLASS_NAMES),
        "rows": [list(row) for row in MOMENT_MATRIX],
        "row_meanings": ["normalization", "total_weight", "total_weight_squared",
                         "sum_of_squared_inner_block_weights"],
        "unique_class_logical_masses": [str(value) for value in masses],
        "dual_row_coefficients": [str(value) for value in dual],
        "dual_values_on_classes": list(signs),
        "code_specific_coupling_lower_bound": str(budget),
        "scalar_weight_moment_lower_bound": str(scalar_bound),
        "gap_to_scalar_bound": str(budget - scalar_bound),
    }


def main():
    if sys.version_info < (3, 10):
        raise SystemExit("Python 3.10 or later is required.")
    if sys.flags.optimize:
        raise SystemExit("Assertions are scientific checks; do not use -O.")
    group = signed_group()
    cosets, normalizer = base_cosets(group)
    single, double = inner_errors()
    tables = inner_response_tables(cosets, group, normalizer, single, double)
    cases = physical_errors(single, double)
    axes = []
    for axis in "XYZ":
        configurations, counts = class_configurations(axis, cosets)
        for kind, errors in cases.items():
            for error in errors:
                row = tuple(response(configs, error, tables) for configs in configurations)
                assert row == EXPECTED_RESPONSES[kind], (axis, kind, error, row)
                assert sum(mass * value for mass, value in zip(MASSES, row)) == 1
        assert sum(map(abs, MASSES)) == Fraction(521, 8)
        axes.append({
            "logical_axis": axis,
            "maximum_pauli_weight": 13,
            "dense_class_term_counts": dict(zip(CLASS_NAMES, counts)),
            "distinct_native_terms": sum(counts),
            "complete_low_weight_block_signatures_classified": sum(map(len, configurations)),
            "full_error_image_intertwiners_checked": sum(map(len, cases.values())),
            "error_pattern_counts": {kind: len(errors) for kind, errors in cases.items()},
            "signed_logical_membership_certified_by_tensor_factors": True,
            "coupling_budget_per_logical_half_splitting": "521/8",
        })
    print(json.dumps({
        "checks_passed": True,
        "physical_qubits": 25,
        "encoded_qubits": 1,
        "distance_from_written_concatenation_proof": 9,
        "protected_error_weight": 2,
        "base_physical_pauli_labels_enumerated": 1024,
        "base_stabilizer_count": len(group),
        "base_normalizer_count": len(normalizer),
        "base_identity_single_double_error_count": 106,
        "class_response_table": {kind: [str(value) for value in values]
                                 for kind, values in EXPECTED_RESPONSES.items()},
        "moment_certificate": moment_certificate(),
        "axes": axes,
        "sparse_optimal_construction": check_sparse_construction(),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
