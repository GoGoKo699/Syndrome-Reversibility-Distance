#!/usr/bin/env python3
"""Exact three-error gates at interaction order fifteen on 25 qubits.

Literal signed five-qubit Pauli algebra classifies every contributing block
pattern. Tensor-factor character sums check all 64,876 physical Pauli errors
of weight at most three on each logical axis. The 1,480,960-term Hamiltonian
is specified by seven disjoint positive logical-coset classes; it is not
expanded as a dense native list or a 2**25 Hilbert-space matrix. Signed local
factor identities prove membership of the entire construction, with separate
full-code checks of one native word from each block configuration.

Only Python's standard library is required. Do not disable assertions.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
import json
from math import lcm, prod
import sys

import verify_readout_safe_family as family


GENERATORS = ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
CLASS_NAMES = ("A0", "A1", "A2", "A3", "F0", "F1", "B5")
MASSES = tuple(map(Fraction, ("715/16", "-1755/16", "3645/32", "-1485/32",
                              "-675/32", "675/32", "-3/2")))
COUNTS = (10000, 18000, 10800, 2160, 300000, 540000, 600000)
J = Fraction(2867, 8)
OPTIONS = (("I", 0), ("I", 4)) + tuple((a, w) for a in "XYZ" for w in (3, 5))
OPTION_INDEX = {option: index for index, option in enumerate(OPTIONS)}
MOMENT_NAMES = ("1", "W", "W^2", "W^3", "R2", "W*R2", "R3")
MOMENT_COLUMNS = (
    (1, 9, 81, 729, 27, 243, 81),
    (1, 11, 121, 1331, 43, 473, 179),
    (1, 13, 169, 2197, 59, 767, 277),
    (1, 15, 225, 3375, 75, 1125, 375),
    (1, 13, 169, 2197, 43, 559, 145),
    (1, 15, 225, 3375, 59, 885, 243),
    (1, 15, 225, 3375, 45, 675, 135),
)
DUAL = tuple(map(Fraction, ("2867/8", "-2375/24", "7", "-1/6",
                            "47/8", "-1/8", "-1/2")))
PROFILE_VALUES = {
    "identity": ("1",) * 7,
    "single_equivalent": ("13/25", "31/75", "23/75", "1/5", "23/75", "1/5", "1/5"),
    "double_in_one_block_equivalent": ("9/25", "89/225", "97/225", "7/15", "11/75", "41/225", "-1/15"),
    "two_distinct_blocks": ("29/125", "13/125", "1/225", "-1/15", "23/375", "-11/1125", "1/25"),
    "normalizer_triple_in_one_block": ("1/5",) * 6 + ("-1/3",),
    "two_plus_one_blocks": ("17/125", "121/1125", "47/675", "1/45", "1/125", "-1/675", "-1/75"),
    "three_distinct_blocks": ("53/625", "-1/75", "-11/225", "-1/27", "1/125", "-59/5625", "1/125"),
}
PROFILES = {name: tuple(map(Fraction, values)) for name, values in PROFILE_VALUES.items()}
EXPECTED_PROFILE_COUNTS = {
    "identity": 1, "single_equivalent": 375,
    "double_in_one_block_equivalent": 1350, "two_distinct_blocks": 2250,
    "normalizer_triple_in_one_block": 150, "two_plus_one_blocks": 27000,
    "three_distinct_blocks": 33750,
}


def multiply(left, right):
    """Literal P(left)P(right)=i**phase P(word), retaining Hermitian phases."""
    assert len(left) == len(right)
    forward = {("X", "Y"): "Z", ("Y", "Z"): "X", ("Z", "X"): "Y"}
    word, phase = [], 0
    for a, b in zip(left, right):
        if a == "I":
            word.append(b)
        elif b == "I":
            word.append(a)
        elif a == b:
            word.append("I")
        elif (a, b) in forward:
            word.append(forward[a, b])
            phase += 1
        else:
            word.append(forward[b, a])
            phase -= 1
    return "".join(word), phase % 4


def anticommutes(left, right):
    return sum(a != "I" and b != "I" and a != b
               for a, b in zip(left, right)) % 2


def weight(word):
    return len(word) - word.count("I")


def signed_cosets():
    group = {"IIIII": 1}
    for generator in GENERATORS:
        added = {}
        for word, sign in group.items():
            label, phase = multiply(word, generator)
            assert phase % 2 == 0 and label not in group and label not in added
            added[label] = sign * (-1) ** (phase // 2)
        group.update(added)
    assert len(group) == 16
    for left, left_sign in group.items():
        for right, right_sign in group.items():
            result, phase = multiply(left, right)
            assert not anticommutes(left, right) and phase % 2 == 0
            assert group[result] == left_sign * right_sign * (-1) ** (phase // 2)
    assert all(weight(word) == 4 for word in group if word != "IIIII")
    cosets = {"I": {0: {"IIIII": 1},
                     4: {word: sign for word, sign in group.items() if word != "IIIII"}}}
    seen = set(group)
    for axis in "XYZ":
        cosets[axis] = {3: {}, 5: {}}
        for stabilizer, sign in group.items():
            word, phase = multiply(axis * 5, stabilizer)
            assert phase % 2 == 0 and word not in seen
            assert all(not anticommutes(word, g) for g in GENERATORS)
            seen.add(word)
            cosets[axis][weight(word)][word] = sign * (-1) ** (phase // 2)
        assert tuple(map(len, cosets[axis].values())) == (10, 6)
    normalizer = {"".join(letters) for letters in product("IXYZ", repeat=5)
                  if all(not anticommutes("".join(letters), g) for g in GENERATORS)}
    assert normalizer == seen and len(seen) == 64
    assert min(weight(word) for word in normalizer - set(group)) == 3
    return group, cosets


def configurations(axis, cosets):
    answer = [[] for _ in CLASS_NAMES]
    for outer in cosets[axis][3]:
        active = [i for i, a in enumerate(outer) if a != "I"]
        inactive = [i for i in range(5) if i not in active]
        for promotions in range(4):
            for promoted in combinations(active, promotions):
                answer[promotions].append(tuple(
                    (a, 0 if a == "I" else 5 if i in promoted else 3)
                    for i, a in enumerate(outer)))
        for padded in inactive:
            for promotions in range(2):
                for promoted in combinations(active, promotions):
                    answer[4 + promotions].append(tuple(
                        (a, 4 if i == padded else 0 if a == "I" else 5 if i in promoted else 3)
                        for i, a in enumerate(outer)))
    for outer in cosets[axis][5]:
        answer[6].append(tuple((a, 3) for a in outer))
    assert tuple(map(len, answer)) == (10, 30, 30, 10, 20, 60, 6)
    flat = [config for configs in answer for config in configs]
    assert len(flat) == len(set(flat)) == 166
    # All normalizer words with nontrivial compression must induce an outer
    # logical Pauli. Exhaust all signed outer words and every inner shell.
    complete = set()
    for shell in cosets[axis].values():
        for outer in shell:
            choices = [(0, 4) if a == "I" else (3, 5) for a in outer]
            for weights in product(*choices):
                if sum(weights) <= 15:
                    complete.add(tuple(zip(outer, weights)))
    assert complete == set(flat)
    counts = tuple(sum(prod(len(cosets[a][w]) for a, w in config)
                       for config in configs) for configs in answer)
    assert counts == COUNTS and sum(counts) == 1480960
    for configs, column in zip(answer, MOMENT_COLUMNS):
        # Uniform configuration and factor averages give equal mass to every
        # native word of a class because each configuration has this same size.
        assert len({prod(len(cosets[a][w]) for a, w in config)
                    for config in configs}) == 1
        for config in configs:
            weights = [w for _, w in config]
            w, r2, r3 = sum(weights), sum(x * x for x in weights), sum(x ** 3 for x in weights)
            assert (1, w, w * w, w ** 3, r2, w * r2, r3) == column
    return answer


def _local_errors():
    words = []
    for size in range(4):
        for sites in combinations(range(5), size):
            for axes in product("XYZ", repeat=size):
                word = ["I"] * 5
                for site, axis in zip(sites, axes):
                    word[site] = axis
                words.append("".join(word))
    assert len(words) == len(set(words)) == 376
    return words


def character_tables(cosets):
    table = {}
    for error in _local_errors():
        column = []
        for axis, w in OPTIONS:
            words = cosets[axis][w]
            numerator = sum(1 - 2 * anticommutes(word, error) for word in words)
            assert 30 % len(words) == 0
            column.append(numerator * (30 // len(words)))
        table[error] = tuple(column)
    assert table["IIIII"] == (30,) * 8
    assert all(table[e][1] == -2 for e in table if weight(e) in (1, 2))
    return table


def physical_errors(table):
    records, seen = [], set()
    physical_counts = Counter()
    for size in range(4):
        for sites in combinations(range(25), size):
            for axes in product("XYZ", repeat=size):
                word = ["I"] * 25
                for site, axis in zip(sites, axes):
                    word[site] = axis
                word = "".join(word)
                label = family.pauli(word)
                assert label not in seen and weight(word) == size
                seen.add(label)
                factors = tuple(word[5 * b:5 * b + 5] for b in range(5))
                signature = tuple(table[local] for local in factors)
                pattern = tuple(sorted(weight(local) for local in factors if weight(local)))
                normalizer_triple = (pattern == (3,) and all(
                    not anticommutes(next(local for local in factors if weight(local)), g)
                    for g in GENERATORS))
                physical_counts[pattern] += 1
                records.append((signature, pattern, normalizer_triple))
    assert len(seen) == len(records) == 64876
    assert dict(physical_counts) == {(): 1, (1,): 75, (2,): 450, (1, 1): 2250,
                                     (3,): 1350, (1, 2): 27000, (1, 1, 1): 33750}
    return records, physical_counts


def factorized_profile(signature, encoded_configs):
    """Exact finite tensor sums; cacheable equal local character columns."""
    active = [block for block, values in enumerate(signature) if values != (30,) * 8]
    assert len(active) <= 3
    common = 60 * 30 ** 3
    result = []
    for configs in encoded_configs:
        total = 0
        for config in configs:
            term = 1
            for block in active:
                term *= signature[block][config[block]]
            total += term
        assert 60 % len(configs) == 0
        numerator = total * (60 // len(configs)) * 30 ** (3 - len(active))
        result.append(Fraction(numerator, common))
    return tuple(result)


def check_configuration_signs(axis, configs, cosets):
    """Independent full-code spot check per shell configuration, not every term."""
    basis = family.signed_basis(family.code_generators(2), 25)
    logical = family.pauli(axis * 25)
    checked = 0
    for group in configs:
        for config in group:
            outer = "".join(a for a, _ in config)
            sign = cosets[axis][weight(outer)][outer]
            pieces = []
            for local_axis, w in config:
                word = min(cosets[local_axis][w])
                sign *= cosets[local_axis][w][word]
                pieces.append(word)
            label = family.pauli("".join(pieces))
            difference, phase = family.multiply(label, logical, 25)
            phase = (phase + family.stabilizer_phase(difference, basis, 25)) % 4
            assert phase % 2 == 0 and (-1) ** (phase // 2) == sign
            checked += 1
    assert checked == 166
    return checked


def check_axis(axis, cosets, records):
    configs = configurations(axis, cosets)
    signs_checked = check_configuration_signs(axis, configs, cosets)
    encoded = [[tuple(OPTION_INDEX[option] for option in config) for config in group]
               for group in configs]
    names = {profile: name for name, profile in PROFILES.items()}
    cache, profile_counts, triple_counts = {}, Counter(), Counter()
    for signature, physical_pattern, normalizer_triple in records:
        if signature not in cache:
            cache[signature] = factorized_profile(signature, encoded)
        profile = cache[signature]
        assert profile in names
        name = names[profile]
        assert sum(mass * response for mass, response in zip(MASSES, profile)) == 1
        profile_counts[name] += 1
        if physical_pattern == (3,):
            assert (name == "normalizer_triple_in_one_block") == normalizer_triple
            triple_counts[name] += 1
    assert dict(profile_counts) == EXPECTED_PROFILE_COUNTS
    assert triple_counts == {"single_equivalent": 300,
                             "double_in_one_block_equivalent": 900,
                             "normalizer_triple_in_one_block": 150}
    return {"logical_axis": axis, "physical_error_identities_checked": len(records),
            "distinct_response_profiles": len(profile_counts),
            "local_character_signatures_evaluated": len(cache),
            "profile_counts": dict(profile_counts),
            "same_bottom_block_triple_profile_counts": dict(triple_counts),
            "complete_logical_block_configurations": sum(map(len, configs)),
            "full_code_configuration_sign_spot_checks": signs_checked,
            "native_terms_specified_by_disjoint_factor_classes": sum(COUNTS),
            "native_expansion_enumerated": False}


def _determinant(matrix):
    work = [list(map(Fraction, row)) for row in matrix]
    result = Fraction(1)
    for column in range(len(work)):
        pivot = next(row for row in range(column, len(work)) if work[row][column])
        if pivot != column:
            work[pivot], work[column] = work[column], work[pivot]
            result = -result
        value = work[column][column]
        result *= value
        for row in range(column + 1, len(work)):
            ratio = work[row][column] / value
            work[row] = [a - ratio * b for a, b in zip(work[row], work[column])]
    return result


def check_moments():
    matrix = tuple(zip(*MOMENT_COLUMNS))
    assert _determinant(matrix) == 2949120
    assert tuple(sum(x * mass for x, mass in zip(row, MASSES)) for row in matrix) == (1, 0, 0, 0, 0, 0, 0)
    signs = tuple(1 if mass > 0 else -1 for mass in MASSES)
    assert tuple(sum(y * x for y, x in zip(DUAL, column)) for column in MOMENT_COLUMNS) == signs
    assert sum(map(abs, MASSES)) == DUAL[0] == J
    shell_masses = {w: sum(mass for mass, col in zip(MASSES, MOMENT_COLUMNS) if col[1] == w)
                    for w in (9, 11, 13, 15)}
    assert tuple(shell_masses.values()) == tuple(map(Fraction, ("715/16", "-1755/16", "1485/16", "-429/16")))
    scalar = sum(map(abs, shell_masses.values()))
    assert scalar == 274 and J - scalar == Fraction(675, 8)
    native = tuple(mass / count for mass, count in zip(MASSES, COUNTS))
    assert sum(count * value for count, value in zip(COUNTS, native)) == 1
    assert sum(count * abs(value) for count, value in zip(COUNTS, native)) == J
    denominator = lcm(*(value.denominator for value in native))
    assert native == tuple(map(Fraction, ("143/32000", "-39/6400", "27/2560", "-11/512", "-9/128000", "1/25600", "-1/400000")))
    assert denominator == 3200000 and denominator.bit_length() == 22
    return {"moment_rows": MOMENT_NAMES, "moment_columns": MOMENT_COLUMNS,
            "moment_matrix_determinant": 2949120,
            "signed_class_masses": dict(zip(CLASS_NAMES, map(str, MASSES))),
            "dual_row": list(map(str, DUAL)), "dual_column_values": signs,
            "optimal_native_budget_per_logical_half_splitting": str(J),
            "scalar_shell_relaxation": str(scalar),
            "block_moment_improvement": str(J - scalar),
            "native_counts_by_class": dict(zip(CLASS_NAMES, COUNTS)),
            "relative_native_coefficients_by_class": dict(zip(CLASS_NAMES, map(str, native))),
            "native_coefficient_common_denominator": denominator,
            "native_coefficient_denominator_bits": denominator.bit_length(),
            "bare_coefficient_convention": "relative coefficient times signed logical-coset sign"}


def check_fixed_frame_direction():
    """Check the stronger budget in the fixed transversal X,Y,Z frame.

    Each nontrivial normalizer Pauli compresses to a signed member of exactly
    one of these three logical cosets. Thus a pattern with compressed operator
    C_j=a_j*(Omega_X X+Omega_Y Y+Omega_Z Z) costs at least
    |a_j|*sum|Omega|, by the triangle inequality separately in each coset.
    The three supplied axis Hamiltonians have disjoint logical cosets and
    attain the sum. Detectable/scalar terms cannot reduce a native l1 budget.
    This statement fixes the Pauli frame; it is not Bloch-rotation invariance.
    """
    omega = tuple(map(Fraction, ("3/5", "0", "4/5")))
    assert sum(value * value for value in omega) == 1
    l1 = sum(map(abs, omega))
    pattern_axis_coefficients = tuple(tuple(mass * value for value in omega)
                                      for mass in MASSES)
    matrix = tuple(zip(*MOMENT_COLUMNS))
    for axis in range(3):
        assert tuple(sum(row[j] * pattern_axis_coefficients[j][axis]
                         for j in range(7)) for row in matrix) == (omega[axis], 0, 0, 0, 0, 0, 0)
    lower = sum(sum(map(abs, coefficients)) for coefficients in pattern_axis_coefficients)
    attained = sum(abs(value) * sum(map(abs, MASSES)) for value in omega)
    assert lower == attained == J * l1 == Fraction(20069, 40)
    return {"canonical_logical_frame": "transversal X^25,Y^25,Z^25",
            "exact_budget_formula": "(2867/8)*(abs(Omega_X)+abs(Omega_Y)+abs(Omega_Z))",
            "example_axis_coefficients": list(map(str, omega)),
            "example_logical_half_splitting_squared": "1",
            "example_native_budget": str(attained),
            "normalizer_cosets_disjoint": True,
            "bloch_frame_invariance_claimed": False,
            "error_identity_extension": "linearity of the three separately verified full error-image identities"}


def main():
    if sys.version_info < (3, 10):
        raise SystemExit("Python 3.10 or later is required.")
    if sys.flags.optimize:
        raise SystemExit("Assertions are scientific checks; do not use -O.")
    group, cosets = signed_cosets()
    table = character_tables(cosets)
    records, physical_counts = physical_errors(table)
    result = {"checks_passed": True, "n": 25, "d": 9,
              "protected_error_weight": 3, "interaction_order": 15,
              "literal_signed_base_group_size": len(group),
              "complete_inner_normalizer_size": 64,
              "inner_error_character_tables_checked": len(table),
              "physical_error_pattern_counts": {str(pattern): count for pattern, count in physical_counts.items()},
              "profile_response_vectors": {name: list(map(str, values)) for name, values in PROFILES.items()},
              "moment_certificate": check_moments(),
              "fixed_frame_mixed_axis_budget": check_fixed_frame_direction(),
              "axes": [check_axis(axis, cosets, records) for axis in "XYZ"]}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
