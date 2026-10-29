#!/usr/bin/env python3
"""Independent exact checks of the sparse two-error synthesis certificate.

All 69,125 native words are checked against signed stabilizer elimination,
for each logical axis. Integer orbit characters certify all 70,126 physical
Pauli errors of weight at most two, without a Hilbert-space matrix. Literal
translated characters independently check the bitset correlation formula.
Only Python's standard library is required; this program rejects ``-O``.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
import json
from math import lcm
from pathlib import Path
import sys

import growing_two_error_sparse as sparse
from pairwise_coupling import couple
from positive_two_error import CLASS_NAMES, RESPONSES, check_positive_classes
from verify_readout_safe_family import (
    anticommutes, axis_at, base_shells, code_generators, multiply, pauli,
    signed_basis, stabilizer_phase, weight,
)


def _negative(site, level):
    answer, stride = 0, 1
    for _ in range(level):
        site, digit = divmod(site, 5)
        answer += (-digit % 5) * stride
        stride *= 5
    return answer


def _difference(left, right, level):
    answer, stride = 0, 1
    for _ in range(level):
        left, a = divmod(left, 5)
        right, b = divmod(right, 5)
        answer += ((b - a) % 5) * stride
        stride *= 5
    return answer


def _features(level):
    n = 5 ** level
    return (("identity",),) + tuple(("single", a) for a in range(3)) + tuple(
        ("pair", delta, a, b) for delta in range(1, n)
        if delta < _negative(delta, level) for a in range(3) for b in range(3))


def _error_coverage(level=3):
    """Classify every distinct physical label, independently of the loader."""
    n = 5 ** level
    features = _features(level)
    assert features == sparse.error_features(level)
    assert all(_negative(site, level) == sparse.negative_site(site, level)
               for site in range(n))
    index = {feature: row for row, feature in enumerate(features)}
    coverage = Counter({0: 1})
    kinds = Counter(identity=1)
    for site in range(n):
        for a in range(3):
            coverage[index["single", a]] += 1
            kinds["single"] += 1
    for first, second in combinations(range(n), 2):
        delta = _difference(first, second, level)
        reverse = delta > _negative(delta, level)
        if reverse:
            delta = _negative(delta, level)
        kind = ("same_bottom_block" if first // 5 == second // 5 else
                "same_middle_block" if first // 25 == second // 25 else
                "different_middle_blocks")
        for a, b in product(range(3), repeat=2):
            feature = ("pair", delta, b, a) if reverse else ("pair", delta, a, b)
            coverage[index[feature]] += 1
            kinds[kind] += 1
    assert len(index) == 562 and set(coverage) == set(range(562))
    assert coverage[0] == 1 and set(coverage[row] for row in range(1, 562)) == {125}
    assert kinds == {"identity": 1, "single": 375, "same_bottom_block": 2250,
                     "same_middle_block": 11250, "different_middle_blocks": 56250}
    assert sum(coverage.values()) == 1 + 3 * n + 9 * n * (n - 1) // 2 == 70126
    return {"physical_error_labels": sum(coverage.values()),
            "translation_orbits": len(features), "counts_by_pattern": dict(kinds),
            "nonidentity_error_orbit_size": n}


def _literal_translate(label, shifts, level):
    """Site-by-site permutation, independent of the compiler's bit shifts."""
    n = 5 ** level
    axes = ["I"] * n
    for site in range(n):
        digits, target, stride = site, 0, 1
        for shift in shifts:
            digits, digit = divmod(digits, 5)
            target += ((digit + shift) % 5) * stride
            stride *= 5
        axes[target] = axis_at(label, site, n)
    return pauli("".join(axes))


def _literal_orbit_column(label, level):
    n = 5 ** level
    words = [_literal_translate(label, shifts, level)
             for shifts in product(range(5), repeat=level)]
    assert words == list(sparse.translation_orbit(label, level))
    result = []
    axes = [pauli(axis + "I" * (n - 1)) for axis in "XYZ"]
    for feature in _features(level):
        if feature[0] == "identity":
            error = 0
        elif feature[0] == "single":
            error = axes[feature[1]]
        else:
            _, delta, a, b = feature
            word = ["I"] * n
            word[0], word[delta] = "XYZ"[a], "XYZ"[b]
            error = pauli("".join(word))
        result.append(sum(1 - 2 * anticommutes(word, error, n) for word in words))
    return tuple(result)


def _pattern(label, n=125):
    """Check bottom counts and the active ternary support at both outer levels."""
    bottom = [sum(axis_at(label, 5 * block + site, n) != "I" for site in range(5))
              for block in range(n // 5)]
    assert set(bottom) <= {0, 3, 4, 5}
    active = [value % 2 for value in bottom]
    assert sorted(sum(active[start:start + 5]) for start in range(0, 25, 5)) == [0, 0, 3, 3, 3]
    counts = Counter(bottom)
    possibilities = {"A0": (9, 0, 0), "A1": (8, 0, 1),
                     "A2": (7, 0, 2), "F": (9, 1, 0)}
    return next(name for name, expected in possibilities.items()
                if tuple(counts[w] for w in (3, 4, 5)) == expected)


def check_native_certificate():
    data_path = Path(__file__).resolve().parent / "data" / "two_error_125_orbits.json"
    data = json.loads(data_path.read_text())
    assert data["n"] == 125 and data["logical_axis"] == "Z"
    assert data["orbit_count"] == 553 and data["native_pauli_terms"] == 69125
    assert data["interaction_order"] == 31 and Fraction(data["J"]) == Fraction(3491, 8)
    expected = {"A0": Fraction(899, 8), "A1": -Fraction(837, 4),
                "A2": Fraction(1701, 16), "F": -Fraction(135, 16)}
    basis = signed_basis(code_generators(3), 125)
    assert len(basis) == 124
    reports = []
    for axis in "XYZ":
        logical = pauli(axis * 125)
        seeds = sparse.hamiltonian_seeds(3, axis)
        assert len(seeds) == 553
        denominator = lcm(*(mass.denominator for _, _, mass, _ in seeds))
        totals = [0] * 562
        masses = defaultdict(Fraction)
        expanded = {}
        literal_checks = []
        # Cover every class and both coefficient signs with literal characters.
        selected = {next(i for i, atom in enumerate(seeds) if atom[3] == name)
                    for name in CLASS_NAMES}
        for number, (label, sign, mass, name) in enumerate(seeds):
            assert sign in (-1, 1) and mass and (mass > 0) == (expected[name] > 0)
            assert _pattern(label) == name == sparse.class_pattern(label, 3)
            assert weight(label, 125) == {"A0": 27, "A1": 29, "A2": 31, "F": 31}[name]
            masses[name] += mass
            column = sparse.integer_orbit_column(label, 3)
            assert len(column) == 562 and column[0] == 125
            numerator = int(mass * denominator)
            totals = [old + numerator * value for old, value in zip(totals, column)]
            if number in selected:
                assert _literal_orbit_column(label, 3) == column
                literal_checks.append(name)
            orbit = tuple(sparse.translation_orbit(label, 3))
            assert len(orbit) == len(set(orbit)) == 125
            for word in orbit:
                assert word not in expanded, "Different stored orbit atoms overlap."
                # Full bare Pauli membership, including the Hermitian sign.
                difference, phase = multiply(word, logical, 125)
                phase = (phase + stabilizer_phase(difference, basis, 125)) % 4
                assert phase % 2 == 0 and (-1) ** (phase // 2) == sign
                assert weight(word, 125) == weight(label, 125)
                expanded[word] = sign, mass / 125
            # Generator translations establish invariance of each orbit.
            orbit_set = set(orbit)
            for shifts in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
                assert {sparse.translate(word, shifts, 3) for word in orbit} == orbit_set
        assert dict(masses) == expected
        assert all(value == 125 * denominator for value in totals)
        assert len(expanded) == 69125
        native = sparse.load_hamiltonian(3, axis)
        assert native == expanded
        assert sum(coefficient for _, coefficient in native.values()) == 1
        assert sum(abs(coefficient) for _, coefficient in native.values()) == Fraction(3491, 8)
        assert max(weight(label, 125) for label in native) == 31
        native_denominator = lcm(*(coefficient.denominator for _, coefficient in native.values()))
        reports.append({"logical_axis": axis, "orbit_atoms": len(seeds),
                        "distinct_native_words": len(native),
                        "signed_native_coset_memberships": len(native),
                        "exact_integer_response_equations": len(totals),
                        "represented_physical_errors": 70126,
                        "literal_character_columns_checked": sorted(literal_checks),
                        "orbit_mass_common_denominator_bits": denominator.bit_length(),
                        "native_coefficient_common_denominator_bits": native_denominator.bit_length(),
                        "signed_class_masses": {name: str(masses[name]) for name in CLASS_NAMES},
                        "interaction_order": 31, "J_per_logical_half_splitting": "3491/8"})
    return {"all_three_axes_verified": True, "axes": reports,
            "full_error_coverage": _error_coverage()}


def _check_marginals(distributions, joint):
    assert all(mass > 0 for _, mass in joint)
    assert sum(mass for _, mass in joint) == 1
    for factor, law in enumerate(distributions):
        actual = defaultdict(Fraction)
        for values, mass in joint:
            actual[values[factor]] += mass
        assert dict(actual) == dict(law)
    pairs = 0
    for left, right in combinations(range(len(distributions)), 2):
        actual = defaultdict(Fraction)
        for values, mass in joint:
            actual[values[left], values[right]] += mass
        target = {(a, b): x * y for a, x in distributions[left] for b, y in distributions[right]}
        assert dict(actual) == target
        pairs += len(target)
    return pairs


def _tensor_responses(distributions, joint):
    """Check all physical errors, including two errors inside one factor."""
    n = 2 * len(distributions)
    words = [(pauli("".join(values)), mass) for values, mass in joint]
    errors = [()]
    errors.extend(((site, axis),) for site in range(n) for axis in "XYZ")
    errors.extend(((left, a), (right, b)) for left, right in combinations(range(n), 2)
                  for a, b in product("XYZ", repeat=2))
    same_factor = 0
    for locations in errors:
        error_word = ["I"] * n
        for site, axis in locations:
            error_word[site] = axis
        error = pauli("".join(error_word))
        actual = sum(mass * (1 - 2 * anticommutes(word, error, n)) for word, mass in words)
        expected = Fraction(1)
        for factor, law in enumerate(distributions):
            local_error = pauli("".join(error_word[2 * factor:2 * factor + 2]))
            expected *= sum(mass * (1 - 2 * anticommutes(pauli(word), local_error, 2))
                            for word, mass in law)
        assert actual == expected
        if len(locations) == 2 and locations[0][0] // 2 == locations[1][0] // 2:
            same_factor += 1
    return len(errors), same_factor


def check_coupling():
    denominator = (1 << 127) - 1
    laws = [[("IX", Fraction(denominator - 17, denominator)), ("ZY", Fraction(17, denominator))],
            [("YZ", Fraction(1, 7)), ("XI", Fraction(2, 7)), ("ZZ", Fraction(4, 7))],
            [("XX", Fraction(3, 11)), ("IY", Fraction(8, 11))],
            [("ZI", Fraction(2, 13)), ("YY", Fraction(5, 13)), ("XZ", Fraction(6, 13))]]
    reports = []
    for count in (3, 4):
        distributions = laws[:count]
        # Input generators exercise one-pass iterables; equalities use originals.
        joint = couple((iter(law) for law in distributions))
        pair_entries = _check_marginals(distributions, joint)
        equations, same_factor = _tensor_responses(distributions, joint)
        assert joint == couple(distributions), "The constructor must be deterministic."
        reports.append({"factors": count, "input_support_sizes": list(map(len, distributions)),
                        "largest_input_denominator_bits": denominator.bit_length(),
                        "joint_support": len(joint), "pair_marginal_entries_checked": pair_entries,
                        "physical_tensor_response_equations": equations,
                        "two_errors_in_same_factor_equations": same_factor})
    # Deterministic factors and matching quantile/wrap boundaries are separate edges.
    edges = [[[("fixed", Fraction(1))], [("a", Fraction(1, 2)), ("b", Fraction(1, 2))],
              [("c", Fraction(1, 2)), ("d", Fraction(1, 2))],
              [("e", Fraction(1, 3)), ("f", Fraction(2, 3))]],
             [[(0, Fraction(1, 2)), (1, Fraction(1, 2))]] * 4]
    for distributions in edges:
        _check_marginals(distributions, couple(distributions))
    assert couple([[("only", Fraction(1))]]) == [(("only",), Fraction(1))]
    two = couple(laws[:2])
    assert two == [((a, b), x * y) for a, x in laws[0] for b, y in laws[1]]
    return {"exact_pairwise_coupling_verified": True, "unequal_laws": reports,
            "deterministic_factor_and_wrap_boundary_cases": len(edges),
            "one_and_two_factor_cases": True}


@lru_cache(maxsize=None)
def _recursive_profile(level, name, sites):
    if not sites:
        return Fraction(1)
    if level == 2:
        pattern = ("single" if len(sites) == 1 else
                   "same_block" if sites[0] // 5 == sites[1] // 5 else "different_blocks")
        return RESPONSES[pattern][CLASS_NAMES.index(name)]
    child_n = 5 ** (level - 1)
    result = Fraction(0)
    for outer in base_shells()["Z"][3]:
        active = [block for block in range(5) if axis_at(outer, block, 5) != "I"]
        inactive = [block for block in range(5) if block not in active]
        for kinds, pad, probability in sparse.recurrence_schemes(level, name):
            value = probability / 10
            for block, kind in zip(active, kinds):
                local = tuple(site % child_n for site in sites if site // child_n == block)
                value *= _recursive_profile(level - 1, kind, local)
            if pad is not None:
                local = tuple(site % child_n for site in sites if site // child_n == inactive[pad])
                hit_bottom_blocks = len({site // 5 for site in local})
                value *= 1 - Fraction(16 * hit_bottom_blocks, 15 * (child_n // 5))
            result += value
    return result


def _dense_profile(level, name, sites):
    k, bottom_n = 3 ** (level - 1), 5 ** (level - 1)
    p = Fraction(k, bottom_n)
    if not sites:
        return Fraction(1)
    same = len(sites) == 2 and sites[0] // 5 == sites[1] // 5
    low = Fraction(-1, 15) if same else Fraction(1, 5)
    high = Fraction(1, 9) if same else -Fraction(1, 3)
    if len(sites) == 1 or same:
        if name == "F":
            return 1 + p * (low - 1) - Fraction(16, 15 * bottom_n)
        promotions = CLASS_NAMES.index(name)
        return 1 - p + p * (low + Fraction(promotions, k) * (high - low))
    left, right = sites[0] // 5, sites[1] // 5
    height = level - 1
    digits = [(left // 5 ** r % 5, right // 5 ** r % 5) for r in reversed(range(height))]
    divergence = next(r for r, (a, b) in enumerate(digits, 1) if a != b)
    rho = Fraction(3, 5) ** (divergence - 1) * Fraction(3, 10) * Fraction(3, 5) ** (2 * (height - divergence))
    low, high, padding = Fraction(1, 5), -Fraction(1, 3), -Fraction(1, 15)
    if name == "F":
        return (1 + 2 * p * (low - 1) + Fraction(2, bottom_n) * (padding - 1)
                + rho * ((low - 1) ** 2 + 2 * (low - 1) * (padding - 1) / (k - 1)))
    promotions = CLASS_NAMES.index(name)
    mean = low + Fraction(promotions, k) * (high - low)
    joint = (low ** 2 + Fraction(2 * promotions, k) * low * (high - low)
             + Fraction(promotions * (promotions - 1), k * (k - 1)) * (high - low) ** 2)
    return 1 + 2 * p * (mean - 1) + rho * (1 - 2 * mean + joint)


def check_recurrence():
    reports = []
    for level in range(3, 7):
        k_child, k = 3 ** (level - 2), 3 ** (level - 1)
        checks = 0
        locations = [(), (0,), (0, 1)] + [(0, 5 ** r) for r in range(1, level)]
        for name in CLASS_NAMES:
            schemes = list(sparse.recurrence_schemes(level, name))
            assert all(probability > 0 for _, _, probability in schemes)
            assert sum(probability for _, _, probability in schemes) == 1
            for sites in locations:
                assert _recursive_profile(level, name, sites) == _dense_profile(level, name, sites)
                checks += 1
        schemes = list(sparse.recurrence_schemes(level, "A2"))
        same = sum(probability for kinds, _, probability in schemes if kinds[0] == "A2")
        different = sum(probability for kinds, _, probability in schemes
                        if kinds[0] == kinds[1] == "A1")
        assert same / (k_child * (k_child - 1) // 2) == Fraction(2, k * (k - 1))
        assert different / (k_child ** 2) == Fraction(2, k * (k - 1))
        frontier = list(sparse.recurrence_schemes(level, "F"))
        assert sum(q for kinds, pad, q in frontier if kinds[0] == "F") / (k_child - 1) == Fraction(1, k - 1)
        assert [q for _, pad, q in frontier if pad is not None] == [Fraction(1, k - 1)] * 2
        reports.append({"level": level, "n": 5 ** level,
                        "exact_class_response_equations": checks,
                        "uniform_promotion_pair_and_frontier_laws": True})
    return {"recurrence_profiles_verified": True, "levels": reports}


def check_reference_compiler():
    """Actually run the exact final reduction, without the stored 125-qubit gate."""
    seeds = sparse.synthesize_hamiltonian(2, "Z")
    assert 0 < len(seeds) <= len(_features(2)) == 112
    denominator = lcm(*(mass.denominator for _, _, mass, _ in seeds))
    totals = [0] * 112
    masses = defaultdict(Fraction)
    native = {}
    basis = signed_basis(code_generators(2), 25)
    logical = pauli("Z" * 25)
    for label, sign, mass, name in seeds:
        assert sparse.class_pattern(label, 2) == name
        assert weight(label, 25) in (9, 11, 13)
        masses[name] += mass
        numerator = int(mass * denominator)
        column = sparse.integer_orbit_column(label, 2)
        totals = [old + numerator * value for old, value in zip(totals, column)]
        for word in sparse.translation_orbit(label, 2):
            old_sign, old_mass = native.get(word, (sign, Fraction(0)))
            assert old_sign == sign
            native[word] = sign, old_mass + mass / 25
    assert all(total == 25 * denominator for total in totals)
    assert dict(masses) == {"A0": Fraction(143, 8), "A1": -Fraction(117, 4),
                            "A2": Fraction(243, 16), "F": -Fraction(45, 16)}
    assert sum(mass for _, mass in native.values()) == 1
    assert sum(abs(mass) for _, mass in native.values()) == Fraction(521, 8)
    assert len(native) <= 25 * len(seeds)
    for label, (sign, _) in native.items():
        difference, phase = multiply(label, logical, 25)
        phase = (phase + stabilizer_phase(difference, basis, 25)) % 4
        assert phase % 2 == 0 and (-1) ** (phase // 2) == sign
    return {"reference_end_to_end_compiler_executed": True, "level": 2,
            "n": 25, "logical_axis": "Z", "orbit_atoms": len(seeds),
            "distinct_native_words": len(native),
            "signed_native_coset_memberships": len(native),
            "exact_integer_response_equations": len(totals),
            "represented_physical_errors": 2776, "interaction_order": 13,
            "J_per_logical_half_splitting": "521/8"}


def main():
    if sys.version_info < (3, 10):
        raise SystemExit("Python 3.10 or later is required.")
    if sys.flags.optimize:
        raise SystemExit("Assertions are scientific checks; do not use -O.")
    result = {"checks_passed": True,
              "positive_depth_two_classes": check_positive_classes(),
              "pairwise_coupling": check_coupling(),
              "recursive_class_profiles": check_recurrence(),
              "executed_reference_compiler": check_reference_compiler(),
              "native_depth_three_certificate": check_native_certificate()}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
