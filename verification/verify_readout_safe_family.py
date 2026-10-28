#!/usr/bin/env python3
"""Exact certificates for a readout-safe concatenated-five gate family.

The identical standalone program is maintained in both companion repositories.
It enumerates all signed terms of all three logical axes at concatenation levels
one and two. At higher levels it checks the recursive counts and rational
identities; the general family, distance lower bound, and reversibility claim
are established by the written proofs, not by exhaustive large-code searches.

Only Python's standard library is used. Assertions are scientific checks, so
the program rejects Python's -O option.
"""

from fractions import Fraction
from itertools import product
import json
import sys


FIVE_GENERATORS = ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")


def pauli(word):
    n = len(word)
    x = sum((axis in "XY") << site for site, axis in enumerate(word))
    z = sum((axis in "YZ") << site for site, axis in enumerate(word))
    assert set(word) <= set("IXYZ")
    return x | (z << n)


def axis_at(label, site, n):
    return "IXZY"[((label >> site) & 1) + 2 * ((label >> (n + site)) & 1)]


def weight(label, n):
    return ((label & ((1 << n) - 1)) | (label >> n)).bit_count()


def anticommutes(left, right, n):
    mask = (1 << n) - 1
    return (((left & mask) & (right >> n)) ^
            ((left >> n) & (right & mask))).bit_count() % 2


def multiply(left, right, n):
    """Literal Hermitian P(left) P(right) = i**phase P(label)."""
    mask = (1 << n) - 1
    x, z = left & mask, left >> n
    u, v = right & mask, right >> n
    xx, zz = x ^ u, z ^ v
    phase = ((x & z).bit_count() + (u & v).bit_count()
             + 2 * (z & u).bit_count() - (xx & zz).bit_count()) % 4
    return xx | (zz << n), phase


def signed_group(generators, n):
    group = {0: 0}
    for generator in generators:
        added = {}
        for label, phase in group.items():
            out, extra = multiply(label, generator, n)
            assert out not in group and out not in added
            assert (phase + extra) % 2 == 0
            added[out] = (phase + extra) % 4
        group.update(added)
    return group


def signed_basis(generators, n):
    """Eliminate labels while retaining the phases of signed stabilizers."""
    basis = {}
    for original in generators:
        label, phase = original, 0
        while label:
            pivot = label.bit_length() - 1
            if pivot not in basis:
                basis[pivot] = (label, phase)
                break
            old, old_phase = basis[pivot]
            label, extra = multiply(label, old, n)
            phase = (phase + old_phase + extra) % 4
        assert label, "Stabilizer generators must be independent."
        assert phase % 2 == 0
    return basis


def stabilizer_phase(label, basis, n):
    """Eigenphase of a bare Pauli in the signed stabilizer label span."""
    phase = 0
    while label:
        pivot = label.bit_length() - 1
        assert pivot in basis, "Pauli is outside the stabilizer span."
        other, other_phase = basis[pivot]
        label, extra = multiply(label, other, n)
        phase = (phase + other_phase + extra) % 4
    assert phase % 2 == 0
    return phase


def embed(label, inner_n, total_n, block):
    mask = (1 << inner_n) - 1
    offset = block * inner_n
    return ((label & mask) << offset) | ((label >> inner_n) << (total_n + offset))


def lift_word(word, inner_n):
    total_n = 5 * inner_n
    return sum(embed(pauli(axis * inner_n), inner_n, total_n, block)
               for block, axis in enumerate(word) if axis != "I")


def code_generators(level):
    if level == 1:
        return tuple(map(pauli, FIVE_GENERATORS))
    inner_n = 5 ** (level - 1)
    n = 5 * inner_n
    return tuple(embed(g, inner_n, n, block)
                 for block in range(5) for g in code_generators(level - 1)) + tuple(
                     lift_word(word, inner_n) for word in FIVE_GENERATORS)


def base_shells():
    group = signed_group(tuple(map(pauli, FIVE_GENERATORS)), 5)
    assert len(group) == 16
    shells = {}
    for axis in "XYZ":
        logical = pauli(axis * 5)
        shell = {3: {}, 5: {}}
        for stabilizer, phase in group.items():
            label, extra = multiply(logical, stabilizer, 5)
            total_phase = (phase + extra) % 4
            assert total_phase % 2 == 0
            shell[weight(label, 5)][label] = (-1) ** (total_phase // 2)
        assert len(shell[3]) == 10 and len(shell[5]) == 6
        shells[axis] = shell
    return shells


def level_two_shells(axis, base):
    shells = {9: {}, 11: {}}
    for outer, outer_sign in base[axis][3].items():
        active = [(block, axis_at(outer, block, 5))
                  for block in range(5) if axis_at(outer, block, 5) != "I"]
        assert len(active) == 3
        for promoted in (None, 0, 1, 2):
            shell_weight = 9 if promoted is None else 11
            choices = [tuple(base[local_axis][5 if index == promoted else 3].items())
                       for index, (_, local_axis) in enumerate(active)]
            for pieces in product(*choices):
                label, sign = 0, outer_sign
                for (block, _), (inner, inner_sign) in zip(active, pieces):
                    label |= embed(inner, 5, 25, block)
                    sign *= inner_sign
                assert weight(label, 25) == shell_weight
                assert label not in shells[shell_weight], "Repeated physical term."
                shells[shell_weight][label] = sign
    assert len(shells[9]) == 10000 and len(shells[11]) == 18000
    return shells


def check_axis(level, axis, shells, generators, basis):
    n, d = 5 ** level, 3 ** level
    logical = pauli(axis * n)
    assert all(not anticommutes(logical, g, n) for g in generators)
    coefficients = {d: Fraction(d + 2, 2), d + 2: -Fraction(d, 2)}
    records = {}
    all_responses = {}
    for shell_weight, terms in shells.items():
        counts = [{a: 0 for a in "XYZ"} for _ in range(n)]
        for label, sign in terms.items():
            assert weight(label, n) == shell_weight
            # Independent signed generator elimination verifies each physical
            # term acts as +logical after multiplication by its stored sign.
            difference, extra = multiply(label, logical, n)
            phase = (extra + stabilizer_phase(difference, basis, n)) % 4
            assert phase % 2 == 0 and (-1) ** (phase // 2) == sign
            for site in range(n):
                local_axis = axis_at(label, site, n)
                if local_axis != "I":
                    counts[site][local_axis] += 1
        count = len(terms)
        expected = Fraction(count * shell_weight, 3 * n)
        assert expected.denominator == 1
        assert all(value == expected for row in counts for value in row.values())
        responses = {}
        for site in range(n):
            for error_axis in "XYZ":
                # A one-site Pauli reverses precisely the other two axes.
                reversed_count = sum(value for local_axis, value in counts[site].items()
                                     if local_axis != error_axis)
                responses[site, error_axis] = Fraction(count - 2 * reversed_count, count)
        expected_response = 1 - Fraction(4 * shell_weight, 3 * n)
        assert set(responses.values()) == {expected_response}
        all_responses[shell_weight] = responses
        records[str(shell_weight)] = {
            "term_count": count,
            "count_per_site_per_nonidentity_axis": int(expected),
            "normalized_single_error_response": str(expected_response),
            "logical_mass": str(coefficients[shell_weight]),
            "absolute_coefficient_per_term": str(abs(coefficients[shell_weight]) / count),
            "all_signed_logical_coset_memberships_checked": True,
        }
    assert sum(coefficients.values()) == 1
    for site in range(n):
        for error_axis in "XYZ":
            assert sum(coefficients[w] * all_responses[w][site, error_axis]
                       for w in coefficients) == 1
    budget = sum(abs(value) for value in coefficients.values())
    assert budget == d + 1
    assert budget == Fraction((d + 2) + d, (d + 2) - d)
    result = {
        "axis": axis,
        "shells": records,
        "identity_and_single_error_intertwining_certificates": 1 + 3 * n,
        "coupling_budget_per_logical_half_splitting": str(budget),
        "maximum_pauli_weight": d + 2,
        "distinct_term_count": sum(len(shell) for shell in shells.values()),
    }
    if level == 2 and axis == "Z":
        # Distance nine permits correction of more faults than this Hamiltonian
        # protects. Two X faults in the first inner block expose the distinction.
        error = pauli("XX" + "I" * 23)
        responses = {
            w: Fraction(sum((-1) ** anticommutes(label, error, n)
                            for label in terms), len(terms))
            for w, terms in shells.items()
        }
        assert responses == {9: Fraction(9, 25), 11: Fraction(89, 225)}
        response = sum(coefficients[w] * responses[w] for w in coefficients)
        assert response == Fraction(1, 5) and response != 1
        result["unprotected_two_error_witness"] = {
            "pauli": "XX" + "I" * 23,
            "one_based_sites": [1, 2],
            "normalized_shell_responses": {str(w): str(value)
                                           for w, value in responses.items()},
            "logical_rate_after_two_errors": str(response),
            "intended_logical_rate": "1",
        }
    return result


def check_recursion():
    records = []
    low, high = 10, 6
    for level in range(1, 5):
        n, d = 5 ** level, 3 ** level
        if level > 1:
            low, high = 10 * low ** 3, 30 * high * low ** 2
        assert low == 10 ** ((d - 1) // 2)
        assert high == Fraction(d, 5) * low
        low_response = 1 - Fraction(4 * d, 3 * n)
        high_response = 1 - Fraction(4 * (d + 2), 3 * n)
        a, b = Fraction(d + 2, 2), -Fraction(d, 2)
        assert a + b == 1 and a * d + b * (d + 2) == 0
        assert a * low_response + b * high_response == 1
        # A parent site is active with probability 3/5. In the upper shell,
        # exactly one of its three active children uses the child's upper shell.
        if level > 1:
            inner_n, inner_d = n // 5, d // 3
            inner_low = 1 - Fraction(4 * inner_d, 3 * inner_n)
            inner_high = 1 - Fraction(4 * (inner_d + 2), 3 * inner_n)
            assert low_response == Fraction(2, 5) + Fraction(3, 5) * inner_low
            assert high_response == Fraction(2, 5) + Fraction(3, 5) * (
                Fraction(2, 3) * inner_low + Fraction(1, 3) * inner_high)
        records.append({
            "level": level, "n": n, "distance_from_concatenation_proof": d,
            "maximum_pauli_weight": d + 2,
            "lower_shell_count": str(low), "upper_shell_count": str(high),
            "lower_shell_response": str(low_response),
            "upper_shell_response": str(high_response),
            "coupling_budget_per_logical_half_splitting": str(abs(a) + abs(b)),
        })
    return records


def main():
    if sys.version_info < (3, 10):
        raise SystemExit("Python 3.10 or later is required.")
    if sys.flags.optimize:
        raise SystemExit("Assertions are part of the certificates; do not use -O.")
    base = base_shells()
    levels = []
    for level in (1, 2):
        n, d = 5 ** level, 3 ** level
        generators = code_generators(level)
        assert all(not anticommutes(left, right, n)
                   for left in generators for right in generators)
        basis = signed_basis(generators, n)
        assert len(basis) == len(generators) == n - 1
        # Commuting even generators and odd logical X,Y,Z certify the parity
        # separation used in the written proof that r(S) is infinite.
        assert all(weight(g, n) % 2 == 0 for g in generators)
        assert all(weight(pauli(axis * n), n) % 2 == 1 for axis in "XYZ")
        axes = []
        for axis in "XYZ":
            shells = base[axis] if level == 1 else level_two_shells(axis, base)
            axes.append(check_axis(level, axis, shells, generators, basis))
        levels.append({
            "level": level, "n": n, "k": 1,
            "distance_from_concatenation_proof": d,
            "stabilizer_rank": len(basis),
            "even_stabilizer_odd_logical_parity_certificate": True,
            "axes": axes,
        })
    print(json.dumps({"checks_passed": True, "enumerated_levels": levels,
                      "recursive_arithmetic_checks": check_recursion()},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
