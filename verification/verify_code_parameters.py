"""Exhaustively verify the five- and ten-qubit algebraic examples.

All centralizer elements and all stabilizer elements are enumerated. The
distance and r values are computed directly from their definitions, not
inferred from the circuit theorem or the family construction.
"""

from collections import Counter, defaultdict
import json

from pauli_tools import (
    FIVE_GENERATORS, FIVE_MINIMUM_LOGICALS, cnot, from_string,
    local_factors, normalizer_basis, product_five_generators, rank,
    span, support, symplectic, syndrome, to_string, weight,
)


def inspect_code(generators: list[int], n: int) -> dict:
    assert rank(generators) == len(generators)
    assert all(symplectic(a, b, n) == 0 for a in generators for b in generators)
    stabilizers = span(generators)
    stabilizer_set = set(stabilizers)
    normalizer = span(normalizer_basis(generators, n))
    assert len(normalizer) == 1 << (2 * n - len(generators))
    logicals = [p for p in normalizer if p not in stabilizer_set]
    by_support = defaultdict(list)
    for g in stabilizers:
        by_support[support(g, n)].append(g)
    minimum = None
    collision_count = 0
    for logical in logicals:
        for stabilizer in by_support[support(logical, n)]:
            candidate = weight(logical ^ stabilizer, n)
            collision_count += 1
            minimum = candidate if minimum is None else min(minimum, candidate)
    return {
        "n": n,
        "k": n - len(generators),
        "stabilizer_elements": len(stabilizers),
        "normalizer_elements": len(normalizer),
        "distance": min(weight(p, n) for p in logicals),
        "r": "infinity" if minimum is None else minimum,
        "support_collision_pairs": collision_count,
    }


def main() -> dict:
    generators = [from_string(g) for g in FIVE_GENERATORS]
    group = span(generators)
    single_syndromes = [syndrome(from_string("I" * j + axis + "I" * (4 - j)), generators, 5)
                        for j in range(5) for axis in "XYZ"]
    assert set(single_syndromes) == set(range(1, 16))
    assert Counter(weight(g, 5) for g in group) == {0: 1, 4: 15}
    cosets = {}
    for axis in "XYZ":
        logical = from_string(axis * 5)
        words = [logical ^ g for g in group]
        counts = Counter(weight(p, 5) for p in words)
        assert counts == {3: 10, 5: 6}
        assert from_string(FIVE_MINIMUM_LOGICALS[axis]) in words
        cosets[axis] = dict(sorted(counts.items()))

    original = product_five_generators()
    changed = [cnot(g, 0, 5, 10) for g in original]
    product_result = inspect_code(original, 10)
    changed_result = inspect_code(changed, 10)
    assert (product_result["distance"], product_result["r"]) == (3, "infinity")
    assert (changed_result["distance"], changed_result["r"]) == (3, 4)

    g = from_string("XIYXXYXIXZ")
    logical = from_string("YIYXXXZIXX")
    witness = from_string("ZIIIIZYIIY")
    changed_group = set(span(changed))
    assert g in changed_group
    assert syndrome(logical, changed, 10) == 0 and logical not in changed_group
    assert g ^ logical == witness
    assert support(g, 10) == support(logical, 10)
    columns = [syndrome(p, changed, 10) for p in local_factors(witness, 10)]
    assert rank(columns) == 3
    assert all(symplectic(g, p, 10) == 1 for p in local_factors(witness, 10))
    # Check the circuit definition directly, independently of the rank test.
    terms = span(local_factors(witness, 10))
    assert [p for p in terms if syndrome(p, changed, 10) == 0] == [0, witness]
    return {
        "check": "exact_code_parameters",
        "five_qubit_code": inspect_code(generators, 5),
        "five_qubit_logical_weight_counts": cosets,
        "distinct_nonzero_single_qubit_syndromes": len(set(single_syndromes)),
        "product_code": product_result,
        "cnot_code": changed_result,
        "witness": {"g": to_string(g, 10), "L": to_string(logical, 10),
                    "E": to_string(witness, 10), "syndrome_column_rank": rank(columns)},
        "checks_passed": True,
    }


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True))
