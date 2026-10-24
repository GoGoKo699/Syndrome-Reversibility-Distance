"""Exact structural checks for the one-CNOT family at t=1,2,3.

This checks generator rank, commutation, logical cosets, prescribed weights,
an untouched minimum-weight representative, and a transformed circuit.
It does not enumerate large-code distances. Their lower bounds are proved
in science/code_families.md.
"""

import json

from pauli_tools import (
    FIVE_GENERATORS, FIVE_MINIMUM_LOGICALS, cnot, embed, from_string,
    in_span, local_factors, rank, span_coordinates, support, symplectic,
    syndrome, to_string, weight,
)


def lift(word: str, inner: dict, representatives: dict[str, int]) -> int:
    answer = 0
    for block, axis in enumerate(word):
        if axis != "I":
            answer ^= embed(representatives[axis], inner["n"], 5 * inner["n"], block * inner["n"])
    return answer


def concatenated_five(level: int) -> dict:
    if level == 1:
        x, z = from_string("XXXXX"), from_string("ZZZZZ")
        return {"n": 5, "S": [from_string(g) for g in FIVE_GENERATORS],
                "logical": {"X": x, "Z": z, "Y": x ^ z},
                "minimum": {a: from_string(s) for a, s in FIVE_MINIMUM_LOGICALS.items()}}
    inner = concatenated_five(level - 1)
    n = 5 * inner["n"]
    generators = [embed(g, inner["n"], n, block * inner["n"])
                  for block in range(5) for g in inner["S"]]
    generators += [lift(g, inner, inner["logical"]) for g in FIVE_GENERATORS]
    return {
        "n": n,
        "S": generators,
        "logical": {a: lift(a * 5, inner, inner["logical"]) for a in "XYZ"},
        "minimum": {a: lift(word, inner, inner["minimum"])
                    for a, word in FIVE_MINIMUM_LOGICALS.items()},
    }


def check_level(level: int) -> dict:
    code = concatenated_five(level)
    n, generators, minimum = code["n"], code["S"], code["minimum"]
    distance = 3 ** level
    assert len(generators) == rank(generators) == n - 1
    assert all(symplectic(a, b, n) == 0 for a in generators for b in generators)
    assert all(weight(g, n) % 2 == 0 for g in generators)
    for axis in "XYZ":
        logical = minimum[axis]
        assert syndrome(logical, generators, n) == 0
        assert weight(logical, n) == distance
        assert in_span(logical ^ code["logical"][axis], generators)
        assert not in_span(logical, generators)
    assert support(minimum["X"], n) & 1 == 0
    assert to_string(minimum["Z"], n)[0] == "Z"

    product_generators = [embed(g, n, 2 * n, offset)
                          for offset in (0, n) for g in generators]
    changed = [cnot(g, 0, n, 2 * n) for g in product_generators]
    witness = cnot(embed(minimum["Z"], n, 2 * n, n), 0, n, 2 * n)
    untouched = embed(minimum["X"], n, 2 * n, 0)
    assert cnot(untouched, 0, n, 2 * n) == untouched
    assert weight(untouched, 2 * n) == distance
    assert syndrome(untouched, changed, 2 * n) == 0 and not in_span(untouched, changed)
    assert weight(witness, 2 * n) == distance + 1
    assert syndrome(witness, changed, 2 * n) == 0 and not in_span(witness, changed)
    factors = local_factors(witness, 2 * n)
    columns = [syndrome(p, changed, 2 * n) for p in factors]
    assert rank(columns) == len(columns) - 1
    # Produce a stabilizer whose local anticommutations are all one.
    rows = [sum(symplectic(g, p, 2 * n) << j for j, p in enumerate(factors)) for g in changed]
    selected = span_coordinates((1 << len(factors)) - 1, rows)
    g = 0
    for j, generator in enumerate(changed):
        if (selected >> j) & 1:
            g ^= generator
    logical = g ^ witness
    assert support(g, 2 * n) == support(logical, 2 * n)
    assert syndrome(logical, changed, 2 * n) == 0 and not in_span(logical, changed)
    return {"t": level, "n": 2 * n, "k": 2, "distance_proved_in_text": distance,
            "stabilizer_rank": rank(changed), "untouched_logical_weight": weight(untouched, 2 * n),
            "even_logical_circuit_weight": weight(witness, 2 * n),
            "circuit_syndrome_rank": rank(columns),
            "support_collision_product_weight": weight(g ^ logical, 2 * n),
            "control": 0, "target": n}


def main() -> dict:
    return {"check": "concatenated_family_structure",
            "levels": [check_level(t) for t in (1, 2, 3)],
            "checks_passed": True}


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True))
