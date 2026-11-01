"""Exact protected logical Pauli strings on fixed product codes.

For a specified nonidentity logical word on k Q_m blocks, total-degree
truncation combines the active blocks' increment hierarchies. Exact response
evaluation avoids the enormous dense native expansion. The optional lazy
native iterator yields contributions, with repeated words left uncollected.
"""

from collections import defaultdict
from fractions import Fraction
from itertools import combinations_with_replacement
import json

import three_error_recursion as single
from verify_readout_safe_family import embed


def active_blocks(logical_word):
    assert isinstance(logical_word, str) and logical_word
    assert set(logical_word) <= set("IXYZ")
    result = tuple((block, axis) for block, axis in enumerate(logical_word) if axis != "I")
    assert result, "The target logical Pauli must be nonidentity."
    return result


def target_distance(level, logical_word):
    """Minimum native Pauli weight in this logical coset, not code distance."""
    assert level >= 2
    return len(active_blocks(logical_word)) * 3 ** level


def interaction_order(level, logical_word, protection=3):
    assert 0 <= protection <= 3
    return target_distance(level, logical_word) + 2 * protection


def budget(level, logical_word, protection=3):
    """Exact budget of the constructed unit logical Pauli Hamiltonian."""
    assert 0 <= protection <= 3
    distance = target_distance(level, logical_word)
    increments = (Fraction(1), Fraction(distance),
                  Fraction(distance * (4 * distance + 13), 8),
                  Fraction(distance * (4 * distance * distance + 39 * distance + 107), 24))
    return sum(increments[:protection + 1], Fraction(0))


def increment_components(logical_word, protection=3):
    """Yield active-block increment tuples with total degree <=protection.

    Multisets generate only the retained weak compositions, rather than an
    exponential Cartesian list that would mostly be discarded.
    """
    assert 0 <= protection <= 3
    count = len(active_blocks(logical_word))
    for degree in range(protection + 1):
        for selected in combinations_with_replacement(range(count), degree):
            indices = [0] * count
            for block in selected:
                indices[block] += 1
            yield tuple(indices)


def _validate_error(level, logical_word, error):
    assert level >= 2 and len(error) <= 3
    active_blocks(logical_word)
    n = len(logical_word) * 5 ** level
    assert len({site for site, _ in error}) == len(error)
    assert all(0 <= site < n and axis in ("X", "Y", "Z") for site, axis in error)


def gate_response(level, logical_word, protection, error):
    """Return the exact scalar in G_P E V = scalar * E V P.

    The actual increment tensors are evaluated, including inactive-block
    errors. No assertion of transparency is substituted for this calculation.
    """
    assert 0 <= protection <= 3
    _validate_error(level, logical_word, error)
    child_n = 5 ** level
    active = active_blocks(logical_word)
    restrictions = tuple(tuple((site % child_n, axis) for site, axis in error
                               if site // child_n == block) for block, _ in active)
    answer = Fraction(0)
    for indices in increment_components(logical_word, protection):
        term = Fraction(1)
        for (_, axis), local_error, index in zip(active, restrictions, indices):
            term *= single.increment_response(level, index, axis, local_error)
            if not term:
                break
        answer += term
    return answer


def coefficient_denominator_bound(level, logical_word):
    """Safe common denominator for normalized axis-string coefficients."""
    return single.coefficient_denominator_bound(level) ** len(active_blocks(logical_word))


def seed_pair_coefficients(protection=3):
    """Collected positive-class pair masses at m=2 for two active blocks.

    This 7-by-7 sparse rational table is independent of the two prescribed
    logical axes. Its class averages are uniform, so dividing an entry by
    the two class sizes gives the relative coefficient of each native word.
    """
    assert 0 <= protection <= 3
    increments = single.increment_coefficients(2)
    answer = defaultdict(Fraction)
    for left in range(protection + 1):
        for right in range(protection + 1 - left):
            for i, a in enumerate(increments[left]):
                for j, b in enumerate(increments[right]):
                    if a and b:
                        answer[single.CLASS_NAMES[i], single.CLASS_NAMES[j]] += a * b
    return {pair: mass for pair, mass in sorted(answer.items()) if mass}


def native_contributions(level, logical_word, protection=3):
    """Lazy full Hamiltonian specification, without collecting duplicates.

    The coefficient is the bare Hermitian Pauli coefficient. Enumerating this
    iterator can be extremely expensive; compact response evaluation and
    resource certificates do not require a dense native expansion.
    """
    assert level >= 2 and 0 <= protection <= 3
    active = active_blocks(logical_word)
    child_n = 5 ** level
    n = len(logical_word) * child_n
    for indices in increment_components(logical_word, protection):

        def tensor(position, label, coefficient):
            if position == len(active):
                yield label, coefficient
                return
            block, axis = active[position]
            for word, local_coefficient in single.increment_native_contributions(
                    level, indices[position], axis):
                yield from tensor(position + 1,
                                  label | embed(word, child_n, n, block),
                                  coefficient * local_coefficient)

        yield from tensor(0, 0, Fraction(1))


if __name__ == "__main__":
    import sys
    if sys.flags.optimize:
        raise SystemExit("Assertions are scientific checks; do not use -O.")
    examples = []
    for level, word in ((2, "ZZ"), (3, "ZZ"), (2, "XYZ"), (2, "ZIZ")):
        examples.append({"level": level, "logical_word": word,
                         "physical_qubits": len(word) * 5 ** level,
                         "ordinary_code_distance": 3 ** level,
                         "target_coset_minimum_weight": target_distance(level, word),
                         "interaction_order": interaction_order(level, word),
                         "native_budget": str(budget(level, word))})
    print(json.dumps({"construction": "total-degree truncated active-block increments",
                      "examples": examples}, indent=2))
