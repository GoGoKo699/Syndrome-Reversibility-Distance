"""Exact truncated-tensor construction of growing three-error gates.

The compact recursion retains total increment degree at most the protection
order. Responses are evaluated over exact rational tensor factors, without a
native Pauli expansion. An optional lazy iterator specifies the complete
Hamiltonian; above the base it yields uncollected contributions and can be
very expensive. No polynomial native-enumeration runtime is claimed.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import product
import json

import verify_three_error_gates as base
from verify_readout_safe_family import embed, pauli


CLASS_NAMES = base.CLASS_NAMES
POLARITIES = (1, -1, 1, -1, -1, 1, -1)
RANKS = (0, 1, 2, 3, 2, 3, 3)


def class_masses(level, protection=3):
    """Signed logical masses in the seven bottom-weight patterns.

    At higher levels a pattern's positive distribution need not be the same
    in different increments. These numbers are masses, not such distributions.
    """
    assert level >= 2 and 0 <= protection <= 3
    d = 3 ** level
    if protection == 0:
        values = (1, 0, 0, 0, 0, 0, 0)
    elif protection == 1:
        values = (Fraction(d + 2, 2), -Fraction(d, 2), 0, 0, 0, 0, 0)
    elif protection == 2:
        values = (Fraction((d + 2) * (d + 4), 8), -Fraction(d * (d + 4), 4),
                  Fraction(d * (2 * d + 9), 16), 0, -Fraction(5 * d, 16), 0, 0)
    else:
        values = (Fraction((d + 2) * (d + 4) * (d + 6), 48),
                  -Fraction(d * (d + 4) * (d + 6), 16),
                  Fraction(d * (2 * d + 9) * (d + 6), 32),
                  -Fraction(d * (2 * d + 15) * (d + 6), 96),
                  -Fraction(5 * d * (d + 6), 32),
                  Fraction(5 * d * (d + 6), 32), -Fraction(d, 6))
    return dict(zip(CLASS_NAMES, map(Fraction, values)))


def increment_coefficients(level=2):
    """Return the signed class-mass differences for increments 0,1,2,3."""
    gates = [class_masses(level, t) for t in range(4)]
    return tuple(tuple(gates[t][name] - (gates[t - 1][name] if t else 0)
                       for name in CLASS_NAMES) for t in range(4))


def increment_budgets(level):
    d = 3 ** level
    return (Fraction(1), Fraction(d), Fraction(d * (4 * d + 13), 8),
            Fraction(d * (4 * d * d + 39 * d + 107), 24))


def budget(level, protection=3):
    assert level >= 2 and 0 <= protection <= 3
    return sum(increment_budgets(level)[:protection + 1], Fraction(0))


def coefficient_denominator_bound(level):
    """A common denominator for every collected native coefficient.

    This is a safe bound, not a minimal denominator. It applies simultaneously
    to all four protection orders and increments on each logical Pauli axis.
    """
    assert level >= 2
    denominator = 3200000
    for _ in range(2, level):
        denominator = 10 * denominator ** 3
    return denominator


def moment_columns(level):
    d = 3 ** level
    offsets = ((0, 0, 0), (2, 16, 98), (4, 32, 196), (6, 48, 294),
               (4, 16, 64), (6, 32, 162), (6, 18, 54))
    return tuple((1, d + a, (d + a) ** 2, (d + a) ** 3,
                  3 * d + b, (d + a) * (3 * d + b), 9 * d + c)
                 for a, b, c in offsets)


def dual(level):
    d = 3 ** level
    return (budget(level), -Fraction(12 * d * d + 111 * d + 404, 24),
            Fraction(d + 5, 2), -Fraction(1, 6), Fraction(d + 38, 8),
            -Fraction(1, 8), -Fraction(1, 2))


@lru_cache(maxsize=None)
def _base_data():
    _, cosets = base.signed_cosets()
    table = base.character_tables(cosets)
    configurations = {axis: base.configurations(axis, cosets) for axis in "XYZ"}
    encoded = {axis: tuple(tuple(tuple(base.OPTION_INDEX[option] for option in config)
                                for config in group) for group in configurations[axis])
               for axis in "XYZ"}
    return cosets, table, configurations, encoded


@lru_cache(maxsize=16384)
def base_profile(axis, error):
    """Evaluate all seven base classes on an actual <=3-site Pauli error."""
    _, table, _, encoded = _base_data()
    word = ["I"] * 25
    for site, local_axis in error:
        word[site] = local_axis
    signature = tuple(table["".join(word[5 * block:5 * block + 5])]
                      for block in range(5))
    return base.factorized_profile(signature, encoded[axis])


def _validate_error(level, error):
    assert level >= 2 and len(error) <= 3
    assert len({site for site, _ in error}) == len(error)
    assert all(0 <= site < 5 ** level and axis in "XYZ" for site, axis in error)


@lru_cache(maxsize=None)
def increment_terms(level, index, axis):
    """Parent tensor descriptors: (root word, root sign, three indices).

    Each descriptor has scalar coefficient root_sign/10. Active children
    are the nonidentity root factors in increasing physical block order.
    """
    assert level >= 3 and 0 <= index <= 3 and axis in "XYZ"
    cosets, _, _, _ = _base_data()
    return tuple((outer, sign, indices)
                 for outer, sign in sorted(cosets[axis][3].items())
                 for indices in product(range(index + 1), repeat=3)
                 if sum(indices) == index)


@lru_cache(maxsize=65536)
def increment_response(level, index, axis, error):
    """Exact scalar f in Delta_index E V = f E V axis.

    This evaluator computes the tensor sum. It deliberately does not return
    zero merely because error weight is below the increment index; that
    vanishing property is independently checked by the verifier.
    """
    assert 0 <= index <= 3 and axis in "XYZ"
    _validate_error(level, error)
    if level == 2:
        return sum((mass * response for mass, response in zip(
            increment_coefficients(2)[index], base_profile(axis, error))), Fraction(0))
    child_n = 5 ** (level - 1)
    restrictions = tuple(tuple((site % child_n, local_axis)
                               for site, local_axis in error if site // child_n == block)
                         for block in range(5))
    answer = Fraction(0)
    for outer, _, indices in increment_terms(level, index, axis):
        active = tuple((block, local_axis) for block, local_axis in enumerate(outer)
                       if local_axis != "I")
        response = Fraction(1, 10)
        for (block, local_axis), child_index in zip(active, indices):
            response *= increment_response(level - 1, child_index, local_axis,
                                           restrictions[block])
        answer += response
    return answer


def gate_response(level, protection, axis, error):
    _validate_error(level, error)
    assert 0 <= protection <= 3 and axis in "XYZ"
    return sum((increment_response(level, index, axis, error)
                for index in range(protection + 1)), Fraction(0))


def _base_native(coefficients, axis):
    """Yield actual bare Pauli contributions for uniform depth-two classes."""
    cosets, _, configurations, _ = _base_data()
    for mass, configs in zip(coefficients, configurations[axis]):
        if not mass:
            continue
        for config in configs:
            outer = "".join(local_axis for local_axis, _ in config)
            sign = cosets[axis][base.weight(outer)][outer]
            choices = [tuple(cosets[local_axis][w].items()) for local_axis, w in config]
            count = 1
            for terms in choices:
                count *= len(terms)
            coefficient = mass / (len(configs) * count)
            for pieces in product(*choices):
                word, actual_sign = [], sign
                for piece, factor_sign in pieces:
                    word.append(piece)
                    actual_sign *= factor_sign
                yield pauli("".join(word)), actual_sign * coefficient


def increment_native_contributions(level, index, axis="Z"):
    """Lazy full specification; duplicate words are not collected above m=2."""
    assert level >= 2 and 0 <= index <= 3 and axis in "XYZ"
    if level == 2:
        yield from _base_native(increment_coefficients(2)[index], axis)
        return
    n, child_n = 5 ** level, 5 ** (level - 1)
    for outer, sign, indices in increment_terms(level, index, axis):
        active = tuple((block, local_axis) for block, local_axis in enumerate(outer)
                       if local_axis != "I")

        def tensor(position, label, coefficient):
            if position == 3:
                yield label, coefficient
                return
            block, local_axis = active[position]
            for child_word, child_coefficient in increment_native_contributions(
                    level - 1, indices[position], local_axis):
                yield from tensor(position + 1,
                                  label | embed(child_word, child_n, n, block),
                                  coefficient * child_coefficient)

        yield from tensor(0, 0, Fraction(sign, 10))


def native_contributions(level, protection=3, axis="Z"):
    """Iterate the entire specified gate, using exact bare coefficients.

    Above depth two, contributions may repeat a physical word. Sign coherence
    means repeats never cancel; summing repeated coefficients gives the native
    Hamiltonian. Materializing the full list is intentionally optional.
    """
    assert level >= 2 and 0 <= protection <= 3 and axis in "XYZ"
    if level == 2:
        yield from _base_native(tuple(class_masses(2, protection).values()), axis)
    else:
        for index in range(protection + 1):
            yield from increment_native_contributions(level, index, axis)


if __name__ == "__main__":
    import sys
    if sys.flags.optimize:
        raise SystemExit("Assertions are scientific checks; do not use -O.")
    print(json.dumps({"construction": "total-degree truncated increment tensor sum",
                      "examples": [{"level": level, "n": 5 ** level, "d": 3 ** level,
                                    "interaction_order": 3 ** level + 6,
                                    "native_budget_per_axis": str(budget(level))}
                                   for level in range(2, 7)]}, indent=2))
