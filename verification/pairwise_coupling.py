"""Exact pairwise-independent coupling of up to four finite distributions.

For three or four factors, two uniform torus coordinates provide the forms
U, V, (U+V) mod 1, and (U+2V) mod 1. Any two forms are jointly uniform:
the corresponding integer torus map has nonzero determinant. Applying the
input quantile functions therefore preserves every one- and two-factor law.

Rational quantile breakpoints cut the unit square into rational polygons.
Their exact areas give the coupling without expanding a common-denominator
grid. There are O(K**2) cells for K total input support size, since only O(K)
lines occur and the number of factors is at most four. Values may be arbitrary
Python objects; only their support indices are used as dictionary keys.

This helper constructs a joint distribution, not independent quantum noise.
It is useful when a tensor-product response depends on at most two factors.
"""

from bisect import bisect_right
from fractions import Fraction
from itertools import product


def _area(polygon):
    if len(polygon) < 3:
        return Fraction(0)
    return abs(sum((left[0] * right[1] - right[0] * left[1]
                    for left, right in zip(polygon, polygon[1:] + polygon[:1])),
                   Fraction(0))) / 2


def _clip(polygon, slope, boundary, keep_above):
    """Intersect a convex polygon with one exact rational half-plane."""
    if not polygon:
        return []
    result = []
    previous = polygon[-1]
    previous_value = previous[0] + slope * previous[1] - boundary
    previous_inside = previous_value >= 0 if keep_above else previous_value <= 0
    for current in polygon:
        current_value = current[0] + slope * current[1] - boundary
        current_inside = current_value >= 0 if keep_above else current_value <= 0
        if current_inside != previous_inside:
            fraction = previous_value / (previous_value - current_value)
            crossing = tuple(a + fraction * (b - a)
                             for a, b in zip(previous, current))
            if not result or result[-1] != crossing:
                result.append(crossing)
        if current_inside and (not result or result[-1] != current):
            result.append(current)
        previous, previous_value, previous_inside = current, current_value, current_inside
    if len(result) > 1 and result[0] == result[-1]:
        result.pop()
    return result


def _prepare(distributions):
    prepared = []
    for distribution in distributions:
        entries = [(value, Fraction(mass)) for value, mass in distribution]
        assert entries and all(mass > 0 for _, mass in entries)
        assert sum((mass for _, mass in entries), Fraction(0)) == 1
        intervals, endpoint = [], Fraction(0)
        for index, (_, mass) in enumerate(entries):
            intervals.append((endpoint, endpoint + mass, index))
            endpoint += mass
        prepared.append((entries, intervals))
    assert 1 <= len(prepared) <= 4, "One to four distributions are supported."
    return prepared


def couple(distributions):
    """Return ``[(tuple(values), positive Fraction mass), ...]`` exactly.

    Every input is a finite iterable of ``(value, positive rational mass)``
    pairs whose masses sum to one. The output preserves each input marginal
    and every pairwise product marginal, including unequal and large-denominator
    inputs. With at most four factors its support size is quadratic in the
    total input support size; it does not enumerate a denominator-sized grid.

    Ordering is deterministic for fixed input order. Duplicate input values
    are permitted and retain their separate support indices internally.
    """
    prepared = _prepare(distributions)
    if len(prepared) == 1:
        return [((value,), mass) for value, mass in prepared[0][0]]
    if len(prepared) == 2:
        return [((left, right), a * b)
                for (left, a), (right, b) in product(prepared[0][0], prepared[1][0])]

    # First two quantiles use U and V, giving a rectangular starting grid.
    cells = []
    for left, right in product(prepared[0][1], prepared[1][1]):
        x0, x1, i = left
        y0, y1, j = right
        cells.append(([(x0, y0), (x1, y0), (x1, y1), (x0, y1)], (i, j)))

    line_count = len(prepared[0][0]) + len(prepared[1][0]) - 2
    for factor in range(2, len(prepared)):
        slope = factor - 1
        # 0 <= U+slope*V <= slope+1. A period strip retains its
        # original quantile index, so separated pieces are merged later.
        bands = [(lower + period, upper + period, index)
                 for period in range(slope + 1)
                 for lower, upper, index in prepared[factor][1]]
        right_edges = [upper for _, upper, _ in bands]
        line_count += (slope + 1) * len(prepared[factor][0]) - 1
        updated = []
        for polygon, indices in cells:
            coordinates = [x + slope * y for x, y in polygon]
            minimum, maximum = min(coordinates), max(coordinates)
            first = bisect_right(right_edges, minimum)
            for band in range(first, len(bands)):
                lower, upper, index = bands[band]
                if lower >= maximum:
                    break
                cut = _clip(polygon, slope, lower, True)
                cut = _clip(cut, slope, upper, False)
                if _area(cut):
                    updated.append((cut, indices + (index,)))
        cells = updated
        assert len(cells) <= 1 + line_count * (line_count + 1) // 2

    masses = {}
    for polygon, indices in cells:
        masses[indices] = masses.get(indices, Fraction(0)) + _area(polygon)
    assert masses and all(mass > 0 for mass in masses.values())
    assert sum(masses.values(), Fraction(0)) == 1
    return [(tuple(prepared[factor][0][index][0]
                   for factor, index in enumerate(indices)), masses[indices])
            for indices in sorted(masses)]
