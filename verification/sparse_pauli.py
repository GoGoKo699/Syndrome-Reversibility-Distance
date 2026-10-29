"""Exact recursive sparse Pauli orbit construction, using the standard library.

The code uses at most three positive seed words in each Pauli-weight class and
at most four signed seed words in each final logical Hamiltonian. Coordinatewise
base-five translations spread those seeds over at most 4*n physical Pauli terms.
No complete higher-level logical coset or numerical optimization is constructed.
"""

from fractions import Fraction
from itertools import product

from verify_readout_safe_family import axis_at, base_shells, embed, weight


def interval_coupling(distributions):
    """Couple probability lists with at most sum(k_i)-len(lists)+1 atoms.

    Lists contain (value, positive Fraction mass). Intersecting their ordered
    partitions of [0,1] preserves every one-list marginal exactly.
    """
    distributions = [list(items) for items in distributions]
    assert distributions and all(distributions)
    for items in distributions:
        assert all(mass > 0 for _, mass in items)
        assert sum(mass for _, mass in items) == 1
    indices = [0] * len(distributions)
    remaining = [items[0][1] for items in distributions]
    result, total = [], Fraction(0)
    while total < 1:
        mass = min(remaining)
        assert mass > 0
        result.append((tuple(items[index][0] for items, index in
                             zip(distributions, indices)), mass))
        total += mass
        for j, items in enumerate(distributions):
            remaining[j] -= mass
            if remaining[j] == 0 and total < 1:
                indices[j] += 1
                assert indices[j] < len(items)
                remaining[j] = items[indices[j]][1]
    assert len(result) <= sum(map(len, distributions)) - len(distributions) + 1
    return result


class _ExactConicBasis:
    """Streaming independent integer columns and exact nonnegative weights.

    The inverse uses a nonsingular square row restriction of the kept columns.
    Zero-weight basis columns are retained internally to handle degeneracy;
    they are omitted from the output. The single-error construction uses
    dimensions three and four; the two-error orbit compiler also uses this
    general reducer in larger response dimensions.
    """

    def __init__(self, dimension):
        self.dimension = dimension
        self.columns, self.rows, self.inverse = [], [], []
        self.masses, self.values = [], []

    def insert(self, value, column, mass):
        assert mass > 0 and len(column) == self.dimension
        k = len(self.columns)
        selected = [column[row] for row in self.rows]
        coordinates = [sum((entry * number for entry, number in zip(row, selected)),
                           Fraction(0)) for row in self.inverse]
        pivot, residual = None, None
        for row in range(self.dimension):
            if row in self.rows:
                continue
            residual = column[row] - sum(
                (coordinates[j] * old[row] for j, old in enumerate(self.columns)),
                Fraction(0))
            if residual:
                pivot = row
                break
        if pivot is not None:
            bottom = [sum((old[pivot] * self.inverse[j][i]
                           for j, old in enumerate(self.columns)), Fraction(0))
                      for i in range(k)]
            inverse = [
                [self.inverse[i][j] + coordinates[i] * bottom[j] / residual
                 for j in range(k)] + [-coordinates[i] / residual]
                for i in range(k)]
            inverse.append([-number / residual for number in bottom]
                           + [Fraction(1) / residual])
            self.inverse = inverse
            self.rows.append(pivot)
            self.columns.append(tuple(column))
            self.masses.append(mass)
            self.values.append(value)
            return
        assert k > 0
        step = min([mass] + [self.masses[j] / (-number)
                            for j, number in enumerate(coordinates) if number < 0])
        self.masses = [old + step * number for old, number in
                       zip(self.masses, coordinates)]
        assert all(number >= 0 for number in self.masses)
        new_mass = mass - step
        if not new_mass:
            return
        leaving = next(j for j, number in enumerate(coordinates)
                       if number < 0 and self.masses[j] == 0)
        new_row = [number / coordinates[leaving]
                   for number in self.inverse[leaving]]
        self.inverse = [new_row if i == leaving else
                        [entry - coordinates[i] * correction
                         for entry, correction in zip(row, new_row)]
                        for i, row in enumerate(self.inverse)]
        self.columns[leaving] = tuple(column)
        self.masses[leaving] = new_mass
        self.values[leaving] = value


def exact_conic_reduction(items, dimension):
    """Preserve a nonnegative integer-column sum with <=dimension atoms.

    Each item is (value, integer column, positive Fraction mass). The reducer
    retains only columns supplied by the input and verifies the full exact sum.
    """
    items = list(items)
    target = [Fraction(0)] * dimension
    basis = _ExactConicBasis(dimension)
    for value, column, mass in items:
        mass = Fraction(mass)
        assert mass > 0
        for i, number in enumerate(column):
            target[i] += number * mass
        basis.insert(value, column, mass)
    result = [(value, column, mass) for value, column, mass in
              zip(basis.values, basis.columns, basis.masses) if mass]
    assert len(result) <= dimension
    assert all(sum(column[i] * mass for _, column, mass in result) == target[i]
               for i in range(dimension))
    return result


def axis_counts(label, n):
    return tuple(sum(axis_at(label, site, n) == axis for site in range(n))
                 for axis in "XYZ")


def _add(distribution, label, sign, mass):
    if label in distribution:
        old_sign, old_mass = distribution[label]
        assert old_sign == sign
        distribution[label] = sign, old_mass + mass
    else:
        distribution[label] = sign, mass


def _reduce_class(distribution, n):
    weights = {weight(label, n) for label in distribution}
    assert len(weights) == 1
    w = next(iter(weights))
    assert w > 0 and sum(mass for _, mass in distribution.values()) == 1
    result = exact_conic_reduction(
        [((label, sign), axis_counts(label, n), mass)
         for label, (sign, mass) in sorted(distribution.items())], 3)
    answer = {label: (sign, mass) for (label, sign), _, mass in result}
    assert sum(mass for _, mass in answer.values()) == 1
    assert all(sum(axis_counts(label, n)[i] * mass for label, (_, mass) in answer.items())
               == Fraction(w, 3) for i in range(3))
    return answer


def sparse_seed_classes(level, cache=None):
    """Return axis -> weight -> seed label -> (logical-coset sign, mass).

    Each positive normalized class has <=3 seeds and mean global axis counts
    (w/3,w/3,w/3). Parent generation uses one fixed outer weight-three logical
    representative and interval coupling, not independent tensor products.
    """
    assert level >= 1
    if cache is None:
        cache = {}
    if level in cache:
        return cache[level]
    base = base_shells()
    if level == 1:
        result = {axis: {w: _reduce_class(
            {label: (sign, Fraction(1, len(terms))) for label, sign in terms.items()}, 5)
                        for w, terms in classes.items()}
                  for axis, classes in base.items()}
        cache[level] = result
        return result
    child = sparse_seed_classes(level - 1, cache)
    inner_n, inner_d = 5 ** (level - 1), 3 ** (level - 1)
    n, d = 5 * inner_n, 3 * inner_d
    result = {}
    for axis in "XYZ":
        outer, outer_sign = min(base[axis][3].items())
        active = [(block, axis_at(outer, block, 5))
                  for block in range(5) if axis_at(outer, block, 5) != "I"]
        assert len(active) == 3
        classes = {}
        for parent_weight, promoted in ((d, None), (d + 2, 0)):
            distributions = []
            for index, (_, local_axis) in enumerate(active):
                child_weight = inner_d + (2 if promoted == index else 0)
                distributions.append([((label, sign), mass)
                                      for label, (sign, mass) in sorted(
                                          child[local_axis][child_weight].items())])
            candidates = {}
            coupled = interval_coupling(distributions)
            assert len(coupled) <= 7
            for values, mass in coupled:
                label, sign = 0, outer_sign
                for (block, _), (inner_label, inner_sign) in zip(active, values):
                    label |= embed(inner_label, inner_n, n, block)
                    sign *= inner_sign
                assert weight(label, n) == parent_weight
                _add(candidates, label, sign, mass)
            classes[parent_weight] = _reduce_class(candidates, n)
        result[axis] = classes
    cache[level] = result
    return result


def sparse_seed_hamiltonian(level, axis, cache=None):
    """Return <=4 signed seeds: label -> (coset sign, relative coefficient)."""
    assert axis in "XYZ"
    n, d = 5 ** level, 3 ** level
    classes = sparse_seed_classes(level, cache)[axis]
    items = []
    for w, terms in sorted(classes.items()):
        shell_sign = 1 if w == d else -1
        scale = Fraction(d + 2, 2) if w == d else Fraction(d, 2)
        for label, (coset_sign, mass) in sorted(terms.items()):
            column = (shell_sign,) + tuple(shell_sign * c for c in axis_counts(label, n))
            items.append(((label, coset_sign, shell_sign), column, scale * mass))
    reduced = exact_conic_reduction(items, 4)
    answer = {label: (coset_sign, shell_sign * mass)
              for (label, coset_sign, shell_sign), _, mass in reduced}
    assert sum(coefficient for _, coefficient in answer.values()) == 1
    assert sum(abs(coefficient) for _, coefficient in answer.values()) == d + 1
    assert all(sum(axis_counts(label, n)[i] * coefficient
                   for label, (_, coefficient) in answer.items()) == 0 for i in range(3))
    return answer


def translation_orbit(label, level):
    """List the n=(5**level) coordinatewise base-five translated Pauli words.

    Repeated orbit words are deliberately retained: this is averaging over a
    transitive group, including its stabilizer multiplicities.
    """
    n = 5 ** level
    local_axes = [axis_at(label, site, n) for site in range(n)]
    digits = [tuple((site // 5 ** digit) % 5 for digit in range(level))
              for site in range(n)]
    for shifts in product(range(5), repeat=level):
        shifted = ["I"] * n
        for site, old_digits in enumerate(digits):
            target = sum(((old + step) % 5) * 5 ** digit
                         for digit, (old, step) in enumerate(zip(old_digits, shifts)))
            shifted[target] = local_axes[site]
        x = sum((axis in "XY") << site for site, axis in enumerate(shifted))
        z = sum((axis in "YZ") << site for site, axis in enumerate(shifted))
        yield x | (z << n)


def sparse_hamiltonian(level, axis, cache=None):
    """Expand <=4 seed orbits into <=4*n exact native Pauli coefficients.

    Values are (positive-logical-coset sign, signed coefficient relative to it).
    The literal bare Pauli coefficient is their product.
    """
    n, d = 5 ** level, 3 ** level
    answer = {}
    for label, (sign, coefficient) in sparse_seed_hamiltonian(level, axis, cache).items():
        for translated in translation_orbit(label, level):
            _add(answer, translated, sign, coefficient / n)
    answer = {label: entry for label, entry in answer.items() if entry[1]}
    assert len(answer) <= 4 * n
    assert sum(coefficient for _, coefficient in answer.values()) == 1
    assert sum(abs(coefficient) for _, coefficient in answer.values()) == d + 1
    for site in range(n):
        for local_axis in "XYZ":
            assert sum(coefficient for label, (_, coefficient) in answer.items()
                       if axis_at(label, site, n) == local_axis) == 0
    return answer


def incidence_compressed_hamiltonian(level, axis, cache=None):
    """Reduce the explicit orbit output to <=3*n+1 native Pauli terms.

    This optional full-incidence reduction uses only O(n) input columns. Its
    rational coefficients may be more complicated than the default orbit form.
    """
    n, d = 5 ** level, 3 ** level
    source = sparse_hamiltonian(level, axis, cache)
    items = []
    for label, (coset_sign, coefficient) in sorted(source.items()):
        epsilon = 1 if coefficient > 0 else -1
        column = (epsilon,) + tuple(
            epsilon * int(axis_at(label, site, n) == local_axis)
            for site in range(n) for local_axis in "XYZ")
        items.append(((label, coset_sign, epsilon), column, abs(coefficient)))
    reduced = exact_conic_reduction(items, 3 * n + 1)
    answer = {label: (coset_sign, epsilon * mass)
              for (label, coset_sign, epsilon), _, mass in reduced}
    assert len(answer) <= 3 * n + 1
    assert sum(coefficient for _, coefficient in answer.values()) == 1
    assert sum(abs(coefficient) for _, coefficient in answer.values()) == d + 1
    return answer
