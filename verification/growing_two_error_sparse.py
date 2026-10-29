"""Exact sparse two-error synthesis and the stored 125-qubit Hamiltonian.

The deterministic compiler is a deliberately expensive reference algorithm.
It couples child classes pairwise, averages coordinate translations, and
reduces positive orbit columns exactly. The practical depth-three certificate
was found by a numerical support search on the same recursive candidates;
its stored rational coefficients require no numerical package to verify.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
import json
from pathlib import Path

from pairwise_coupling import couple
from positive_two_error import CLASS_NAMES, seed_classes
from sparse_pauli import _ExactConicBasis
from verify_readout_safe_family import (
    axis_at, base_shells, embed, pauli, signed_group, FIVE_GENERATORS,
)


def negative_site(site, level):
    return sum((-(site // 5 ** digit) % 5) % 5 * 5 ** digit
               for digit in range(level))


@lru_cache(maxsize=None)
def error_features(level):
    """One representative per coordinate-translation orbit of errors <=2."""
    n = 5 ** level
    deltas = tuple(d for d in range(1, n) if d < negative_site(d, level))
    return (("identity",),) + tuple(("single", a) for a in range(3)) + tuple(
        ("pair", delta, a, b) for delta in deltas
        for a in range(3) for b in range(3))


@lru_cache(maxsize=None)
def _wrap_masks(level):
    n = 5 ** level
    return tuple(tuple(sum(1 << site for site in range(n)
                           if (site // 5 ** digit) % 5 >= 5 - shift)
                       for shift in range(5)) for digit in range(level))


def translate_bits(bits, shifts, level):
    """Permute site bits by adding base-five coordinates, without carries."""
    for digit, shift in enumerate(shifts):
        if shift:
            stride = 5 ** digit
            wrap = _wrap_masks(level)[digit][shift]
            bits = ((bits & ~wrap) << (shift * stride)) | (
                (bits & wrap) >> ((5 - shift) * stride))
    return bits


def translate(label, shifts, level):
    n = 5 ** level
    mask = (1 << n) - 1
    return translate_bits(label & mask, shifts, level) | (
        translate_bits(label >> n, shifts, level) << n)


def translation_orbit(label, level):
    """Yield the n translated words, retaining stabilizer multiplicities."""
    return (translate(label, shifts, level)
            for shifts in product(range(5), repeat=level))


def integer_orbit_column(label, level):
    """Integer sums of all n characters, using exact bitset correlations."""
    n = 5 ** level
    x, z = label & ((1 << n) - 1), label >> n
    anti = z, x ^ z, x
    counts = tuple(bits.bit_count() for bits in anti)
    answer = [n] + [n - 2 * count for count in counts]
    for delta in range(1, n):
        if delta >= negative_site(delta, level):
            continue
        shifts = tuple((-(delta // 5 ** digit)) % 5 for digit in range(level))
        moved = tuple(translate_bits(bits, shifts, level) for bits in anti)
        answer.extend(n - 2 * counts[a] - 2 * counts[b]
                      + 4 * (anti[a] & moved[b]).bit_count()
                      for a in range(3) for b in range(3))
    assert len(answer) == len(error_features(level))
    return tuple(answer)


def cycle_axes(label, n, repetitions=1):
    for _ in range(repetitions):
        x, z = label & ((1 << n) - 1), label >> n
        label = (x ^ z) | (x << n)
    return label


def class_pattern(label, level):
    """Identify the class from actual physical weights in bottom blocks."""
    n, k = 5 ** level, 3 ** (level - 1)
    counts = {3: 0, 4: 0, 5: 0}
    for block in range(n // 5):
        w = sum(axis_at(label, 5 * block + site, n) != "I" for site in range(5))
        assert w in (0, 3, 4, 5)
        if w:
            counts[w] += 1
    possibilities = {"A0": (k, 0, 0), "A1": (k - 1, 0, 1),
                     "A2": (k - 2, 0, 2), "F": (k, 1, 0)}
    pattern = tuple(counts[w] for w in (3, 4, 5))
    return next(name for name, expected in possibilities.items() if pattern == expected)


def class_masses(level):
    d = 3 ** level
    return dict(zip(CLASS_NAMES, (Fraction((d + 2) * (d + 4), 8),
                                  -Fraction(d * (d + 4), 4),
                                  Fraction(d * (2 * d + 9), 16),
                                  -Fraction(5 * d, 16))))


def recurrence_schemes(level, name):
    """Yield (three child class names, optional root pad, probability).

    The root pad is specified by one of the two inactive-child indices; all
    localized bottom-block stabilizers inside that child are averaged later.
    """
    assert level >= 3 and name in CLASS_NAMES
    k_child = 3 ** (level - 2)
    if name == "A0":
        yield ("A0",) * 3, None, Fraction(1)
    elif name == "A1":
        for child in range(3):
            kinds = ["A0"] * 3
            kinds[child] = "A1"
            yield tuple(kinds), None, Fraction(1, 3)
    elif name == "A2":
        for child in range(3):
            kinds = ["A0"] * 3
            kinds[child] = "A2"
            yield tuple(kinds), None, Fraction(k_child - 1, 3 * (3 * k_child - 1))
        for left, right in combinations(range(3), 2):
            kinds = ["A0"] * 3
            kinds[left] = kinds[right] = "A1"
            yield tuple(kinds), None, Fraction(2 * k_child, 3 * (3 * k_child - 1))
    else:
        for child in range(3):
            kinds = ["A0"] * 3
            kinds[child] = "F"
            yield tuple(kinds), None, Fraction(k_child - 1, 3 * k_child - 1)
        for inactive in range(2):
            yield ("A0",) * 3, inactive, Fraction(1, 3 * k_child - 1)


def expand_classes(orbits, level, axis="Z"):
    """Expand positive class orbit seeds to normalized native distributions."""
    n = 5 ** level
    repetitions = {"Z": 0, "X": 1, "Y": 2}[axis]
    answer = {}
    for name, seeds in orbits.items():
        atoms = {}
        for label, sign, mass in seeds:
            label = cycle_axes(label, n, repetitions)
            for word in translation_orbit(label, level):
                old_sign, old_mass = atoms.get(word, (sign, Fraction(0)))
                assert old_sign == sign
                atoms[word] = sign, old_mass + mass / n
        answer[name] = [((word, sign), mass)
                        for word, (sign, mass) in sorted(atoms.items())]
    return answer


def recursive_candidates(level, name, child_orbits):
    """Stream exact parent atoms using pairwise-independent child factors.

    This avoids full three-/four-factor Cartesian products. It can still be
    expensive: its inputs are expanded native distributions, not seed lists.
    """
    n, child_n = 5 ** level, 5 ** (level - 1)
    child = {axis: expand_classes(child_orbits, level - 1, axis) for axis in "XYZ"}
    group = signed_group(tuple(map(pauli, FIVE_GENERATORS)), 5)
    localized = [((embed(word, 5, child_n, block), (-1) ** (phase // 2)),
                  Fraction(1, 15 * (child_n // 5)))
                 for block in range(child_n // 5) for word, phase in group.items() if word]
    outer_terms = base_shells()["Z"][3]
    for outer, outer_sign in sorted(outer_terms.items()):
        active = [(block, axis_at(outer, block, 5)) for block in range(5)
                  if axis_at(outer, block, 5) != "I"]
        inactive = [block for block in range(5) if axis_at(outer, block, 5) == "I"]
        for kinds, padded, probability in recurrence_schemes(level, name):
            distributions = [child[axis][kind] for (_, axis), kind in zip(active, kinds)]
            blocks = [block for block, _ in active]
            if padded is not None:
                distributions.append(localized)
                blocks.append(inactive[padded])
            for values, mass in couple(distributions):
                label, sign = 0, outer_sign
                for block, (word, child_sign) in zip(blocks, values):
                    label |= embed(word, child_n, n, block)
                    sign *= child_sign
                assert class_pattern(label, level) == name
                yield label, sign, probability * mass / len(outer_terms)


def synthesize_classes(level, cache=None):
    """Deterministic polynomial-time reference synthesis, potentially slow.

    Returns positive orbit seeds, at most 4+9*(n-1)/2 per class. Only depth
    two is loaded from stored data. At larger depths every candidate is
    generated by the recurrence and reduced using exact rational arithmetic.
    No practical large-depth runtime guarantee is made.
    """
    assert level >= 2
    if cache is None:
        cache = {}
    if level in cache:
        return cache[level]
    if level == 2:
        cache[level] = seed_classes("Z")
        return cache[level]
    child = synthesize_classes(level - 1, cache)
    classes = {}
    dimension = len(error_features(level))
    for name in CLASS_NAMES:
        basis = _ExactConicBasis(dimension)
        for label, sign, mass in recursive_candidates(level, name, child):
            basis.insert((label, sign), integer_orbit_column(label, level), mass)
        atoms = tuple((label, sign, mass)
                      for (label, sign), mass in zip(basis.values, basis.masses) if mass)
        assert len(atoms) <= dimension and sum(mass for _, _, mass in atoms) == 1
        classes[name] = atoms
    cache[level] = classes
    return classes


def synthesize_hamiltonian(level, axis="Z", cache=None):
    """Reference end-to-end compiler, returning at most C2(n) signed orbits.

    This calls the potentially expensive class compiler. It does not load the
    practical depth-three Hamiltonian. Optimal class masses are preserved by
    the complete error-response columns and the four-class moment theorem.
    """
    assert level >= 2 and axis in "XYZ"
    classes = synthesize_classes(level, cache)
    dimension = len(error_features(level))
    basis = _ExactConicBasis(dimension)
    scales = class_masses(level)
    for name in CLASS_NAMES:
        polarity = 1 if scales[name] > 0 else -1
        for label, sign, mass in classes[name]:
            column = tuple(polarity * value for value in integer_orbit_column(label, level))
            basis.insert((label, sign, polarity, name), column, abs(scales[name]) * mass)
    repetitions = {"Z": 0, "X": 1, "Y": 2}[axis]
    result = tuple((cycle_axes(label, 5 ** level, repetitions), sign, polarity * mass, name)
                   for (label, sign, polarity, name), mass in zip(basis.values, basis.masses)
                   if mass)
    assert len(result) <= dimension
    assert {name: sum(mass for _, _, mass, cls in result if cls == name)
            for name in CLASS_NAMES} == scales
    return result


def hamiltonian_seeds(level=3, axis="Z"):
    """Load the practical exact certificate: (label,coset sign,relative mass,class)."""
    assert level == 3 and axis in "XYZ"
    data = json.loads((Path(__file__).resolve().parent / "data" /
                       "two_error_125_orbits.json").read_text())
    assert data["n"] == 125 and data["logical_axis"] == "Z"
    repetitions = {"Z": 0, "X": 1, "Y": 2}[axis]
    masses = class_masses(level)
    answer = []
    for entry in data["orbits"]:
        label = cycle_axes(pauli(entry["pauli"]), 125, repetitions)
        name = entry["class"]
        mass = Fraction(entry["mass"])
        assert mass > 0
        answer.append((label, entry["coset_sign"], mass if masses[name] > 0 else -mass, name))
    return tuple(answer)


def load_hamiltonian(level=3, axis="Z"):
    """Return label -> (positive logical-coset sign,relative coefficient).

    The bare Hermitian Pauli coefficient is the product of the two values.
    """
    n = 5 ** level
    answer = {}
    for label, sign, mass, _ in hamiltonian_seeds(level, axis):
        for word in translation_orbit(label, level):
            old_sign, old_mass = answer.get(word, (sign, Fraction(0)))
            assert old_sign == sign
            answer[word] = sign, old_mass + mass / n
    return answer
