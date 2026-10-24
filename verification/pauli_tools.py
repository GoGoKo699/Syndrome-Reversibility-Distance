"""Exact binary Pauli arithmetic for the verification scripts.

Bits 0,...,n-1 store X components; bits n,...,2*n-1 store Z
components. Physical qubit 0 is the leftmost printed character and the
least significant computational-basis bit. Overall phases are discarded.
This module does not construct physical syndrome projectors.
"""

from __future__ import annotations


FIVE_GENERATORS = ("XZZXI", "IXZZX", "XIXZZ", "ZXIXZ")
FIVE_MINIMUM_LOGICALS = {"X": "IYYIX", "Z": "ZYIIY", "Y": "YXIIX"}


def from_string(word: str) -> int:
    n = len(word)
    if any(char not in "IXYZ" for char in word):
        raise ValueError("A Pauli word must contain only I, X, Y, Z")
    x = sum((char in "XY") << j for j, char in enumerate(word))
    z = sum((char in "YZ") << j for j, char in enumerate(word))
    return x | (z << n)


def to_string(pauli: int, n: int) -> str:
    alphabet = "IXZY"
    return "".join(
        alphabet[((pauli >> j) & 1) | (((pauli >> (n + j)) & 1) << 1)]
        for j in range(n)
    )


def support(pauli: int, n: int) -> int:
    return (pauli & ((1 << n) - 1)) | (pauli >> n)


def weight(pauli: int, n: int) -> int:
    return support(pauli, n).bit_count()


def symplectic(a: int, b: int, n: int) -> int:
    mask = (1 << n) - 1
    return (((a & mask) & (b >> n)).bit_count()
            + ((a >> n) & (b & mask)).bit_count()) % 2


def syndrome(pauli: int, generators: list[int], n: int) -> int:
    return sum(symplectic(pauli, g, n) << j for j, g in enumerate(generators))


def span(generators: list[int]) -> list[int]:
    result = [0]
    for generator in generators:
        result += [element ^ generator for element in result]
    return result


def echelon(vectors: list[int]) -> dict[int, int]:
    basis: dict[int, int] = {}
    for vector in vectors:
        while vector:
            pivot = vector.bit_length() - 1
            if pivot in basis:
                vector ^= basis[pivot]
            else:
                basis[pivot] = vector
                break
    return basis


def rank(vectors: list[int]) -> int:
    return len(echelon(vectors))


def in_span(vector: int, generators: list[int]) -> bool:
    basis = echelon(generators)
    while vector:
        pivot = vector.bit_length() - 1
        if pivot not in basis:
            return False
        vector ^= basis[pivot]
    return True


def span_coordinates(target: int, vectors: list[int]) -> int:
    """Return a bit mask of vectors whose XOR is target, or raise ValueError."""
    basis: dict[int, tuple[int, int]] = {}
    for j, vector in enumerate(vectors):
        coefficient = 1 << j
        while vector:
            pivot = vector.bit_length() - 1
            if pivot in basis:
                old_vector, old_coefficient = basis[pivot]
                vector ^= old_vector
                coefficient ^= old_coefficient
            else:
                basis[pivot] = (vector, coefficient)
                break
    coefficient = 0
    while target:
        pivot = target.bit_length() - 1
        if pivot not in basis:
            raise ValueError("Target is not in the span")
        vector, selected = basis[pivot]
        target ^= vector
        coefficient ^= selected
    return coefficient


def nullspace(rows: list[int], width: int) -> list[int]:
    reduced = rows.copy()
    pivots = []
    at = 0
    for j in range(width):
        found = next((i for i in range(at, len(reduced)) if (reduced[i] >> j) & 1), None)
        if found is None:
            continue
        reduced[at], reduced[found] = reduced[found], reduced[at]
        for i in range(len(reduced)):
            if i != at and (reduced[i] >> j) & 1:
                reduced[i] ^= reduced[at]
        pivots.append(j)
        at += 1
    answer = []
    for free in range(width):
        if free in pivots:
            continue
        vector = 1 << free
        for pivot, row in zip(pivots, reduced):
            if (row >> free) & 1:
                vector |= 1 << pivot
        answer.append(vector)
    return answer


def normalizer_basis(generators: list[int], n: int) -> list[int]:
    mask = (1 << n) - 1
    rows = [(g >> n) | ((g & mask) << n) for g in generators]
    return nullspace(rows, 2 * n)


def embed(pauli: int, source_n: int, target_n: int, offset: int) -> int:
    mask = (1 << source_n) - 1
    return ((pauli & mask) << offset) | ((pauli >> source_n) << (target_n + offset))


def cnot(pauli: int, control: int, target: int, n: int) -> int:
    """Phase-free conjugation by CNOT(control, target)."""
    if (pauli >> control) & 1:
        pauli ^= 1 << target
    if (pauli >> (n + target)) & 1:
        pauli ^= 1 << (n + control)
    return pauli


def product_five_generators() -> list[int]:
    return [from_string(g + "IIIII") for g in FIVE_GENERATORS] + [
        from_string("IIIII" + g) for g in FIVE_GENERATORS
    ]


def local_factors(pauli: int, n: int) -> list[int]:
    return [from_string("I" * j + char + "I" * (n - j - 1))
            for j, char in enumerate(to_string(pauli, n)) if char != "I"]
