import math
from collections import defaultdict


def construct_lcs_grid(
    u: str, v: str, first_col: list[int], first_row: list[int]
) -> defaultdict:
    """Constructs the lcs grid given strings u,v, and the initial offsets for the first row and col, to be used for the four russians technique."""

    m = len(u)
    n = len(v)

    dp = defaultdict(int)

    # init
    for i in range(m + 1):
        dp[(i, 0)] = first_col[i]
    for j in range(n + 1):
        dp[(0, j)] = first_row[j]

    # fill table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if u[i - 1] == v[j - 1]:
                dp[(i, j)] = max(dp[(i, j)], dp[(i - 1, j - 1)] + 1)
            else:
                dp[(i, j)] = max(dp[(i - 1, j)], dp[(i, j - 1)])

    return dp


def offset_vectors_generator(t, current_vector=[]):
    """
    Generates all possible offset vectors of length t, where each element is 0 or 1.
    """
    if len(current_vector) == t:
        yield tuple(current_vector)
        return

    for offset in range(2):
        yield from offset_vectors_generator(t, current_vector + [offset])


def precompute_grids(t: int) -> defaultdict:
    """Precomputes all possible t x t grids, of which there are 4^t x 4^t."""

    grids = defaultdict(defaultdict)
    bases = ("A", "T", "C", "G", "$")

    # Generate all possible strings of length t
    def generate_strings(length, current_string=""):
        if length == 0:
            yield current_string
            return
        for base in bases:
            yield from generate_strings(length - 1, current_string + base)

    strings_t = list(generate_strings(t))

    for u in strings_t:
        for v in strings_t:
            for first_row_offsets in offset_vectors_generator(t):
                for first_col_offsets in offset_vectors_generator(t):
                    first_row = [sum(first_row_offsets[:i]) for i in range(t + 1)]
                    first_col = [sum(first_col_offsets[:i]) for i in range(t + 1)]
                    grids[(u, v, first_row_offsets, first_col_offsets)] = (
                        construct_lcs_grid(u, v, first_col, first_row)
                    )

    return grids


def four_russians(u: str, v: str) -> int:
    """
    Computes LCS between strings u,v using four russians algorithm in O(n^2/log(n)).
    String lengths must be equal and powers of 2 for implementation convinience.
    Strings must also only contain letters ATCG (4 bases).
    """

    m = len(u)
    n = len(v)

    # Pad strings to nearest power of 2
    next_power_of_2 = 1
    while next_power_of_2 < max(m, n):
        next_power_of_2 *= 2

    u = u.ljust(next_power_of_2, "$")
    v = v.ljust(next_power_of_2, "$")
    n = next_power_of_2

    t = int(math.log(n, 4)) if n > 0 else 0
    if t == 0:
        # Handle small n case separately
        return construct_lcs_grid(u, v, [0] * (n + 1), [0] * (n + 1))[(m, n)]
    grids = precompute_grids(t)

    F = defaultdict(int)  # {(i,j): score} main dp grid.

    # Initialize first row and col of F
    for i in range(n + 1):
        F[(i, 0)] = 0
        F[(0, i)] = 0

    assert n / t == n // t, "number of subgrids must fit exactly into the overall grid."
    subgrids = n // t

    for i in range(subgrids):
        for j in range(subgrids):
            u_chunk = u[i * t : (i + 1) * t]
            v_chunk = v[j * t : (j + 1) * t]

            # These are the offsets for the first row and col of the subgrid
            first_row_offsets = tuple(
                F[(i * t, j * t + k + 1)] - F[(i * t, j * t + k)] for k in range(t)
            )
            first_col_offsets = tuple(
                F[(i * t + k + 1, j * t)] - F[(i * t + k, j * t)] for k in range(t)
            )

            # The precomputed grid gives us the LCS scores for the subgrid, with the given offsets.
            lcs_grid = grids[(u_chunk, v_chunk, first_row_offsets, first_col_offsets)]

            # The top-left corner of the subgrid in the main F grid.
            offset = F[(i * t, j * t)]

            # Update the main F grid with the last row and col from the precomputed grid.
            for k in range(t + 1):
                F[(i * t + k, (j + 1) * t)] = offset + lcs_grid.get((k, t), 0)
                F[((i + 1) * t, j * t + k)] = offset + lcs_grid.get((t, k), 0)

    return F[(m, n)]


if __name__ == "__main__":
    tests = [
        ("GACGTAGCATAAGCGC", "TGCAACGTATAACGGG", 11),
        ("TCAGTACTAGTTATCAGTCTAGTCAGCTACTA", "GTCAGTACTAGTTATCAGTCTAGTCAGCTACT", 31),
    ]

    for s1, s2, expected_score in tests:
        # Basic LCS
        n = len(s1)
        first_row = [0 for _ in range(n + 1)]
        first_col = [0 for _ in range(n + 1)]
        score = construct_lcs_grid(s1, s2, first_row, first_col)[(n, n)]
        print("BASIC LCS ALGORITHM")
        print("################")
        print(f"s1 = {s1}")
        print(f"s2 = {s2}")
        print(f"expected score = {expected_score}")
        print(f"score = {score}")

        # Print spacing

        print()
        print("################")
        print("################")
        print("################")
        print()

        # Four russians
        score = four_russians(s1, s2)
        print("FOUR RUSSIANS ALGORITHM")
        print("################")
        print(f"s1 = {s1}")
        print(f"s2 = {s2}")
        print(f"expected score = {expected_score}")
        print(f"score = {score}")

        print()
        print()
