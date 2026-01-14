from collections import defaultdict
from utils import print_dict_grid


def pair_matches(x: str, y: str) -> bool:
    """Checks if given pair is RNA complimentary"""
    return tuple(sorted([x, y])) == ("A", "U") or tuple(sorted([x, y])) == ("C", "G")


def print_rna_secondary_structure(sequence, pairs):
    # start with all dots
    structure = ["."] * len(sequence)

    # fill in the brackets for each pair
    for i, j in pairs:
        # We use min/max to ensure i is the opening and j is the closing
        start, end = min(i, j), max(i, j)
        structure[start] = "("
        structure[end] = ")"

    print("".join(structure))


def traceback(i: int, j: int, F: defaultdict, s: str) -> list[tuple[int, int]]:
    if i >= j:
        return []  # Base case: no more pairs possible

    # 1. Check if j is unpaired
    if F[(i, j)] == F.get((i, j - 1), 0):
        return traceback(i, j - 1, F, s)

    # 2. Check if i is unpaired
    if F[(i, j)] == F.get((i + 1, j), 0):
        return traceback(i + 1, j, F, s)

    # 3. Check if i,j is a pair
    if pair_matches(s[i], s[j]):
        if F[(i, j)] == F.get((i + 1, j - 1), 0) + 1:
            return [(i, j)] + traceback(i + 1, j - 1, F, s)

    # 4. Check for bifurcation
    for k in range(i + 1, j):
        if F[(i, j)] == F.get((i, k), 0) + F.get((k + 1, j), 0):
            return traceback(i, k, F, s) + traceback(k + 1, j, F, s)

    return []


def nussinov(s: str):
    """
    Deduces the RNA secondary structure from an RNA sequence, using the nussinov algorithm.

    Algorithm works by maximising the number of complimentary base pairs in a range, because base pair bonding is an exothermic reaction, so the more bonds, the lower the resting energy state, which molecules in nature tend to do.
    Optimal structures constructed through four cases on optimal substructures, which is the bases of the dp recursion.

    Coord (i,j) corresponds to interval [i,j].
    Init with main diagonal with score zero because cannot bond to yourself, and also zero the diagonal below main, because (i,i-1) represents an empty range, which trivially cannot contain bonds.

    Fill in the top-right half, expanding from main diagonal.
    Expanding from main diagonal corresponds with starting with many small intervals, gradually increasing the interval length.
    Optimal score is the top right cell.

    Filling in grid is mostly O(n^2), but bifurication case makes it O(n^3).
    """

    F = defaultdict(int)  # {(i, j): score}
    n = len(s)

    # init
    for i in range(n):
        F[(i, i)] = 0
        F[(i + 1, i)] = 0

    # fill main grid
    for offset in range(n):
        for i in range(n - offset):
            j = i + offset
            bifurication_case = 0
            for k in range(i + 1, j):
                bifurication_case = max(bifurication_case, F[(i, k)] + F[((k + 1, j))])
            match_reward = 1 if pair_matches(s[i], s[j]) else -float("inf")
            F[(i, j)] = max(
                F[(i + 1, j)],
                F[(i, j - 1)],
                F[(i + 1, j - 1)] + match_reward,
                bifurication_case,
            )

    return F


if __name__ == "__main__":
    s = "AAUCUGUUACGCA"
    F = nussinov(s)
    n = len(s)
    pairs = traceback(0, n - 1, F, s)
    print_dict_grid(F)
    print(s)
    print_rna_secondary_structure(s, pairs)
