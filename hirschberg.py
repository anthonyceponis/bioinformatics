from collections import defaultdict
from needleman_wunsch import needleman_wunsch


def hirschberg(x: str, y: str) -> tuple[str, str]:
    """
    Compute global alignment using Hirschberg-style recursion.
    Simplified for LCS-style scoring:
        match = 1
        mismatch = 0
        gap = 0
    """

    # base case
    if len(y) <= 1 or len(x) <= 1:
        _, x, y = needleman_wunsch(x, y, match=1, mismatch=0, gap=0)
        return x, y

    # split on x
    x_m = len(x) // 2

    prefix_dp = defaultdict(lambda: defaultdict(int))
    suffix_dp = defaultdict(lambda: defaultdict(int))

    # fill prefix dp: x[:x_m] vs y
    for i in range(len(y)):
        for j in range(x_m):
            if y[i] == x[j]:
                prefix_dp[i + 1][j + 1] = prefix_dp[i][j] + 1
            else:
                prefix_dp[i + 1][j + 1] = max(
                    prefix_dp[i][j + 1],
                    prefix_dp[i + 1][j]
                )

    # fill suffix dp backwards: x[x_m:] vs suffixes of y
    for i in range(len(y) - 1, -1, -1):
        for j in range(len(x) - 1, x_m - 1, -1):
            if y[i] == x[j]:
                suffix_dp[i][j] = suffix_dp[i + 1][j + 1] + 1
            else:
                suffix_dp[i][j] = max(
                    suffix_dp[i + 1][j],
                    suffix_dp[i][j + 1]
                )

    # find middle node / best split of y
    middlenode_i = 0
    middlenode_score = -1

    for i in range(len(y) + 1):
        curr = prefix_dp[i][x_m] + suffix_dp[i][x_m]

        if curr > middlenode_score:
            middlenode_score = curr
            middlenode_i = i

    # recurse
    x1, y1 = hirschberg(x[:x_m], y[:middlenode_i])
    x2, y2 = hirschberg(x[x_m:], y[middlenode_i:])

    return x1 + x2, y1 + y2


if __name__ == "__main__":
    ax, ay = hirschberg("ACGT", "AGT")
    print(ax)
    print(ay)
