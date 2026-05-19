def needleman_wunsch(x, y, match=1, mismatch=-1, gap=-1) -> tuple[int, str, str]:
    m, n = len(x), len(y)

    # DP table
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Backpointer table
    back = [[None] * (n + 1) for _ in range(m + 1)]

    # Initialise first column: align x chars with gaps
    for i in range(1, m + 1):
        dp[i][0] = dp[i - 1][0] + gap
        back[i][0] = "up"

    # Initialise first row: align y chars with gaps
    for j in range(1, n + 1):
        dp[0][j] = dp[0][j - 1] + gap
        back[0][j] = "left"

    # Fill table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            diag_score = dp[i - 1][j - 1] + (match if x[i - 1] == y[j - 1] else mismatch)
            up_score = dp[i - 1][j] + gap      # x[i-1] aligned with gap
            left_score = dp[i][j - 1] + gap    # gap aligned with y[j-1]

            best = max(diag_score, up_score, left_score)
            dp[i][j] = best

            if best == diag_score:
                back[i][j] = "diag"
            elif best == up_score:
                back[i][j] = "up"
            else:
                back[i][j] = "left"

    # Backtrack
    aligned_x = []
    aligned_y = []

    i, j = m, n
    while i > 0 or j > 0:
        move = back[i][j]

        if move == "diag":
            aligned_x.append(x[i - 1])
            aligned_y.append(y[j - 1])
            i -= 1
            j -= 1

        elif move == "up":
            aligned_x.append(x[i - 1])
            aligned_y.append("-")
            i -= 1

        elif move == "left":
            aligned_x.append("-")
            aligned_y.append(y[j - 1])
            j -= 1

    aligned_x.reverse()
    aligned_y.reverse()

    return dp[m][n], "".join(aligned_x), "".join(aligned_y)
