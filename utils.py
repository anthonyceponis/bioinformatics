def print_dict_grid(grid_dict, default=".", padding=2):
    # 1. Extract all row and column indices
    rows = [r for r, c in grid_dict.keys()]
    cols = [c for r, c in grid_dict.keys()]

    # 2. Handle empty dictionary case
    if not rows:
        print("Grid is empty.")
        return

    # 3. Define boundaries
    min_r, max_r = min(rows), max(rows)
    min_c, max_c = min(cols), max(cols)

    # 4. Print the grid
    for r in range(min_r, max_r + 1):
        row_str = []
        for c in range(min_c, max_c + 1):
            # Get value or default, then convert to string
            val = str(grid_dict.get((r, c), default))
            # Center the value for alignment
            row_str.append(val.center(padding))

        print(" ".join(row_str))
