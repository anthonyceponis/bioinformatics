class Tree:
    def __init__(
        self,
        x: float,
        left=None,
        right=None,
        label="",
    ):
        self.x = x
        self.label = label
        self.left = left
        self.right = right

        # store references to leaves for quick access
        self.leaves = [] if left is None else left.leaves
        if right is not None:
            self.leaves += right.leaves

        # current node is a leaf!
        if len(self.leaves) == 0:
            self.leaves = [x]


def print_tree(node, indent="", is_left=True):
    if node is None:
        return

    # Process the right child first (displayed on top)
    if node.right:
        print_tree(node.right, indent + ("    " if is_left else "│   "), False)

    # Format the current node label/height
    # If it's a leaf, show the label. If internal, show the height (x).
    display_val = f"[{node.label}]" if node.label else f"({node.x:.2f})"

    print(indent + ("└── " if is_left else "┌── ") + display_val)

    # Process the left child (displayed on bottom)
    if node.left:
        print_tree(node.left, indent + ("│   " if is_left else "    "), True)


def average_distance(
    c1: Tree, c2: Tree, original_distance_matrix: list[list[int]]
) -> float:
    total_dist = 0

    for x in c1.leaves:
        for y in c2.leaves:
            total_dist += original_distance_matrix[x][y]

    return total_dist / (len(c1.leaves) * len(c2.leaves))


def upgma(clusters: list[Tree], original_distance_matrix: list[list[int]]):

    n = len(clusters)

    if n == 1:
        return clusters

    # Construct the cluster distance matrix
    distance_matrix = [
        [
            average_distance(clusters[i], clusters[j], original_distance_matrix)
            for i in range(n)
        ]
        for j in range(n)
    ]

    # Take clusters and cluster distance matrix, find two closest clusters, and merge.
    smallest_distance = float("inf")
    c1 = 0
    c2 = 0
    for i in range(n):
        for j in range(n):
            if distance_matrix[i][j] < smallest_distance:
                smallest_distance = distance_matrix[i][j]
                c1 = i
                c2 = j

    new_cluster = Tree(smallest_distance, clusters[c1], clusters[c2])
    clusters = [clusters[i] for i in range(n) if i != c1 and i != c2]
    clusters.append(new_cluster)

    return upgma(clusters, original_distance_matrix)


def visualise_tree(tree: Tree):
    return None


if __name__ == "__main__":
    labels = ["i", "j", "k", "l"]
    clusters = [Tree(i, None, None, labels[i]) for i in range(len(labels))]
    original_distance_matrix = [
        [0, 3, 4, 3],  # Distances from i to i, j, k, l
        [3, 0, 4, 5],  # Distances from j to i, j, k, l
        [4, 4, 0, 2],  # Distances from k to i, j, k, l
        [3, 5, 2, 0],  # Distances from l to i, j, k, l
    ]
    tree = upgma(clusters, original_distance_matrix)
    print_tree(tree)
