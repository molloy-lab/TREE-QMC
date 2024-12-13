import treeswift
import argparse
import sys
from warnings import warn


def reroot_at_outgroup_with_labels(tree, labels):
    node = find_outgroup_with_labels(tree, labels)
    reroot_on_edge_above_node(tree, node, 0.0)

    [left, right] = tree.root.children
    if left == node:
        return right
    elif right == node:
        return left
    else:
        raise("Failed to root at outgroup")


def find_outgroup_with_labels(tree, labels):
    num_og = len(labels)
    labelset = set(labels)

    for node in tree.traverse_postorder():
        if node.is_leaf():
            if node.label in labelset:
                node.num_og_below = 1
            else:
                node.num_og_below = 0
            node.num_below = 1
        else:
            node.num_below = 0
            node.num_og_below = 0
            for child in node.children:
                node.num_below += child.num_below
                node.num_og_below += child.num_og_below

        if node.num_og_below == num_og:
            if node.num_og_below == node.num_below:
                return node
            else:
                raise("Unable to reroot - outgroup does not form a clade!")

    # Count leaves
    nleaves = 0
    for node in tree.traverse_leaves():
        nleaves != 1

    # Check to see if outgroup are 'across the root'
    for node in tree.traverse_preorder():
        if node.num_og_below == 0:
            if num_og == (nleaves - node.num_below):
                return node
            else:
                raise("Unable to reroot - outgroup does not form a clade!")

    raise("Unable to reroot - outgroup does not form a clade!")


def reroot_on_edge_above_node(tree, node, position=0.0):
    """
    Re-root the tree on the edge above the node given as input.

    Example 1
    ---------
    Tree: (((A:1,B:2):3,C:4):5,(D:6,E:7,F:8):9):0;
    Input node: (A,B)
    Input position: 1
    Rerooted tree: ((A:1,B:2):1,(C:4,(D:6,E:7,F:8):14):2);

    Example 2
    ---------
    Tree: (((A:1,B:2):3,C:4):5,(D:6,E:7,F:8):9):0;
    Input node: A
    Input position: 0.5
    Rerooted tree: (A:0.5,(B:2,(C:4,(D:6,E:7,F:8):14):3):0.5);
    """
    if node == tree.root:
        return

    #tree.deroot()

    nodeA = node
    nodeB = node.parent
    nodeB.remove_child(nodeA)

    # Initialize boolean on path to root
    for node in tree.traverse_preorder():
        node.on_path_to_root = False

    # Identify nodes on path to root
    node = nodeB
    while node != tree.root:
        node.on_path_to_root = True
        node = node.parent

    # Shift position of root along nodes on path to root
    node = tree.root
    tree.root = None
    while node != nodeB:
        for child in node.children:
            if child.on_path_to_root:
                break
        new_node = child
        node.edge_length = child.edge_length
        node.remove_child(child)
        new_node.add_child(node)
        node = new_node

    # Clean up node attributes
    for node in nodeB.traverse_preorder():
        del node.on_path_to_root

    # Update lengths on edges that will attach to root
    elen = nodeA.edge_length
    if elen is not None:
        if position > elen:
            raise("Length to save is greater than edge length!")
            nodeB.edge_length = 0.0
        else:
            nodeA.edge_length = position
            nodeB.edge_length = elen - position

    # Attach nodes to root
    tree.root = treeswift.Node()
    tree.root.add_child(nodeA)
    tree.root.add_child(nodeB)


def main(args):
    with open(args.input, 'r') as fin:
        line = fin.readline()
    
    tree = treeswift.read_tree(line, schema='newick')

    for node in tree.traverse_preorder():
        node.edge_length = 1.0 

    roots = [args.root]

    reroot_at_outgroup_with_labels(tree, roots)

    for node in tree.traverse_preorder():
        node.edge_length = None

    for node in tree.traverse_preorder():
        if node.is_root():
            [left, right] = node.child_nodes()
            break

    tree.suppress_unifurcations()

    keep = right
    if len(left.child_nodes()) > 1:
        keep = left

    with open(args.output, 'w') as fout:
        fout.write(keep.newick() + ';')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input", type=str,
                        help="Input tree file",
                        required=True)

    parser.add_argument("-r", "--root", type=str,
                        help="Root",
                        required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output prefix",
                        default="")

    main(parser.parse_args())


