import sys
import dendropy

def printSubtree(v):
	if v in children:
		print v, children[v][0], children[v][1]
		for c in children[v]:
			printSubtree(c)

with open(sys.argv[1], "r") as ins:
	print len([None for line in ins])

with open(sys.argv[1], "r") as ins:
	for line in ins:
		children = {}
		parent = {}
		name = {}

		tree = dendropy.Tree.get_from_string(line, schema="newick", rooting="force-rooted", suppress_internal_node_taxa=True, suppress_leaf_node_taxa=True)
		if len(tree.seed_node.child_edges()) > 2:
			tree.reroot_at_edge(tree.seed_node.child_edges()[0])

		for p in tree.internal_nodes():
			children[hash(p)] = []
			for c in p.child_node_iter():
				children[hash(p)].append(hash(c))
				parent[hash(c)] = hash(p)

		for v in tree.leaf_nodes():
			#name[hash(v)] = v.label.replace(' ', '_')
			name[hash(v)] = v.label.split(' ')[0]
		print len(name)
		print len(children)
		print hash(tree.seed_node)
		for n in name:
			print name[n], n
		printSubtree(hash(tree.seed_node))
		#print(tree.as_ascii_plot())