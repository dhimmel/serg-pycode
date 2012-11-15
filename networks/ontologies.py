from operator import attrgetter

class AcyclicDirectedGraph(object):
    
    def __init__(self):
        """An Acyclic Directed Graph object most likely used for modelling
        ontologies. Specific ontologies should implement a subclass
        """
        self.nodes = set()
        self.root_nodes = set()

    def create_relationships(self, parents, children):
        """Make all children children of all parents and vice versa.
        parents and childrens are iterables of nodes."""
        for parent in parents:
            for child in children:
                parent.children.add(child)
                child.parents.add(parent)

    def snip_node(self, node):
        """Snip out a node from the ontology. Children of the snipped
        node become children of the snipped node's parents."""
        for child in node.children:
            child.parents.remove(node)
            child.parents |= node.parents
            if not child.parents:
                self.root_nodes.add(child)
        for parent in node.parents:
            parent.children.remove(node)
            parent.children |= node.children
        self.nodes.remove(node)
        if node in self.root_nodes:
            self.root_nodes.remove(node)
        del node
    
    def remove_node(self, node, delete_disconnected_children=False):
        """Remove a node from the ontology including all relationships
        to that node."""
        for child in node.children:
            child.parents.remove(node)
            if delete_disconnected_children and not child.parents:
                self.remove_node(child, delete_disconnected_children=True)
        for parent in node.parents:
            parent.children.remove(node)
        self.nodes.remove(node)
        if node in self.root_nodes:
            self.root_nodes.remove(node)
        del node

    
    def get_all_descendants(self, node, include_self=False):
        descendants = set()
        incomplete_nodes = set([node]) if include_self else node.children
        while incomplete_nodes:
            node = incomplete_nodes.pop()
            if node in descendants:
                continue
            descendants.add(node)
            incomplete_nodes |= node.children
        return descendants
    
    def descendants_of_depth(self, node, depth, include_intermediaries=False):
        descendants = set()
        current_level = set([node])
        for i in range(depth):
            next_level = set()
            for node in current_level:
                next_level |= node.children
            if include_intermediaries:
                descendants |= current_level
            current_level = next_level
        descendants |= current_level
        return descendants
        
    
    def keep_only(self, keeper_nodes, keeper_node_format='code'):
        """Snips all nodes not specified in keeper_nodes.
        keeper_node_format indicates the attribute name from a Node object
        which an element from keeper_nodes must match for node preservation.
        """
        keeper_nodes = set(keeper_nodes)
        for node in list(self.nodes):
            if attrgetter(keeper_node_format)(node) not in keeper_nodes:
                self.snip_node(node)
        self.compute_dicts()
    
    def add_annotation_type(self, type_):
        for node in self.nodes:
            node.direct_annots[type_] = set()
            node.prop_annots[type_] = set()

    def __repr__(self):
        s = ("Number of Nodes: " + str(len(self.nodes)) + "\n" +
             "Number of Root Nodes: " + str(len(self.root_nodes)))
        return s
    
    def __iter__(self):
        return iter(self.nodes)

    def search_by_attr(self, to_match, attr, attr_modifier, plural_attr, matching_fxn, cutoff=None):
        """Finds node from ontology that best matches. 
        Minimizes score of matching function."""
        bests = []
        best_score = None
        for node in self.nodes:
            potentials = attrgetter(attr)(node)
            if not plural_attr:
                potentials = [potentials]
            for potential in potentials:
                potential = attr_modifier(potential)
                score = matching_fxn(to_match, potential)
                current = (score, node)
                if not bests:
                    bests.append(current)
                    best_score = score
                elif score < best_score:
                    bests = [current]
                    best_score = score
                elif score == best_score:
                    if not any(best[1] is node for best in bests):
                        bests.append(current)
        if cutoff:
            bests = filter(lambda x: x[0] <= cutoff, bests)
        return bests
    
    def write_structure(self, path):
        self.indent = 0
        with open(path, 'w') as f:
            for node in self.root_nodes:
                self.write_single_node_structure(node, f)
                    
    def write_single_node_structure(self, node, f):
        s = '  ' * self.indent + '--' + str(node) + '\n'
        f.write(s)
        for child in node.children:
            self.indent += 1
            self.write_single_node_structure(child, f)
            self.indent -= 1

   
class AcyclicDirectedGraphNode(object):
    
    def __init__(self):
        """All Acyclic Directed Graph nodes must have
        a parrents and children attributes.
        """
        self.parents = set()
        self.children = set()
        self.direct_annots = dict()
        self.prop_annots = dict()

    def propogate_annots(self, type_):
        """Returns annotation of the specified type_ to itself 
        and all of its descendants. Propogates upward.
        """
        if not self.prop_annots[type_]:
            self.prop_annots[type_] |= self.direct_annots[type_]
            for child in self.children:
                self.prop_annots[type_] |= child.propogate_annots(type_)
        return self.prop_annots[type_]
