from bioservices.kegg import KEGG
import networkx as nx


class KEGGPathways:
    """
    KEGG PATHWAY Database API
    """

    def __init__(self, organism="Homo sapiens"):
        self.database = KEGG()
        self.organism = self.get_organism_code(organism.lower())

    def search_by_gene(self, gene_name: str):
        """

        Args:
            gene_name: gene name (ex. 'BRCA2')

        Returns:
            Dictionary with ids of all pathways containing given gene as keys and their full names as values.

        """
        try:
            pathways = self.database.get_pathway_by_gene(gene_name, self.organism)
            return pathways if pathways else {}
        except AttributeError:
            return {}

    def get_pathway(self, pathway_id: str, self_loops: bool = False):
        """

        Args:
            pathway_id: KEGG pathway id (ex. 'hsa04110')
            self_loops: information about whether or not include self loops in returned graph

        Returns:
            `networkx.DiGraph` object: Directed graph depicting pathway, with a comma-separated string
            containing gene names as graph nodes and directed edges representing interactions between genes.
            Each edge has weight 'type', which is a list of interaction types between two nodes.

        """

        G = nx.DiGraph()
        try:
            pathway = self.database.parse_kgml_pathway(pathway_id)
        except TypeError:
            # incorrect pathway_id
            pathway = None

        if pathway:
            names = {}
            for entry in pathway['entries']:
                # only intra-pathway interactions taken into account
                if entry['gene_names']:
                    names[entry['id']] = {'name': entry['gene_names'], 'type': entry['type']}

            for rel in pathway['relations']:
                if rel['entry1'] in names.keys() and rel['entry2'] in names.keys():
                    e1 = names[rel['entry1']]['name']
                    e2 = names[rel['entry2']]['name']
                    G.add_node(e1, type=names[rel['entry1']]['type'])
                    G.add_node(e2, type=names[rel['entry2']]['type'])
                    if G.has_edge(e1, e2):
                        G[e1][e2]['type'] = G[e1][e2]['type'] + [rel['name']]
                    else:
                        # assumption of interaction direction entry1 -> entry2 #TODO: validate
                        if e1 != e2 or (e1 == e2 and self_loops):
                            G.add_edge(e1, e2, type=[rel['name']])

        not_gene_nodes = []
        for node in G.nodes():
            # only interactions between genes
            if G.node[node]['type'] != 'gene':
                for in_edge in G.in_edges(node):
                    for out_edge in G.out_edges(node):
                        if in_edge[0] != out_edge[1] or (in_edge[0] == out_edge[1] and self_loops):
                            G.add_edge(in_edge[0], out_edge[1], type=['indirect'])
                not_gene_nodes.append(node)
        G.remove_nodes_from(not_gene_nodes)

        return G

    def fetch_organism_codes(self):
        """

        Returns:
            Dictionary with organisms as keys, and KEGG organism codes as values
            {   'homo sapiens' : 'hsa',
                'human' : 'hsa',
                ...
            }

        """
        codes = {}
        for line in self.database.list('organism').split('\n'):
            if line:
                code = line.split('\t')[1]
                org = line.split('\t')[2]
                if '(' in org:
                    org = [x.strip().lower() for x in org[:-1].split('(')]
                    for o in org:
                        codes[o] = code
                else:
                    codes[org] = code
        return codes

    def get_organism_code(self, org: str):
        """

        Args:
            org: organism name (ex. 'Homo sapiens', 'human') - lowercase and uppercase optional

        Returns:
            str: KEGG organism code

        """
        codes = self.fetch_organism_codes()
        try:
            return codes[org]
        except KeyError:
            print('Invalid organism name.')
            raise

    def get_gene_code(self, gen: str):
        """

        Args:
            gen: gene name (ex. 'FGR', 'NIPAL1')

        Returns:
            KEGG gene code

        """
        code_gen = self.database.find(self.organism, gen)

        if code_gen == str('\n'):
            code_gen = str()
            print('Invalid gene name: '+str(gen))
        return code_gen
