from bioservices.kegg import KEGG
import networkx as nx


class KEGGPathways:

    def __init__(self, organism="Homo sapiens"):
        self.database = KEGG()
        self.organism = self.get_organism_code(organism)

    def search_by_gene(self, gene_name):
        """

        Args:
            gene_name: gene name (ex 'BRCA2')

        Returns:
            Dictionary with all ids of all pathways cointaining given gene as keys and their full names as values.

        """
        try:
            pathways = self.database.get_pathway_by_gene(gene_name, self.organism)
            return pathways if pathways else {}
        except AttributeError:
            return {}

    def search_by_pathway(self, pathway_id):
        """

        Args:
            pathway_id: KEGG pathway id (ex. 'hsa04110')

        Returns:
            G: Directed graph (networkx.DiGraph object) depicting pathway, with a comma-separated string
            containing gene names as graph nodes and directed edges representing interactions between genes.
            Each edge has weight 'type', which is a list of interaction types between two nodes.

        """
        G = nx.DiGraph()
        pathway = self.database.parse_kgml_pathway(pathway_id)
        names = {x['id']: x['name'] for x in pathway['entries']}
        for x in pathway['relations']:
            e1 = names[x['entry1']]
            e2 = names[x['entry2']]
            if e1 != 'undefined' and e2 != 'undefined':
                # assumption of interaction direction entry1 -> entry2 #TODO: validate
                record1 = self.database.get(e1).split()
                name1 = ''.join(record1[record1.index('NAME') + 1: record1.index('DEFINITION')])
                record2 = self.database.get(e2).split()
                name2 = ''.join(record2[record2.index('NAME') + 1: record2.index('DEFINITION')])
                G.add_nodes_from([name1, name2])
                if G.has_edge(name1, name2):
                    G[name1][name2]['type'] = G[name1][name2]['type'] + [x['name']]
                else:
                    G.add_edge(name1, name2, type=[x['name']])
        return G

    def fetch_organism_codes(self):
        """

        Returns:
            Dictionary with organisms as keys, and KEGG organism codes as values:
            {   'Homo sapiens' : 'hsa',
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
                    org = [x.strip() for x in org[:-1].split('(')]
                    for o in org:
                        codes[o] = code
                else:
                    codes[org] = code
        return codes

    def get_organism_code(self, org):
        """

        Args:
            org: organism name (ex. 'Homo sapiens', 'human')

        Returns:
            KEGG organism code

        """
        codes = self.fetch_organism_codes()
        try:
            return codes[org]
        except KeyError:
            print('Invalid organism name.')
            raise
