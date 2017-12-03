import gzip
from pathlib import Path
from typing import Mapping, Sequence
from urllib.request import urlretrieve

import os

from declarative_parser.parser import Argument, Parser, action
from models import Gene
from utils import jit

REMOTE = 'https://github.com/kn-bibs/pathways-data/raw/master/gsea/msigdb/'
DATA_DIR = Path('data')


def gzip_open_text(path, mode='r'):
    return gzip.open(path, mode + 't')


class GeneSet:

    def __init__(self, name, genes, url=None):
        self.name = name
        self.genes = {Gene(name) for name in genes}
        self.url = url
        self.enrichment = None

        # Following is meant to enable fast __contains__ test:
        # granted, one might want to simply use memory address
        # retrieved with id(gene) but it is not guaranteed to
        # (and will not) work with multiprocessing.
        # On the other hand a set of integers is easily pick-able,
        # as is the integer-holding 'id' attribute of Gene(s).
        # Trivia: set(list-comprehension) is faster than set(generator).
        self.gene_ids = set([gene.id for gene in self.genes])

    def restrict_to_genes(self, genes: Sequence[Gene]):
        """Clear itself of genes which are not in the the provided list.

        Returns: set of genes removed from the dataset
        """
        excessive_genes = self.genes.difference(genes)

        for gene in excessive_genes:
            self.genes.discard(gene)
            self.gene_ids.discard(gene.id)

        return excessive_genes

    def __contains__(self, gene: Gene):
        return gene.id in self.gene_ids

    def __len__(self):
        return len(self.genes)

    def __lt__(self, other):
        return self.enrichment < other.enrichment

    def __repr__(self):
        return f'<GeneSet: {self.name} with {len(self.genes)}>'


class MolecularSignatureDatabase:

    def __init__(self, gene_sets: Mapping[str, GeneSet], label=None):
        self.label = label
        self.gene_sets = gene_sets


class GMTSignatureDatabase(MolecularSignatureDatabase):

    def __init__(self, path, label=None):
        self.path = DATA_DIR / path
        gene_sets = self.load()
        super().__init__(gene_sets, label)

    def load(self, opener=open):

        if self.path.suffix.endswith('.gz'):
            opener = gzip_open_text

        gene_sets = {}

        with opener(self.path) as f:
            for line in f:
                name, url, *genes = line.strip().split('\t')
                gene_sets[name] = GeneSet(name, genes, url)

        return gene_sets


class RemoteDatabase(GMTSignatureDatabase):

    def __init__(self, set_name, identifiers='symbols', version=6.1, remote=REMOTE, label=None):
        path = f'{version}/{set_name}.v{version}.{identifiers}.gmt.gz'
        self.raw_path = path
        self.remote = remote
        GMTSignatureDatabase.__init__(self, path, label)

    def fetch(self):
        url = self.remote + str(self.raw_path)
        os.makedirs(self.path.parent, exist_ok=True)
        urlretrieve(url, self.path)

    def load(self, opener=open):
        if not self.path.exists():
            self.fetch()
        return super().load(opener)


DATABASE_PRESETS = {
    'H': (
        'h.all',
        'hallmark gene sets'
    )
}


class DatabaseParser(Parser):
    """Help"""

    __skip_if_absent__ = False

    name_or_path = Argument(
        default='H',
        help='Name of Molecular Signature Database to use. '
             'By default hallmark genes will be used. '
             f'Use one of following: {", ".join(DATABASE_PRESETS)} '
             'or provide a path to custom database (.gmt file). '
             # http://software.broadinstitute.org/gsea/msigdb_license_terms.jsp
             'Please note that due to licencing restrictions data '
             'from following sources: '
             'KEGG, BioCarta and AAAS/STKE Cell Signaling Database '
             'are currently unavailable.'
        ,
    )
    version = Argument(
        default=6.1
    )
    identifiers = Argument(
        choices=['symbols', 'entrez'],
        default='symbols'
    )
    remote = Argument(
        # we may want to create a mirror when things go serious to
        # do not overuse BI and to gain independence
        default=REMOTE
    )

    @action
    def show_gene_sets(namespace):
        try:
            parser = DatabaseParser(**vars(namespace))
            db = parser.produce(None).database
        except Exception as e:
            print(e)
            raise
        for name, gene_set in db.gene_sets.items():
            print(name, gene_set.url)

    def produce(self, unknown_args):
        n = self.namespace
        name = n.name_or_path

        # try to get one of pre-defined database by name
        if name in DATABASE_PRESETS:
            set_name, label = DATABASE_PRESETS[name]

            n.database = RemoteDatabase(
                set_name, label=label,
                version=n.version,
                identifiers=n.identifiers,
                remote=n.remote
            )
        else:
            # or to load a database from path
            n.database = GMTSignatureDatabase(name)
        return n
