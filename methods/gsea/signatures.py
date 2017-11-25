import gzip
from pathlib import Path
from urllib.request import urlretrieve

import os

from command_line.parser import Argument, Parser, action
from models import Gene

REMOTE = 'https://github.com/kn-bibs/pathways-data/raw/master/gsea/msigdb/'
DATA_DIR = Path('data')


def gzip_open_text(path, mode='r'):
    return gzip.open(path, mode + 't')


class GeneSet:

    def __init__(self, name, genes, url=None):
        self.name = name
        self.genes = {Gene(name) for name in genes}
        self.genes = genes
        self.url = url
        self.enrichment = None

    def __contains__(self, item):
        # TODO test
        return item in self.genes

    def __len__(self):
        return len(self.genes)

    def __lt__(self, other):
        return self.enrichment < other.enrichment

    def __repr__(self):
        return f'<GeneSet: {self.name} with {len(self.genes)}>'


class MolecularSignatureDatabase:

    def __init__(self, path, label=None):
        self.path = DATA_DIR / path
        self.label = label
        self.gene_sets = {}
        self.load()

    def load(self, opener=open):
        if self.path.suffix.endswith('.gz'):
            opener = gzip_open_text

        with opener(self.path) as f:
            for line in f:
                name, url, *genes = line.strip().split('\t')
                self.gene_sets[name] = GeneSet(name, genes, url)


class RemoteDatabase(MolecularSignatureDatabase):

    def __init__(self, set_name, identifiers='symbols', version=6.1, remote=REMOTE, label=None):
        path = f'{version}/{set_name}.v{version}.{identifiers}.gmt.gz'
        self.raw_path = path
        self.remote = remote
        super().__init__(path, label)

    def fetch(self):
        url = self.remote + str(self.raw_path)
        os.makedirs(self.path.parent, exist_ok=True)
        urlretrieve(url, self.path)

    def load(self):
        if not self.path.exists():
            self.fetch()
        super().load()


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
            n.database = MolecularSignatureDatabase(name)
        return n
