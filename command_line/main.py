import argparse

from methods import Method
from models import Phenotype, Experiment

from .types import Slice, one_of, Indices
from .parser import Parser, Argument


class MethodParser(Parser):

    def __init__(self, method, **kwargs):
        # TODO: add arguments from method to dir(self) and use dir insetead fo var in Parser.init
        super().__init__(**kwargs)
        self.method = method


class PhenotypeFactory(Parser):
    """Provide {name} samples"""

    # at least one file is always required
    # TODO do not require name
    files = Argument(type=argparse.FileType('r'), nargs='+')

    name = Argument(help='Your custom name for this set of samples.')
    # TODO support multiple files:
    samples = Argument(
        #action='append',
        nargs='*',
        help='Names of samples (columns) to be extracted from the file. '
             'Sample names are determined from the first non-empty row. '
             'Samples for each of files should be separated by space.'
    )

    # we want to handle either ":4", "5:" or even "1,2,3"
    columns = Argument(
        type=one_of(Slice, Indices),
        # user may (but do not have to) specify columns
        # to be extracted from given file.
        nargs='*',
        help='Columns to be extracted from files: '
             'either a comma delimited list of 0-based numbers (e.g. 0,2,3) '
             'or a range defined using Python slice notation (e.g. 3:10). '
             'Columns for each of files should be separated by space.'
        # TODO assert if columns: len(files) == len(columns)
    )

    def produce(self, unknown_args=None):
        opts = self.namespace
        name = opts.name or self.name

        if opts.files:
            # load all files
            sample_collections = []
            for i, file_obj in enumerate(opts.files):
                sample_collections.append(
                    Phenotype.from_file(
                        f'Sample collection, part {i} of {name}',
                        file_obj,
                        columns_selector=opts.columns.get_iterator if opts.columns else None,
                        samples=opts.samples
                    )
                )

            opts.phenotype = sum(sample_collections, Phenotype(name))
        return opts


class SingleFileExperimentFactory(Parser):
    """Provide both case and control samples within a single file."""

    # exactly one file is required
    files = Argument(
        type=argparse.FileType('r'),
        nargs=1,
        # required=True,
        help='file with samples for both control and cases.'
    )
    case = Argument(
        type=one_of(Slice, Indices),
        # required=True,
        help='columns from which case samples should be extracted.'
    )
    control = Argument(
        type=one_of(Slice, Indices),
        # required=True,
        help='columns from which control samples should be extracted.',
    )

    def produce(self, unknown_args=None):

        opts = self.namespace

        if opts.files:
            opts.control = PhenotypeFactory(name='control', files=opts.files, columns=opts.control).produce()
            # reuse the same file
            for f in opts.files:
                f.seek(0)
            opts.case = PhenotypeFactory(name='case', files=opts.files, columns=opts.case).produce()

        return opts


class CLIExperiment(Parser):
    """Use both: case and control or data to create an Experiment."""

    pull_to_namespace_above = True

    control = PhenotypeFactory()
    case = PhenotypeFactory()
    data = SingleFileExperimentFactory()

    def produce(self, unknown_args=None):

        opts = self.namespace
        if opts.data.files:
            if opts.control.files or opts.case.files:
                raise ValueError('Cannot handle data and case/control at once')
            #self.data.produce()

            opts.case = self.data.namespace.case
            opts.control = self.data.namespace.control
        else:
            del opts.data

        opts.experiment = Experiment(opts.case.phenotype, opts.control.phenotype)

        return opts


class CLI(Parser):
    """The main parser, the one exposed directly to the user."""

    methods = {
        method.name: method
        for method in Method.members
    }
    method_name = Argument(choices=methods, required=True, name='method', short='m')
    experiment = CLIExperiment()

    def produce(self, unknown_args):
        options = self.namespace

        # first - take an appropriate method class
        method = self.methods[options.method]

        # initialize parser for this method
        # (different methods require different arguments)
        method_parser = MethodParser(method=method)

        # parse arguments
        method_options, unknown_args = method_parser.parse_known_args(unknown_args)

        # ant initialize the method with these arguments
        options.method = method(**vars(method_options))

        return options

    def parse(self, args):

        if '-h' in args or '--help' in args:
            self.parser.parse_args(args)

        options, unknown_args = self.parse_known_args(args)
        assert not unknown_args

        return options

