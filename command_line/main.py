import argparse

from command_line.method_parser import MethodParser
from methods import Method
from models import Phenotype, Experiment
from .parser import Parser, Argument
from .types import Slice, one_of, Indices


class PhenotypeFactory(Parser):
    """Provide {name} samples"""

    # at least one file is always required
    # TODO do not require name
    files = Argument(type=argparse.FileType('r'), nargs='+')

    name = Argument(help='Your custom name for this set of samples.')

    samples = Argument(
        type=lambda x: x.split(','),
        nargs='*',
        help='Names of samples (columns) to be extracted from the file. '
             'Sample names are determined from the first non-empty row. '
             'Use comma to separate samples. '
             'Samples for each of files should be separated by space.'
    )

    # we want to handle either ":4", "5:" or even "1,2,3"
    columns = Argument(
        type=lambda x: [one_of(Slice, Indices)(y) for y in x.split(' ')],
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
                        columns_selector=opts.columns[i].get_iterator if opts.columns else None,
                        samples=opts.samples[i] if opts.samples else None
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

            opts.case = self.data.namespace.case
            opts.control = self.data.namespace.control
        elif opts.case.files and opts.control.files:
            # that's nice :)
            pass
        else:
            raise ValueError('Neither data nor (case & control) have been provided!')

        del opts.data

        opts.experiment = Experiment(opts.case.phenotype, opts.control.phenotype)

        return opts


class CLI(Parser):
    """The main parser, the one exposed directly to the user."""

    method_name = Argument(choices=Method.members, name='method', short='m', optional=False)
    experiment = CLIExperiment()

    @staticmethod
    def create_method(name):
        # first - take an appropriate method class
        method = Method.members[name]

        # initialize parser for this method
        # (different methods require different arguments)
        method_parser = MethodParser(method=method)

        return method_parser

    def parse(self, args):
        help_args = {'-h', '--help'}

        if help_args.intersection(args):
            args_without_help = [
                arg
                for arg in args
                if arg not in help_args
            ]

            # do we have method name specified?
            if len(args_without_help):
                method_name = args_without_help[0]
                method = Method.members[method_name]
                method_parser = MethodParser(method=method)

                return method_parser.parse(args_without_help[1:] + ['-h'])

        return super().parse(args)

    def produce(self, unknown_args):
        options = self.namespace

        method_parser = self.create_method(options.method)

        # parse arguments
        method_options, remaining_unknown_args = method_parser.parse_known_args(unknown_args)

        for argument in unknown_args[:]:
            if argument not in remaining_unknown_args:
                unknown_args.remove(argument)

        # ant initialize the method with these arguments
        options.method = method_parser.method(**vars(method_options))

        return options
