import argparse

from methods import Method
from models import Phenotype, Experiment

from .parser import Parser, Argument
from .types import Slice, one_of, Indices, dsv
from .method_parser import MethodParser


class PhenotypeFactory(Parser):
    """Provide {name} samples. Requires a file (or files) with samples.

     The files should come in Delimiter Separate Values format
     (like .csv or .tsv). The default delimiter is a tab character.
     First column of each file should contain gene identifiers.

     To use only a subset of samples from files(s) specify column numbers
     (--columns) or sample names (--samples) of desired samples.
     """

    files = Argument(
        type=argparse.FileType('r'),
        # at least one file is always required
        nargs='+',
        optional=False
    )

    name = Argument(help='Your custom name for this set of samples.')

    samples = Argument(
        type=dsv(str),
        nargs='*',
        as_many_as=files,
        help='Names of samples (columns) to be extracted from the file. '
             'Sample names are determined from the first non-empty row. '
             'Use comma to separate samples. '
             'Samples for each of files should be separated by space.'
    )

    columns = Argument(
        # we want to handle either ":4", "5:" or even "1,2,3"
        type=one_of(Slice, Indices),
        # user may (but do not have to) specify columns
        # to be extracted from given file(s).
        nargs='*',
        as_many_as=files,
        help='Columns to be extracted from files: '
             'either a comma delimited list of 0-based numbers (e.g. 0,2,3) '
             'or a range defined using Python slice notation (e.g. 3:10). '
             'Columns for each of files should be separated by space.'
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
                        samples=opts.samples[i] if opts.samples else None,
                        reverse_selection=getattr(opts, 'reverse', False)
                    )
                )

            opts.phenotype = sum(sample_collections, Phenotype(name))
        return opts


class SingleFileExperimentFactory(Parser):
    """Provide both case and control samples within a single file."""

    # exactly one file is required
    files = Argument(
        type=argparse.FileType('r'),
        nargs=1,    # transforms result into a single-element list
        optional=False,
        help='file with samples for both control and cases.'
    )
    case = Argument(
        type=one_of(Slice, Indices),
        nargs=1,
        help='columns from which case samples should be extracted.'
    )
    control = Argument(
        type=one_of(Slice, Indices),
        nargs=1,
        help='columns from which control samples should be extracted.',
    )

    def produce(self, unknown_args=None):

        opts = self.namespace

        def produce_phenotype(created_group, other_group):
            reverse = hasattr(opts, 'reverse_' + created_group)
            get_columns_from = created_group
            if reverse:
                get_columns_from = other_group

            return PhenotypeFactory(
                name=created_group,
                files=opts.files,
                columns=getattr(opts, get_columns_from),
                reverse=reverse
            ).produce()

        if opts.files:
            if not (opts.case and opts.control):
                if opts.case:
                    opts.reverse_control = True
                elif opts.control:
                    opts.reverse_case = True
                else:
                    raise ValueError(
                        'Neither --case nor --control provided: '
                        'please specify which columns should be used as control '
                        'and which should be used as case.'
                    )

            phenotypes = {'control': produce_phenotype('control', 'case')}
            # reuse the same file(s)
            for f in opts.files:
                f.seek(0)
            phenotypes['case'] = produce_phenotype('case', 'control')

            for name, phenotype in phenotypes.items():
                setattr(opts, name, phenotype)

        return opts


class CLIExperiment(Parser):
    """Use both: case and control or data to create an Experiment."""

    pull_to_namespace_above = True

    control = PhenotypeFactory()
    case = PhenotypeFactory()
    data = SingleFileExperimentFactory()

    def produce(self, unknown_args=None):

        opts = self.namespace
        if opts.data:
            if opts.control or opts.case:
                raise ValueError('Cannot handle data and case/control at once')

            opts.case = self.data.namespace.case
            opts.control = self.data.namespace.control
        elif opts.case and opts.control:
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

            if len(args_without_help) != 0:

                name = args_without_help[0]

                # in case of a conflict, help for both (for a sub-parser
                # and for a method) should be displayed.

                methods = {
                    name: MethodParser(method=method)
                    for name, method in Method.members.items()
                }

                def match_parser(subparsers):
                    return subparsers.get(name, None)

                all_subparsers = [methods, self.subparsers, self.lifted_parsers]

                for parser in filter(bool, map(match_parser, all_subparsers)):
                    return parser.parse(args_without_help[1:] + ['-h'])

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

