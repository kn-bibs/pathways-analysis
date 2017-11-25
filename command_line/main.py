import argparse
from pathlib import Path

from methods import Method
from models import SampleCollection, Experiment

from declarative_parser import Parser, Argument
from declarative_parser.types import Slice, one_of, Indices, dsv, Range
from declarative_parser.constructor_parser import ConstructorParser


class SampleCollectionFactory(Parser):
    """Provide {parser_name} samples. Requires a file (or files) with samples.

     The files should come in Delimiter Separated Values format
     (like .csv or .tsv). The default delimiter is a tab character.
     The first column of each file should contain gene identifiers.

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
             'Use a comma to separate samples. '
             'Samples for each of files should be separated by space.'
    )

    columns = Argument(
        # we want to handle either ":4", "5:" or even "1,2,3"
        type=one_of(Slice, Indices, Range),
        # user may (but do not have to) specify columns
        # to be extracted from given file(s).
        nargs='*',
        as_many_as=files,
        help='Columns to be extracted from files: '
             'either a comma delimited list of 0-based numbers (e.g. 0,2,3) '
             'or a range defined using Python slice notation (e.g. 3:10). '
             'Columns for each of files should be separated by space.'
    )

    delimiter = Argument(
        default='\t',
        help='Delimiter of the provided file(s). Default: tabulation mark.'
    )

    header = Argument(
        nargs='*',
        type=one_of(int, str),
        as_many_as=files,
        default=lambda file_object: 0,
        help='Defines how the sample names should be created. '
             'Provide a number to specify which line should be used '
             'to extract names for samples. Please remember that '
             'empty lines will be skipped. If your file has no row '
             'with sample names, provide a string to be used as a '
             'prefix for naming consecutive samples. '
             'For example, `--header cancer` will lead to naming '
             'all relevant samples like: cancer_1, cancer_2, etc. '
             'Default: create sample names from first non-empty '
             'line in the file.'
    )

    description_column = Argument(
        short='d',
        action='store_true',
        help='Enable this switch, if there is a column with columns '
             'descriptions (the column has to be on position two, '
             'i.e. immediately after gene identifiers). By default '
             'it is assumed that there is no such column.'
    )

    constructors_by_ext = {
        'tsv': SampleCollection.from_file,
        'csv': SampleCollection.from_csv_file,
        'gct': SampleCollection.from_gct_file
    }

    deduce_format = Argument(
        type=bool,
        default=True,
        help='Deduce file format and automatically set the best '
             'parsing parameters. The format will be inferred from '
             'extension of the provided file(s). '
             f'Following formats are supported: {constructors_by_ext}. '
             'Default: true.'
    )

    def produce(self, unknown_args=None):
        opts = self.namespace
        name = opts.name or self.name

        if opts.files:
            # load all files
            sample_collections = []

            if callable(opts.header):
                opts.header = [opts.header(f) for f in opts.files]

            for i, file_obj in enumerate(opts.files):

                use_header = isinstance(opts.header[i], int)

                constructor = SampleCollection.from_file

                if opts.deduce_format:
                    extension = Path(file_obj.name).suffix[1:]
                    if extension in self.constructors_by_ext:
                        constructor = self.constructors_by_ext[extension]

                sample_collections.append(
                    constructor(
                        f'Sample collection, part {i} of {name}',
                        file_obj,
                        columns_selector=opts.columns[i].get_iterator if opts.columns else None,
                        samples=opts.samples[i] if opts.samples else None,
                        reverse_selection=getattr(opts, 'reverse', False),
                        delimiter=opts.delimiter,
                        header_line=opts.header[i] if use_header else None,
                        use_header=use_header,
                        prefix=opts.header[i] if not use_header else None,
                        description_column=opts.description_column
                    )
                )

            opts.sample_collection = sum(sample_collections, SampleCollection(name))
        return opts


class SingleFileExperimentFactory(Parser):
    """Provide both: case and control samples from a single file.

    This is just a shortcut for specifying the same file for both:
    case and control samples sets. You have to provide --case or
    --control (or both) to specify which columns contain controls.

    If you specify only one of --case and --control, it will be
    assumed that all other columns should be used for the other
    set of samples (if you use `--case 0,1,2` and your file has
    five columns with samples, then columns three and four will
    be used to create control samples).

    To enable more advanced features, please use `control`&`case`
    options (instead of the currently selected `data` sub-parser).
    """

    # exactly one file is required
    files = Argument(
        type=argparse.FileType('r'),
        nargs=1,    # transforms result into a single-element list
        optional=False,
        help='file with samples for both control and cases.'
    )
    case = Argument(
        type=one_of(Slice, Indices, Range),
        nargs=1,
        help='columns from which case samples should be extracted.'
    )
    control = Argument(
        type=one_of(Slice, Indices, Range),
        nargs=1,
        help='columns from which control samples should be extracted.',
    )

    def produce(self, unknown_args=None):

        opts = self.namespace

        def produce_collection_of_samples(created_group, other_group):
            reverse = hasattr(opts, 'reverse_' + created_group)
            get_columns_from = created_group
            if reverse:
                get_columns_from = other_group

            return SampleCollectionFactory(
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
                        'and which should be used as the case.'
                    )

            collections = {
                'control': produce_collection_of_samples('control', 'case')
            }
            # reuse the same file(s)
            for f in opts.files:
                f.seek(0)
            collections['case'] = produce_collection_of_samples('case', 'control')

            for name, sample_collection in collections.items():
                setattr(opts, name, sample_collection)

        return opts


class CLIExperiment(Parser):
    """Use both: case and control or data to create an Experiment."""

    __pull_to_namespace_above__ = True
    __skip_if_absent__ = False

    control = SampleCollectionFactory()
    case = SampleCollectionFactory()
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
        elif opts.case:
            raise ValueError('Control has not been provided!')
        elif opts.control:
            raise ValueError('Case has not been provided!')
        else:
            raise ValueError('Neither data nor (case & control) have been provided!')

        del opts.data

        opts.experiment = Experiment(opts.case.sample_collection, opts.control.sample_collection)

        return opts


class CLI(Parser):
    """The main parser, the one exposed directly to the user."""

    method_name = Argument(choices=Method.members, name='method', optional=False)
    experiment = CLIExperiment()
    __parsing_order__ = 'breadth-first'

    @staticmethod
    def create_method(name):
        # first - take an appropriate method class
        method = Method.members[name]

        # initialize parser for this method
        # (different methods require different arguments)
        method_parser = ConstructorParser(constructor=method)

        return method_parser

    def parse_args(self, args):
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
                    name: ConstructorParser(constructor=method)
                    for name, method in Method.members.items()
                }

                def match_parser(subparsers):
                    return subparsers.get(name, None)

                all_subparsers = [methods, self.subparsers, self.lifted_parsers]

                for parser in filter(bool, map(match_parser, all_subparsers)):
                    return parser.parse_args(args_without_help[1:] + ['-h'])

        return super().parse_args(args)

    def produce(self, unknown_args):
        options = self.namespace

        method_parser = self.create_method(options.method)

        # parse arguments
        method_options, remaining_unknown_args = method_parser.parse_known_args(unknown_args)

        for argument in unknown_args[:]:
            if argument not in remaining_unknown_args:
                unknown_args.remove(argument)

        # and initialize the method with these arguments
        options.method = method_parser.constructor(**vars(method_options))

        return options

