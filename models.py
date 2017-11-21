from functools import lru_cache
from typing import Callable, Mapping, Sequence, List
from warnings import warn

import pandas as pd


class Gene:
    """Stores gene's identifier and description (multiton).

    At a time there can be only one gene with given identifier,
    i.e. after the first initialization, all subsequent attempts
    to initialize a gene with the same identifier will return
    exactly the same object. This is so called multiton pattern.

    Example:

        >>> x = Gene('TP53')
        >>> y = Gene('TP53')
        >>> assert x is y   # passes, there is only one gene
    """

    instances = {}

    def __new__(cls, *args, **kwargs):
        if not args:
            # for pickling the requirements are lessened
            # ONLY for pickling
            return super(Gene, cls).__new__(cls)
        name = args[0]
        if name not in cls.instances:
            cls.instances[name] = super(Gene, cls).__new__(cls)
            # TODO: test for description & name
            cls.instances[name].__init__(*args, **kwargs)
        return cls.instances[name]

    def __init__(self, name, description=None):
        self.name = name
        self.description = description

    def __repr__(self):
        return f'<Gene: {self.name}>'


class Sample:
    """Sample contains expression values for genes."""

    def __init__(self, name, data: Mapping[Gene, float]):
        self.name = name
        self.data = data

    @property
    def genes(self):
        return self.data.keys()

    @classmethod
    def from_names(cls, name, data: Mapping[str, float]):
        """Create a sample from a gene_name: value mapping.

        Args:
            name: name of sample
            data: mapping (e.g. dict) where keys represent gene names
        """
        return cls(name, {Gene(gene_name): value for gene_name, value in data.items()})

    @classmethod
    def from_array(cls, name, panda_series: pd.Series, descriptions=False):
        """Create a sample from pd.Series or equivalent.

        Args:
            name: name of the sample
            panda_series:
                series object where columns represent values of genes and
                names are either gene identifiers of tuples:
                ``(gene_identifier, description)``
            descriptions:
                are descriptions present in names of the series object?
        """
        gene_maker = Gene

        if descriptions:
            gene_maker = lambda data: Gene(*data)

        return cls(name, {
            gene_maker(key): value
            for key, value in panda_series.to_dict().items()
        })

    def as_array(self):
        """

        Returns: one-dimensional labeled array with Gene objects as labels

        """
        return pd.Series(self.data)

    def __eq__(self, other):
        return self.name == other.name and self.data == other.data

    def __repr__(self):
        return f'<Sample "{self.name}" with {len(self.data)} genes>'


def first_line(file_object, skip_rows):
    line = None

    while not (line and skip_rows < 0):
        line = file_object.readline()
        if line:
            skip_rows -= 1

    # return to the beginning
    file_object.seek(0)

    return line


# TODO class variable with set of genes + method(s) for checking data integrity
class SampleCollection:
    """A collection of samples of common origin or characteristic.

    An example sample_collection can be:
        (Breast_cancer_sample_1, Breast_cancer_sample_2) named "Breast cancer".

        The common origin/characteristics for "Breast cancer" sample_collection could be
        "a breast tumour", though samples had been collected from two donors.

    Another example are controls:
        (Control_sample_1, Control_sample_2) named "Control".

        The common characteristic for these samples is that both are controls.
    """

    def __init__(self, name, samples=None):
        self.samples: List[Sample] = samples or []
        self.name = name

    @property
    def labels(self):
        return [sample.name for sample in self.samples]

    @property
    @lru_cache(maxsize=1)
    def genes(self):
        """Return all genes present in the collection of samples.

        Raises AssertionError if there is inconsistency in genes in samples.
        """
        genes = self.samples[0].genes
        assert all(sample.genes == genes for sample in self.samples[1:])
        return genes

    @lru_cache(maxsize=None)
    def of_gene(self, gene):
        return tuple(
            sample.data[gene]
            for sample in self.samples
        )

    def as_array(self):
        """
        Returns: :class:`pandas.DataFrame` object with data for all samples.
        """
        return {s.name: pd.DataFrame(s) for s in self.samples}

    def __add__(self, other):
        return SampleCollection(self.name, self.samples + other.samples)

    @classmethod
    def from_file(
            cls, name, file_object,
            columns_selector: Callable[[Sequence[int]], Sequence[int]]=None,
            samples=None, delimiter: str='\t', index_col: int=0,
            use_header=True, reverse_selection=False, prefix=None,
            header_line=0, description_column=None
    ):
        """Create a sample_collection (collection of samples) from csv/tsv file.

        Args:
            name:
                a name of the sample_collection (or group of samples) which will
                identify it (like "Tumour_1" or "Control_in_20_degrees")

            file_object: a file (containing gene expression)
                of the following structure:
                    - names of samples separated by a tab in the first row,
                    - gene symbol/name followed by gene expression values
                      for every sample in remaining rows;

                an additional column "description" is allowed between genes
                column and sample columns, though it has to be explicitly
                declared with `description_column` argument.

            columns_selector:
                a function which will select (and return) a subset of
                provided column identifiers (do not use with `samples`)

            samples:
                a list of names of samples to extract from the file
                (do not use with `columns_selector`)

            reverse_selection:
                if you want to use all columns but the selected ones
                (or all samples but the selected) set this to True

            delimiter: the delimiter of the columns

            index_col: column to use as the gene names

            use_header: does the file have a header?

            prefix: prefix for custom samples naming schema

            header_line: number of non-empty line with sample names

            description_column:
                is column with description of present in the file
                (on the second position, after gene identifiers)?
        """
        if file_object.tell() != 0:
            warn(f'Passed file object: {file_object} was read before.')
            raise Exception()

        line = first_line(file_object, header_line)
        header_items = [item.strip() for item in line.split('\t')]
        gene_columns = [index_col]

        if description_column:
            description_column = 1
            gene_columns.append(description_column)
        else:
            if any('description' == name.lower() for name in header_items):
                warn(
                    'First line of your file contains "description" column, '
                    'but you did not provide "--description_column" argument.'
                )

        # a reasonable assumption is that the columns with samples
        # start after columns with gene symbol and gene description
        column_shift = max(gene_columns) + 1

        if columns_selector:
            # sniff how many columns do we have in the file
            columns_count = line.count(delimiter)

            # TODO test for this change
            all_sample_columns = list(range(column_shift, columns_count))

            # generate identifiers (numbers) for all columns
            # and take the requested subset
            columns = columns_selector(all_sample_columns)

            if reverse_selection:
                columns = list(columns)
                columns = [c for c in all_sample_columns if c not in columns]

            # https://github.com/pandas-dev/pandas/issues/9098#issuecomment-333677100
            columns = gene_columns + list(columns)
        else:
            columns = None

        if not use_header:
            if samples:
                raise ValueError(
                    'To select samples by their name, you need a file with '
                    'samples names in the header. If you use such file, '
                    'please set `use_header=True`, otherwise skip `samples` '
                    'in your arguments.'
                )
            if header_line:
                warn(
                    '`header_line` has no effect when '
                    '`use_header` is set to `False`'
                )

        # we could leave it to pandas, but it shows an ugly,
        # not very helpful message. It is better to show the
        # user where exactly the problem occurs.
        if samples:

            available_samples = [
                name
                for name in header_items[column_shift:]
            ]

            lacking_samples = set(samples) - set(available_samples)

            if lacking_samples:
                raise ValueError(
                    f'Samples {lacking_samples} are not available in {file_object.name} file.\n'
                    f'Following samples were found: {", ".join(available_samples)}.'
                )

            if index_col:
                # TODO https://github.com/pandas-dev/pandas/issues/9098
                warn(
                    'Using "samples" with "index_col" 0 may cause an '
                    'unexpected behaviour due to an upstream issue in '
                    'pandas package (pandas-dev/pandas/issues/9098) '
                    'for pandas in versions older than 0.21.'
                )

            additional_column_names = [
                header_items[index]
                for index in gene_columns
            ]

            # https://github.com/pandas-dev/pandas/issues/9098#issuecomment-333677100
            samples = additional_column_names + list(samples)

        # just to reassure that the pointer is on the beginning
        if file_object.tell() != 0:
            warn('Passed file object was read before.')

        if samples and columns:
            warn(
                'Please, provide either columns or samples, '
                'not both. We will use columns this time.'
            )

        try:
            data = pd.read_table(
                file_object,
                delimiter=delimiter,
                # None - do not use, 0 - use first row
                header=header_line if use_header else None,
                index_col=gene_columns,
                usecols=columns or samples,
                prefix=f'{prefix}_' if prefix else ''
            )
        except Exception as e:
            from traceback import print_tb
            from traceback import print_stack
            print_tb(e)
            print(e)
        descriptions = description_column is not None

        samples = [
            Sample.from_array(sample_name, sample_data, descriptions=descriptions)
            for sample_name, sample_data in data.items()
        ]

        return cls(name, samples)

    @classmethod
    def from_gct_file(cls, name, file_object, **kwargs):
        """Parse file in Gene Cluster Text file format, as defined on:

        software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats
        User is allowed to provide settings different from the standard.
        """
        version = file_object.readline()
        rows_count, samples_count = map(int, file_object.readline().split('\t'))

        default_values = {
            'description_column': True,
            'header_line': 2
        }

        if version != '#1.2\n':
            warn('Unsupported version of GCT file')

        file_object.seek(0)

        for key, value in default_values.items():
            kwargs[key] = value

        self = cls.from_file(
            name, file_object,
            **kwargs
        )

        # if user did not choose a subset of samples
        if not (kwargs['columns_selector'] or kwargs['samples']):
            # check if the samples numbers are ok
            if len(self.samples) != samples_count:
                warn(
                    f'Samples count ({len(self.samples)}) '
                    f'does not match with the {samples_count} '
                    f'declared in {name} file.'
                )

        if rows_count != len(self.samples[0].genes):
            warn(
                f'Number of rows ({rows_count}) does not match '
                f'with the {len(self.samples[0].genes)} '
                f'declared in {name} file'
            )

        return self

    @classmethod
    def from_csv_file(cls, **kwargs):
        if 'delimiter' in kwargs:
            if kwargs['delimiter'] != ',':
                raise Warning(
                    'You are using not comma delimiter for what looks like csv file. '
                    'Is this really the thing you want to do?'
                )
        else:
            kwargs['delimiter'] = ','
        return cls.from_file(**kwargs)


# TODO class variable with set of genes + method(s) for checking data integrity
class Experiment:

    def __init__(self, case: SampleCollection, control: SampleCollection):
        self.control = control
        self.case = case

    def get_all(self):
        return self.control + self.case

    # TODO: are there many ways to compute fold-change?
    def get_fold_change(self, sample_from_case, use_log=False):
        assert sample_from_case in self.case.samples
        # TODO: implement inline
        calc_fold_change(sample_from_case, self.control, use_log=use_log)
        """
        def fold_change(case, base, log2=False):
            fold_changes = case.copy()
            for (idx, row) in base.iterrows():
                fold_changes.loc[[idx]] /= (np.mean(row) or 0.01)  # TODO for now arbitrary value 0.01 when 0's are found

            if log2:
                fold_changes = np.log2(fold_changes)  # TODO Runtime Warning when 0's are encountered

            return fold_changes
        """


class Study:
    def __init__(self, cases: Sequence[SampleCollection], control: SampleCollection):
        for case in cases:
            self.experiments = Experiment(case, control)

