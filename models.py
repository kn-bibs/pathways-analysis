from functools import lru_cache
from metrics import ratio_of_classes
from numpy import log2
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
    __slots__ = ('name', 'description', 'id')

    def __new__(cls, *args, **kwargs):
        if not args:
            # for pickling the requirements are lessened
            # ONLY for pickling
            return super(Gene, cls).__new__(cls)

        name = args[0]

        if name not in cls.instances:
            gene = super(Gene, cls).__new__(cls)
            gene.__init__(*args, **kwargs)
            gene.id = len(cls.instances) - 1
            cls.instances[name] = gene

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
        Returns:
            one-dimensional labeled array with Gene objects as labels

        """
        return pd.Series(self.data)

    def __eq__(self, other):
        return self.name == other.name and self.data == other.data

    def __repr__(self):
        return f'<Sample "{self.name}" with {len(self.data)} genes>'

    def exclude_genes(self, gene_list: list):
        for gene in gene_list:
            assert isinstance(gene, Gene)
            if gene in self.data.keys():
                del self.data[gene]


def first_line(file_object, skip_rows=0):
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

    def __init__(self, name: str, samples=None):
        self.samples: List[Sample] = samples or []
        self.name = name
        # integrity check
        # Raises AssertionError if there is inconsistency in genes in samples.
        # genes = self.samples[0].genes
        # assert all(sample.genes == genes for sample in self.samples[1:])

    @property
    def labels(self):
        return [sample.name for sample in self.samples]

    @property
    @lru_cache(maxsize=1)
    def genes(self):
        """
        Returns:
            all genes present in the collection of samples.
        """
        genes = self.samples[0].genes
        return genes

    @lru_cache(maxsize=None)
    def of_gene(self, gene):
        return tuple(
            sample.data[gene]
            for sample in self.samples
        )

    def as_array(self):
        """
        Returns:
            `pandas.DataFrame`: two-dimensional labeled array with Gene objects as row labels,
            storing data from all samples
        """
        df = pd.DataFrame()
        for sample in self.samples:
            if df.empty:
                df = sample.as_array().to_frame(sample.name)
            else:
                kwargs = {sample.name: sample.as_array().values}
                df = df.assign(**kwargs)
        return df

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
                None - do not use, 0 - use first row

            description_column:
                is column with description of present in the file
                (on the second position, after gene identifiers)?
        """
        if file_object.tell() != 0:
            warn(f'Passed file object: {file_object} was read before.')
            raise Exception()

        line = first_line(file_object, header_line or 0)
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

            all_sample_columns = list(range(column_shift, columns_count + 1))

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
        if not any(key in kwargs for key in ['samples', 'columns_selector']):
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
    def from_csv_file(cls, name, file_object, **kwargs):
        if 'delimiter' in kwargs:
            if kwargs['delimiter'] != ',':
                warn(
                    'You are using not comma delimiter for what looks like csv file. '
                    'Is this really the thing you want to do?'
                )
        else:
            kwargs['delimiter'] = ','
        return cls.from_file(name, file_object, **kwargs)

    def exclude_genes(self, gene_list: list):
        for sample in self.samples:
            sample.exclude_genes(gene_list)


# TODO class variable with set of genes + method(s) for checking data integrity
class Experiment:
    """
    Stores all user's experiment data.
    """

    def __init__(self, case: SampleCollection, control: SampleCollection):
        self.control = control
        self.case = case

    def get_all(self):
        return self.control + self.case

    def calculate_fold_change(self):
        """

        Returns:
            `pandas.DataFrame` object: two-dimensional labeled array with Gene objects as row labels, storing
            fold change and log transformed fold change values for every gene - fold change of the expression
            level of given gene in the sample under study to the normal level (average in a control group)

        """
        fc = {}
        for (idx, row) in self.get_all().as_array().iterrows():
            control = [row[label] for label in self.control.labels]
            case = [row[label] for label in self.case.labels]
            ratio = ratio_of_classes(case, control)
            fc[idx] = [ratio, log2(ratio)]
        return pd.DataFrame.from_dict(fc, orient="index", columns=['FC', 'logFC'])

    def exclude_genes(self, gene_list: list):
        self.get_all().exclude_genes(gene_list)
