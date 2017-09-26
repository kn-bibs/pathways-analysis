import numpy as np
import pandas as pd


class Gene:

    instances = {}

    def __new__(cls, name, *args, **kwargs):
        if name not in cls.instances:
            cls.instances[name] = super(Gene, cls).__new__(cls, *args, **kwargs)
        return cls.instances[name]

    def __init__(self, name):
        self.name = name


# TODO description about what is data?
class Sample:
        self.data = data
        self.name = name

    @classmethod
    def from_names(cls, name, data):
        return cls(name, {Gene(_name) : value for _name, value in data.items()})

    def as_array(self):
        """

        Returns: one-dimensional labeled array with Gene objects as labels

        """
        return pd.Series(self.data)


# TODO class variable with set of genes + method(s) for checking data integrity
class Phenotype:
    def __init__(self, name, samples=None):
        self.samples = samples or []
        self.name = name

    def as_array(self):
        """

        Returns: pandas DataFrame object for all samples, which can be passed to ttest_ind

        """
        return {s.name: pd.DataFrame(s) for s in self.samples}

    def __add__(self, other):
        return self.samples + other.samples

# TODO class variable with set of genes + method(s) for checking data integrity
# TODO unify file reading with argument_parser
class Experiment:

    def __init__(self, case: Phenotype, control: Phenotype):
        self.control = control
        self.case = case

    def get_all(self):
        return self.control + self.case

    @classmethod
    def from_tsv(cls, file_name):
        # TODO
        control = Phenotype(name='Some Hardcore Tumor', samples=[Sample()])
        case = Phenotype(name='reference', samples=[Sample()])
        return cls(control, case)
        """
        def read_tsv_data(file_name, control_samples):
            ""
            Args:
                # TODO
                file_name: file containing gene expression of the following structure:
                - names of samples separated by tab in first row
                - gene symbol/name followed by gene expression values for every sample in remaining rows
                control_samples: list of control_samples indexes (1 - based)
            Returns:
                data frames for two sets of samples
        ""
        with open(file_name, 'r') as f:
            case_samples = [x for x in range(len(f.readline().split('\t')[1:])) if x not in control_samples]

        control = pd.read_csv(file_name, delimiter='\t', header=0, index_col=0, usecols=[0] + list(control_samples))
        case = pd.read_csv(file_name, delimiter='\t', header=0, index_col=0, usecols=case_samples)

        return control, case
        """

    @classmethod
    def from_gsea_file(cls):
        pass

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
    def __init__(self, cases: ArrayType, control: Phenotype):
        for case in cases:
            self.experiments = Experiment(case, control)

