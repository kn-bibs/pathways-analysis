import pandas as pd
import numpy as np


def read_csv(file, control_samples):
    """
    :param file: file containing gene expression of the following structure:
        - names of samples separated by tab in first row
        - gene symbol/name followed by gene expression values for every sample in remaining rows
    :param control_samples: list of control_samples indexes (1 - based)
    :return: data frames for two sets of samples
    """
    with open(file, 'r') as f:
        case_samples = [x for x in range(len(f.readline().split('\t')[1:])) if x not in control_samples]

    control = pd.read_csv(file, delimiter='\t', header=0, index_col=0, usecols=[0] + list(control_samples))
    case = pd.read_csv(file, delimiter='\t', header=0, index_col=0, usecols=case_samples)

    return control, case

# con, case = read_csv('sample_expression.tsv', control_samples=range(3, 6))
# print(con.loc[['MYH16']])  # row
# print(con['Tumour_sample_1']) #column


def fold_change(case,base, log2 = False):
    """

    :param case:
    :param base:
    :param log:
    :return:
    """
    fold_changes = case.copy()
    for (idx,row) in base.iterrows():
        fold_changes.loc[[idx]] /= (np.mean(row) or 0.01) # TODO for now arbitrary value 0.01 when 0's are found

    if log2:
        fold_changes = np.log2(fold_changes) #TODO Runtime Warning when 0's are encountered

    return fold_changes

# f_cs = fold_change(case,con)
# print(f_cs)
# f_cs = fold_change(case,con,log2=True)
# print(f_cs)
