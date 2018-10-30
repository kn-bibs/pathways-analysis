"""
KEGG database interaction types.

Weights based on:
'Onto-Tools: new additions and improvements in 2006',
Khatri P, Voichita C, Kattan K, Ansari N, Khatri A, Georgescu C, Tarca AL, Draghici S.
Nucleic Acids Res., 2007, vol. 35 (pg. W206-W211)
"""

interaction_weights = {
    'activation': 1.0,
    'complex': 1.0,
    'compound': 1.0,
    'dissociation': 1.0,
    'glycosylation': 1.0,
    'indirect': 1.0,
    'indirect effect': 1.0,
    'phosphorylation': 1.0,
    'state': 1.0,
    'binding/association': 1.0,
    'dephosphorylation': 1.0,
    'expression': 1.0,
    'inhibition': -1.0,
    'methylation': 1.0,
    'repression': -1.0,
    'ubiquination': 1.0
}

# maximum impact factor
MAX_IF = 10e200
