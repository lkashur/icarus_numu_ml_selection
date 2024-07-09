import numpy as np
import uproot
import pandas as pd
from ROOT import TFile, TTree
from array import array

header = ['run', 'event', 'subrun', 'nu_id', 'image_id', 'id', 'trigger', 'category', 'category_topology', 'category_interaction_mode',
           'true_visible_energy', 'reco_visible_energy', 'true_pi0_leading_photon_energy', 'reco_pi0_leading_photon_energy',
           'true_pi0_subleading_photon_energy', 'reco_pi0_subleading_photon_energy', 'true_pi0_costheta', 'reco_pi0_costheta',
           'true_pi0_leading_photon_start_to_vertex', 'reco_pi0_leading_photon_start_to_vertex', 'true_pi0_subleading_photon_start_to_vertex',
           'reco_pi0_subleading_photon_start_to_vertex', 'true_pi0_mass', 'reco_pi0_mass', 'selected_1muNph', 'crtpmt_match', 'cryostat']

def read_log(path, tag, header, category_selector=None):
    """
    Reads an input log file and extracts lines with the specified tag
    into a Pandas DataFrame.

    Parameters
    ----------
    path: str
        The full path of the input log file.
    tag: str
        The identifier that tags relevant lines in the log file.
    header: list[str]
        The list of column names for the CSV file.
    category_selector: callable
        A function that takes a line and returns True if the line is to be included in the output DataFrame.

    Returns
    -------
    data: Pandas.DataFrame
        The DataFrame containing the requested information.
    """
    input_file = open(path)
    lines = input_file.readlines()
    selected = [x.strip('\n').split(',')[1:] for x in lines if tag in x]
    selected = [x if x[-1] != '' else x[:-1] for x in selected]
    data = pd.DataFrame(selected, columns=header[:len(selected[0])])
    for k in header[:len(data.columns)]:
        data[k] = pd.to_numeric(data[k], errors='coerce', downcast='float')
        if data[k].apply(float.is_integer).all():
            data[k] = data[k].astype(int)
    if category_selector is not None:
        data = data[data.apply(category_selector, axis=1)]
    return data

name = 'output_mc_v09_84_01_01r3'
output = TFile(f'{name}.root', 'recreate')

input_name = f'/exp/icarus/app/users/lkashur/CAFAna/icarus_numu_ml_selection/output_mc_v09_84_01_01r3.log'
for channel in ['1muNph']:
    tree = TTree(f'selected_{channel}', f'selected_{channel}')
    df = read_log(input_name, 'SELECTED', header, lambda x: ((x['crtpmt_match'] == 1) & (x[f'selected_{channel}'] == 1) & (np.abs(x['trigger'] - 1500) < 11)))

    vars = [array('d', [0]) for _ in header]
    for ni, n in enumerate(header):
        tree.Branch(n, vars[ni], f'{n}/D')

    print(df.head())
    for i in range(len(df)):
        for ni, n in enumerate(header):
            vars[ni][0] = df[n].iloc[i]
        tree.Fill()
    tree.Write()

output.Close()
