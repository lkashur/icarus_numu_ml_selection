import numpy as np
import pandas as pd
import awkward as ak
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle, Patch
import seaborn as sns
sns.set_context('talk')
import argparse
import uproot
import toml
from tqdm import tqdm


def read_log(log, tag, cfg):
    input_file = open(log)
    lines = input_file.readlines()
    selected = [x.strip('\n').split(',')[1:] for x in lines if tag in x]
    selected = [x if x[-1] != '' else x[:-1] for x in selected]
    header = cfg['data']['columns']
    df = pd.DataFrame(selected, columns=header[:len(selected[0])])
    for k in header[:len(df.columns)]:
        df[k] = pd.to_numeric(df[k], errors='coerce', downcast='float')
        if df[k].apply(float.is_integer).all():
            df[k] = df[k].astype(int)
    return df

def plot_hist(selected, var, cfg):
    fig, ax = plt.subplots(figsize=(8,6))
    
    var_bins = cfg['variables'][var]['bins']
    
    contents = list()
    centers = list()
    labels = list()
    width = list()
    for i,m in enumerate(cfg['categories']['topos']):
        mask = np.isin(selected['category_topology'], m)
        c, e = np.histogram(selected[var][mask], bins=int(var_bins[0]), range=var_bins[1:])
        contents.append(c)
        centers.append((e[1:] + e[:-1]) / 2.0)
        width.append(np.diff(e))
        labels.append(cfg['categories']['topo_labels'][i])

    colors = cfg['categories']['topo_colors']
    ax.hist(centers, weights=contents, bins=int(var_bins[0]), range=var_bins[1:], label=labels, color=colors, histtype='barstacked', ec='white', lw=0.25)
    h, l = ax.get_legend_handles_labels()

    show_percentage = cfg['variables'][var]['show_percentage']
    if show_percentage:
        l = [f'{l} ({np.sum(contents[li]):.0f}, {np.sum(contents[li]) / np.sum(contents):.02%})'for li, l in enumerate(l)]
    else:
        l = [f'{l} ({np.sum(contents[li]):.0f})'for li, l in enumerate(l)]
    h = h[::-1]
    l = l[::-1]
    ax.legend(h, l, fontsize=14)
    ax.set_xlabel(cfg['variables'][var]['xlabel'])
    ax.set_ylabel(cfg['variables'][var]['ylabel'])
    ax.minorticks_on()
    ax.tick_params(which='major', length=6, width=1.5, direction='in')
    ax.tick_params(which='minor', length=3, width=1.5, direction='in')
    plt.title(cfg['variables'][var]['title'], fontsize=25)
    plt.show()

def main(args):

    # Config
    cfg = toml.load(args.cfg)

    # Load log data as pandas dataframe
    selected = read_log(args.log, 'SELECTED', cfg)

    # Plot
    #plot_hist(selected, 'reco_pi0_costheta', cfg)
    #plot_hist(selected, 'reco_pi0_leading_photon_energy', cfg)
    plot_hist(selected, 'reco_pi0_leading_photon_start_to_vertex', cfg)
    #plot_hist(selected, 'reco_pi0_mass', cfg)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cfg', required=True)
    parser.add_argument('--log', required=True)
    args = parser.parse_args()
    main(args)
