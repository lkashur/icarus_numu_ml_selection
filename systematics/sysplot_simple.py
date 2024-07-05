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
from ROOT import TFile, TEfficiency, TH1D, TGraphAsymmErrors


def read_log(log, tag, cfg, data_or_mc):
    input_file = open(log)
    lines = input_file.readlines()
    tagged = [x.strip('\n').split(',')[1:] for x in lines if tag in x]
    tagged = [x if x[-1] != '' else x[:-1] for x in tagged]
    if data_or_mc == 'mc':
        header = cfg['data']['columns_mc']
    elif data_or_mc == 'data':
        header = cfg['data']['columns_data']
    df = pd.DataFrame(tagged, columns=header[:len(tagged[0])])
    for k in header[:len(df.columns)]:
        df[k] = pd.to_numeric(df[k], errors='coerce', downcast='float')
        if df[k].apply(float.is_integer).all():
            df[k] = df[k].astype(int)
    return df

def get_stat_cov(bins):
    stat_cov = np.diag(bins)
    stat_err = np.sqrt(np.diagonal(stat_cov))
    return stat_cov, stat_err

def plot_hist(cfg, selected_mc, selected_data, var, norm_mode='data'):
    fig, ax = plt.subplots(figsize=(8,6))
    
    show_percentage = cfg['variables'][var]['show_percentage']
    var_bins = cfg['variables'][var]['bins']

    # MC (Nu + Cosmic)
    contents_mc = []
    centers_mc = []
    labels_mc = []
    width_mc = []
    for i,m in enumerate(cfg['categories']['topos']):
        mask = np.isin(selected_mc['category_topology'], m)
        b, e = np.histogram(selected_mc[var][mask], bins=int(var_bins[0]), range=var_bins[1:])
        contents_mc.append(b)
        centers_mc.append((e[1:] + e[:-1]) / 2.0)
        width_mc.append(np.diff(e))
        labels_mc.append(cfg['categories']['topo_labels'][i])
    if norm_mode == 'data':
        contents_mc = [(c / cfg['data']['pot_mc']) * cfg['data']['pot_data'] for c in contents_mc]

    # Data (on-beam)
    bins_data, edges_data = np.histogram(selected_data[var], bins=int(var_bins[0]), range=var_bins[1:])
    centers_data = (edges_data[1:] + edges_data[:-1]) / 2.0
    width_data = np.diff(edges_data)
    label_data = f'Data ({np.sum(bins_data):.0f})' if show_percentage else 'Data' 
    stat_cov_data, stat_err_data = get_stat_cov(bins_data)

    # Plotting
    colors = cfg['categories']['topo_colors']
    ax.hist(centers_mc, weights=contents_mc, bins=int(var_bins[0]), range=var_bins[1:], label=labels_mc, color=colors, histtype='barstacked', ec='white', lw=0.25)
    ax.errorbar(centers_data, bins_data, xerr=width_data/2, yerr=stat_err_data, label=label_data, linestyle='none', capsize=2, ecolor='black')
    h, l = ax.get_legend_handles_labels()
    h_mc, l_mc = h[:-1], l[:-1]

    if show_percentage:
        l_mc = [f'{l} ({np.sum(contents_mc[li]):.0f}, {np.sum(contents_mc[li]) / np.sum(contents_mc):.02%})'for li, l in enumerate(l_mc)]
        
    else:
        l_mc = [f'{l} ({np.sum(contents_mc[li]):.0f})'for li, l in enumerate(l_mc)]

    h = [h[-1]] + h_mc
    l = [l[-1]] + l_mc

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

def calc_flat_efficiency(signal):
    total_signal_events = len(signal)
    selected_signal_events = len(signal[signal['selected_1muNph'] == 1])
    return 100 * selected_signal_events / total_signal_events

def calc_flat_purity(selected):
    total_selected_events = len(selected)
    matched_selected_events = len(selected[selected['category'] == 0])
    return 100 * matched_selected_events / total_selected_events
        
def plot_diff_pur_eff(df, var, var_range, quantiles, metric):
    bcs = []
    bxerr0s = []
    bxerr1s = []
    byerr0s = []
    byerr1s = []
    vals = []
    if var_range[0] == var_range[-1]:
        df = df
    else:
        df = df[(df[var] > var_range[0]) & (df[var] < var_range[1])]
    df['var_q'] = pd.qcut(df[var], q=quantiles)
    
    hedges = sorted([i.left for i in df.var_q.unique().tolist()])
    lastedge = max([i.right for i in df.var_q.unique().tolist()])
    hedges.append(lastedge)
    hedges = np.array(hedges)
    hpass = TH1D('hpass', '', quantiles, hedges)
    htotal = TH1D('htotal', '', quantiles, hedges)
    
    for i,(name,group) in enumerate(df.groupby('var_q')):
        #print(name)
        #print(group)
        bxerr0s.append(name.left)
        bcs.append(name.mid)
        bxerr1s.append(name.right)
        if metric == 'pur':
            vals.append(calc_flat_purity(group) / 100)
            hpass.SetBinContent(i+1, len(group[group['category'] == 0]))
        elif metric == 'eff':
            vals.append(calc_flat_efficiency(group) / 100)
            hpass.SetBinContent(i+1, len(group[group['selected_1muNph'] == 1]))
        htotal.SetBinContent(i+1, len(group))
            
    gr = TGraphAsymmErrors()
    gr.Divide(hpass, htotal, 'cl=0.683 b(1,1) mode')
    for i in range(10):
        byerr0s.append(gr.GetErrorYlow(i))
        byerr1s.append(gr.GetErrorYhigh(i))

    fig, ax = plt.subplots(figsize=(8,6))
    ax.errorbar(bcs, vals, xerr=[np.array(bcs) - np.array(bxerr0s), np.array(bxerr1s) - np.array(bcs)], yerr=[np.array(byerr0s), np.array(byerr1s)], linestyle='')
    ax.set_ylim([0,1])
    plt.show()


def main(args):

    # Config
    cfg = toml.load(args.cfg)

    # Load log data as pandas dataframe
    selected_mc = read_log(args.log_mc, 'SELECTED', cfg, data_or_mc='mc') # stacked hists, purity
    signal = read_log(args.log_mc, 'SIGNAL', cfg, data_or_mc='mc') # efficiency
    selected_data = read_log(args.log_data, 'DATA', cfg, data_or_mc='data')
   
    #print(f'Signal Efficiency: {round(calc_flat_efficiency(signal), 1)}%')
    #print(f'Selection Purity: {round(calc_flat_purity(selected), 1)}%')

    #plot_diff_pur_eff(signal, 'true_pi0_leading_photon_energy', [0, 600], 10, metric='eff')
    #plot_diff_pur_eff(selected, 'reco_pi0_leading_photon_energy', [0, 600], 10, metric='pur')
    #plot_diff_pur_eff(signal, 'true_visible_energy', [-1, -1], 10, metric='eff')
    
    # Plot
    #plot_hist(selected, 'reco_pi0_costheta', cfg)
    #plot_hist(selected, 'reco_pi0_leading_photon_energy', cfg)
    #plot_hist(selected, 'reco_pi0_subleading_photon_energy', cfg)
    #plot_hist(selected, 'reco_pi0_leading_photon_start_to_vertex', cfg)
    #plot_hist(selected, 'reco_pi0_subleading_photon_start_to_vertex', cfg)
    plot_hist(cfg, selected_mc, selected_data, 'reco_pi0_mass')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cfg', required=True)
    parser.add_argument('--log_mc', required=True)
    parser.add_argument('--log_data', required=True)
    args = parser.parse_args()
    main(args)
