import skimage as sk
from skimage import io
import skimage.external.tifffile as tff
import os
import pylab
import numpy as np
from scipy import ndimage as ndi
import scipy as spy
import scipy.spatial.distance
import skimage as sk
import skimage.filters as skf
import skimage.morphology as skm
import skimage.measure as skmes
import skimage.segmentation as sks
import pandas as pd
import matplotlib.patches as mptch
import scipy.io as sio
import scipy as scpy
import glob
import skimage.transform as skt
import networkx as nx
import time
import math
import re
import SimpleITK as sitk
from ipywidgets import interact, fixed
from IPython.display import clear_output
from sklearn.decomposition import PCA
#import community
#import pyintergraph
# from cdlib import algorithms

from skimage.transform import downscale_local_mean
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.metrics.pairwise import cosine_similarity

import matplotlib.pyplot as plt
import scipy.signal as spy_sig

import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

import pickle


### Function definitions
def load_pickle(d_dir,fn):
    pickle_in = open(os.path.join(d_dir,fn),'rb')
    dict_d = pickle.load(pickle_in)
    pickle_in.close()
    return dict_d

def save_pickle(d_dir, fn, v):
    pickle_out = open(os.path.join(d_dir,fn),"wb")
    pickle.dump(v, pickle_out)
    pickle_out.close()

def construct_fn_dict(file_names, key_extr_pat=r'.*cnm_C_key(\d*)_(.*\.tif)\..*'):
    print('No of files: for ', key_extr_pat, '    ', len(file_names))
    dict_fn = {}
    for i in range(0,len(file_names)):
        fn = os.path.basename(file_names[i])
        matchObj = re.match(key_extr_pat, fn, re.M|re.I)
        key = matchObj.group(1)
        key = int(key)
        dict_fn[key] = fn
    dict_fn_sorted = {}
    for k in sorted(dict_fn.keys()):
        dict_fn_sorted[k] = dict_fn[k]
    return dict_fn_sorted

def draw_activity_comp_in_range(act_dat, ax, act_filt_nsp_ids, max_dff_int, begin_tp, end_tp, n_start= 0, n_stop = -1, dff_bar=1, fps=30, lw=.2):
    act_dat1 = np.copy(act_dat)
    if (end_tp < 0):
        end_tp = act_dat1.shape[1]    
    sort_ind =  np.argsort(act_filt_nsp_ids)
    act_dat1 = act_dat1[sort_ind,:]
    act_filt_nsp_ids = act_filt_nsp_ids[sort_ind]
    act_dat1 = act_dat1[:,begin_tp:end_tp]
    n= n_stop 
    if (n_stop == -1):
        n=act_dat1.shape[0]
        n_stop = n
    t=act_dat1.shape[1]
    #max_dff_int = max_dff_int.astype('float')
    max_dff_int = max_dff_int/2
    
    if (max_dff_int==0):
        max_dff_int=0.25
    
    ax.margins(0.008)
#     ax.axis('off')
    lt = len(act_dat1[0,:])
    lt = np.arange(lt)
    lt = lt.astype('float')
    lt = lt/fps
    for i in range(n_start, n_stop):
#         color_str = 'C' + str(act_filt_nsp_ids[i]%9)
        color_str = 'C' + str(i%9)
        ax.plot(lt, act_dat1[i,:]+(n-i-1)*max_dff_int, linewidth =lw, color=color_str)
    ax.vlines(x=-1., ymin=0, ymax=dff_bar, lw=2)
    ax.set_xlabel('Time(s)')
    l1 = np.arange(0,n*max_dff_int, max_dff_int)+max_dff_int/6
    l2 = np.arange(n,0, -1)

def draw_activity_comp_in_range(act_dat, ax, act_filt_nsp_ids, max_dff_int, begin_tp, end_tp, n_start= 0, n_stop = -1, dff_bar=1, fps=30, lw=.2):
    act_dat1 = np.copy(act_dat)
    if (end_tp < 0):
        end_tp = act_dat1.shape[1]    
    sort_ind =  np.argsort(act_filt_nsp_ids)
    act_dat1 = act_dat1[sort_ind,:]
    act_filt_nsp_ids = act_filt_nsp_ids[sort_ind]
    act_dat1 = act_dat1[:,begin_tp:end_tp]
    n= n_stop 
    if (n_stop == -1):
        n=act_dat1.shape[0]
        n_stop = n
    t=act_dat1.shape[1]
    #max_dff_int = max_dff_int.astype('float')
    max_dff_int = max_dff_int/2
    
    if (max_dff_int==0):
        max_dff_int=0.25
    
    ax.margins(0.008)
#     ax.axis('off')
    lt = len(act_dat1[0,:])
    lt = np.arange(lt)
    lt = lt.astype('float')
    lt = lt/fps
    for i in range(n_start, n_stop):
#         color_str = 'C' + str(act_filt_nsp_ids[i]%9)
        color_str = 'C' + str(i%9)
        ax.plot(lt, act_dat1[i,:]+(n-i-1)*max_dff_int, linewidth =lw, color=color_str)
    ax.vlines(x=-1., ymin=0, ymax=dff_bar, lw=2)
    ax.set_xlabel('Time(s)')
    l1 = np.arange(0,n*max_dff_int, max_dff_int)+max_dff_int/6
    l2 = np.arange(n,0, -1)
    
def draw_activity_comp(act_dat, ax, act_filt_nsp_ids, max_dff_int, begin_tp, end_tp, dff_bar=1, fps=30, lw=1.):
    act_dat1 = np.copy(act_dat)
    if (end_tp < 0):
        end_tp = act_dat1.shape[1]    
    sort_ind =  np.argsort(act_filt_nsp_ids)
    act_dat1 = act_dat1[sort_ind,:]
    act_filt_nsp_ids = act_filt_nsp_ids[sort_ind]
    act_dat1 = act_dat1[:,begin_tp:end_tp]    
    n=act_dat1.shape[0]
    t=act_dat1.shape[1]
    #max_dff_int = max_dff_int.astype('float')
    max_dff_int = max_dff_int/2
    
    if (max_dff_int==0):
        max_dff_int=0.25
    
    ax.margins(0.008)
#     ax.axis('off')
    lt = len(act_dat1[0,:])
    lt = np.arange(lt)
    lt = lt.astype('float')
    lt = lt/fps
    for i in range(act_dat1.shape[0]):
        color_str = 'C' + str(act_filt_nsp_ids[i]%9)
        ax.plot(lt, act_dat1[i,:]+(n-i-1)*max_dff_int, alpha=0.75, linewidth =lw, color=color_str)
    ax.vlines(x=-1., ymin=0, ymax=dff_bar, lw=2)
    ax.set_xlabel('Time(s)')
    l1 = np.arange(0,n*max_dff_int, max_dff_int)+max_dff_int/6
    l2 = np.arange(n,0, -1)

def plot_rasterplot_comp_notck(dat_act, nsp_ids, ax, vmin=0, vmax=1, sz_per_neuron=0.1, fs=30, 
                    gauss_sm=False, gauss_ker=3, ylab_tk_sz=8, colmap=False):    
    s = np.copy(dat_act)

    if (s.shape[0]>0):
        if (gauss_sm):
            s = ndi.gaussian_filter1d(s,gauss_ker,axis=1)
        sort_ind =  np.argsort(nsp_ids)
        s = s[sort_ind,:]
        nsp_ids = nsp_ids[sort_ind]

        ticks = list()
        n = s.shape[0]
        ticks_lab = np.unique(nsp_ids)
        for i in ticks_lab:
            ticks.append(np.argwhere(nsp_ids==i)[0])
        if (colmap):
            ax.imshow(s, vmin=vmin,vmax=vmax,aspect='auto',extent=[0,s.shape[1]/fs, s.shape[0],0])
        else:
            ax.imshow(s, vmin=vmin,vmax=vmax,aspect='auto',extent=[0,s.shape[1]/fs, s.shape[0],0], cmap='Greys')
        ax.margins(0.008)
    
        ax.set_yticks([]) 

def draw_graphs_cmap_comp(coms, dat1, dat2, im_lab, norm_corr = False, do_gauss_smooth = False, gauss_ker = 11, 
                     fig_sz = (16, 8), node_size=10, vmin=0.3, vmax=1, cmap_ = plt.cm.Reds, edge_wt=1, 
                          alpha=0.7, cb=True, cb_shr=0.7):
    """
    uses plot_wt_graph_comp to plot the comparative graphs.
    """
    index = list(range(1,coms.shape[0]+1))
    cord = {}
    for i in range(1,coms.shape[0]+1):
        cord[i] = [coms[i-1,1],coms[i-1,0]]
    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=fig_sz, sharex=True, sharey=True)
    ax[0].imshow(im_lab, cmap='Greys', vmin=0, vmax=40)
    ax[1].imshow(im_lab, cmap='Greys', vmin=0, vmax=40)
    if (do_gauss_smooth):
        dat1 = ndi.gaussian_filter1d(dat1,gauss_ker,axis=1)
        dat2 = ndi.gaussian_filter1d(dat2,gauss_ker,axis=1)
    r1 = np.corrcoef(dat1)
    r2 = np.corrcoef(dat2)
#     r1 = np.nan_to_num(r1)
#     r2 = np.nan_to_num(r2)
    if(norm_corr):
        r1 = np.arctanh(r1)
        r2 = np.arctanh(r2)
    np.fill_diagonal(r1,0)
    np.fill_diagonal(r2,0)
    print('index:', len(index), r1.shape, r2.shape)
    plot_wt_graph_cmap(index, r1, ax[0], cord, vmin=vmin, vmax=vmax, node_size=node_size, 
               edge_wt=edge_wt, cb=cb, with_labels=False, cb_shr=cb_shr, cmap_ = cmap_, alpha=alpha)
    plot_wt_graph_cmap(index, r2, ax[1], cord, vmin=vmin, vmax=vmax, node_size=node_size, 
               edge_wt=edge_wt, cb=cb, with_labels=False, cb_shr=cb_shr, cmap_ = cmap_, alpha=alpha)
    ax[0].set_aspect('equal')  
    ax[1].set_aspect('equal')  
    plt.show()
    return fig

def convert_f_to_cs(act_f, p=2, noise_range=[.25, .5]):
    c_arr = np.zeros_like(act_f)
    s_arr = np.zeros_like(act_f)
    for i in range(act_f.shape[0]):
        f = np.copy(act_f[i,:])
        c, bl, c1, g, sn, sp, lam = deconvolution.constrained_foopsi(f, bl=None,  c1=None, g=None,  sn=None, 
                       p=p, method_deconvolution='oasis', bas_nonneg=True,
                       noise_range=noise_range, noise_method='logmexp', lags=5, fudge_factor=1.,
                       verbosity=True, solvers=None, optimize_g=0)
        c_arr[i,:] = c
        s_arr[i,:] = sp
    return c_arr, s_arr

def calc_rise_tm(c, s_zsc, zsc_th=3):
    rise_tm_dat = {}
    rise_tm_pos_dat = {}
    for i in range(c.shape[0]):
        e = s_zsc >= zsc_th
        e_i = np.copy(e[i,:])
        c_i = c[i,:]
        rise_tm_list = []
        rise_tm_pos_list = []
        ind = 0
        while ind < e.shape[1]-10:
            if(e_i[ind]):
                j = ind + 1
                prev = c_i[ind]
                rise_tm = 0
                while (c_i[j] >= prev) & (j < len(c_i)-2):
                    rise_tm = rise_tm + 1
                    prev = c_i[j]
                    j = j + 1
                rise_tm_list.append(rise_tm)
                rise_tm_pos_list.append(j)
                ind = j + 1
            else:
                ind = ind + 1
        rise_tm_dat[i] = rise_tm_list
        rise_tm_pos_dat[i] = rise_tm_pos_list
    return rise_tm_dat, rise_tm_pos_dat

def sort_graph_props_in_DIVs(d_k_in_DIVs, d_g_prop):
    d_g_prop_in_DIVs = {}
    for dk in d_k_in_DIVs.keys():
        print('DIV: ', dk)
        k_in_dk = d_k_in_DIVs[dk]
        list_vals = []
        for k in k_in_dk:
            list_vals.append(d_g_prop[k])
        d_g_prop_in_DIVs[dk] = list_vals
    return d_g_prop_in_DIVs

def calc_pwc_mn(d_k_in_DIVs, d_dff, d_nspIDs, dff_cut=0.1, norm_corr = False):
    d_mn_pwc = {}
    d_mn_pwc_intra = {}
    d_mn_pwc_inter = {}

    for dk in d_k_in_DIVs.keys():
        print('DIV: ', dk)
        k_in_dk = d_k_in_DIVs[dk]
        list_vals = []
        list_vals_intra = []
        list_vals_inter = []

        for k in k_in_dk:
            tmp_dat = np.copy(d_dff[k])
            nsp_ids = d_nspIDs[k]
            nsp_ids_uniq = np.unique(nsp_ids)
            if (0 in nsp_ids_uniq):
                tmp_dat_filt = tmp_dat[nsp_ids>0,:]
                nsp_ids_filt = nsp_ids[nsp_ids>0]
            else:
                tmp_dat_filt = np.copy(tmp_dat)
                nsp_ids_filt = np.copy(nsp_ids)
            tmp_dat_filt[tmp_dat_filt<dff_cut] = 0
            r = np.corrcoef(tmp_dat_filt)
            r = np.nan_to_num(r)
            np.fill_diagonal(r,0)
            if(norm_corr):
                r = np.arctanh(r)
            r[np.isinf(r)]=0
            r_intra, r_inter = extract_intra_inter_nsp_neurons(r, nsp_ids_filt)
            r = r[np.triu_indices(r.shape[0], k = 1)]

#             r[r<0] = 0
#             r_inter[r_inter<0] = 0
#             r_intra[r_intra<0] = 0
            list_vals.append(r.mean())
            list_vals_inter.append(r_inter.mean())
            list_vals_intra.append(r_intra.mean())
        print(list_vals)
        d_mn_pwc[dk] = np.array(list_vals)
        d_mn_pwc_inter[dk] = np.array(list_vals_inter)
        d_mn_pwc_intra[dk] = np.array(list_vals_intra)
    return(d_mn_pwc, d_mn_pwc_inter, d_mn_pwc_intra)

def sort_dict(d):
    d1 = {}
    for k in sorted(d.keys()):
        d1[k] = d[k]
    return d1

def calc_fwhm_spikes(c, s_zsc, zsc_th=3, perc=.2):
    fwhm_dat = {}
    fwhm_pos_back = {}
    fwhm_pos_fwd = {}
    no_spk_dat = {}
    for i in range(c.shape[0]):
        e = s_zsc > zsc_th
        e_i = np.copy(e[i,:])
        c_i = c[i,:]
        fwhm_pos_back_list = []
        fwhm_pos_fwd_list = []
        fwhm_list = []
        no_spk_list = []
        ind = 0
        while ind < e.shape[1]-10:
            if(e_i[ind]):
                j = ind + 1
                prev = c_i[ind]
                while (c_i[j] >= prev) & (j < len(c_i)-2):
                    prev = c_i[j]
                    j = j + 1
                j = j - 1
                prev = perc*(prev - c_i[ind]) + c_i[ind]
                ## Backword ind
                back_i = j - 1
                while (c_i[back_i] >= prev):
                    back_i = back_i - 1
                    if(back_i < 0):
                        break
                back_i = back_i + 1
                fwhm_pos_back_list.append(back_i)
                ##   Forward
                fwd_i = j + 1
                while (c_i[fwd_i] >= prev) & (fwd_i < len(c_i)-10):
                    fwd_i = fwd_i + 1
                fwd_i = fwd_i - 1
                fwhm_pos_fwd_list.append(fwd_i)
                fwhm = fwd_i - back_i + 1
                fwhm_list.append(fwhm)
                no_spk_list.append(1 + np.sum(e_i[back_i:(fwd_i+1)]))
                ind = fwd_i + 1
            else:
                ind = ind + 1
        fwhm_pos_back[i] = fwhm_pos_back_list
        fwhm_pos_fwd[i] = fwhm_pos_fwd_list
        fwhm_dat[i] = fwhm_list
        no_spk_dat[i] = no_spk_list
    return fwhm_pos_back, fwhm_pos_fwd, fwhm_dat, no_spk_dat

def gen_mn_of_mn(d_mn_pwc):
    d = {}
    d_std = {}
    for k in d_mn_pwc:
        d[k] = np.nanmean(d_mn_pwc[k])
        d_std[k] = np.nanstd(d_mn_pwc[k])
    return d, d_std



def gen_poly_fit(d, deg=4):
    x = list(d.keys())
    y = np.array(list(d.values()))
    x = np.array(x)
    m = np.isnan(y)
    y = y[~m]
    x = x[~m]
    x = x[:,np.newaxis]
    poly_reg = PolynomialFeatures(degree=deg)
    x_poly = poly_reg.fit_transform(x)
    pol_reg = LinearRegression()
    pol_reg.fit(x_poly, y)
    return x, pol_reg.predict(poly_reg.fit_transform(x))


def calc_mod_comm(d_k_in_DIVs, d_dff, d_nspIDs, d_com, dff_cut=0.05, norm_corr = False, louvian=True, ccut=0.0, binary=False):
    d_mod = {}
    d_no_comm = {}
    for dk in d_k_in_DIVs.keys():
        print('DIV: ', dk)
        k_in_dk = d_k_in_DIVs[dk]
        list_mod = []
        list_no_comm = []
        for k in k_in_dk:
            tmp_dat = np.copy(d_dff[k])
            nsp_ids = d_nspIDs[k]
            nsp_ids_uniq = np.unique(nsp_ids)
            if (0 in nsp_ids_uniq):
                tmp_dat_filt = tmp_dat[nsp_ids>0,:]
                nsp_ids_filt = nsp_ids[nsp_ids>0]
            else:
                tmp_dat_filt = np.copy(tmp_dat)
                nsp_ids_filt = np.copy(nsp_ids)
            tmp_dat_filt[tmp_dat_filt<dff_cut] = 0
            r = np.corrcoef(tmp_dat_filt)
            r = np.nan_to_num(r)
            np.fill_diagonal(r,0)
            
            coms = d_com[k]
            coms = coms[nsp_ids>0,:]
            index = list(range(1,coms.shape[0]+1))
            g=nx.Graph()
            for i in range(len(index)):
                g.add_node(index[i])
            for i in range(len(index)):
                for j in range(i+1, len(index)):

                    if (r[i,j] > ccut):
                        if(not binary):
                            g.add_edge(index[i],index[j],weight=r[i,j])
                        if(binary):
                            g.add_edge(index[i],index[j])
#             cord = {}
#             for i in range(1,coms.shape[0]+1):
#                 cord[i] = [coms[i-1,1],coms[i-1,0]]
            if (louvian):
                clusters = community.best_partition(g)
                size = float(len(set(clusters.values())))
                mod = community.modularity(clusters, g, weight='weight')
            else:
                gi = pyintergraph.nx2igraph(g)
                clusters = community.best_partition(gi)
                a = clusters.communities
                size = len(a)
                mod = 0
#             c_leiden = algorithms.leiden(gi)
#             a = c_leiden.communities
            list_no_comm.append(size)
            list_mod.append(mod)
        print('Mod: ', list_mod)
        print('No. Comm: ', list_no_comm)
        d_mod[dk] = np.array(list_mod)
        d_no_comm[dk] = np.array(list_no_comm)
    return d_no_comm, d_mod

def gen_graph_collection(d_k_in_DIVs, d_dff, d_nspIDs, d_com, dff_cut=0.05, norm_corr = False, ccut=0.0, binary=False):
    d_g = {}
    for dk in d_k_in_DIVs.keys():
        print('DIV: ', dk)
        k_in_dk = d_k_in_DIVs[dk]
        for k in k_in_dk:
            tmp_dat = np.copy(d_dff[k])
            nsp_ids = d_nspIDs[k]
            nsp_ids_uniq = np.unique(nsp_ids)
            if (0 in nsp_ids_uniq):
                tmp_dat_filt = tmp_dat[nsp_ids>0,:]
                nsp_ids_filt = nsp_ids[nsp_ids>0]
            else:
                tmp_dat_filt = np.copy(tmp_dat)
                nsp_ids_filt = np.copy(nsp_ids)
            tmp_dat_filt[tmp_dat_filt<dff_cut] = 0
            r = np.corrcoef(tmp_dat_filt)
            r = np.nan_to_num(r)
            np.fill_diagonal(r,0)            
            coms = d_com[k]
            coms = coms[nsp_ids>0,:]
            index = list(range(1,coms.shape[0]+1))
            g=nx.Graph()
            for i in range(len(index)):
                g.add_node(index[i])
            for i in range(len(index)):
                for j in range(i+1, len(index)):
                    if (r[i,j] >= ccut):
                        if(not binary):
                            g.add_edge(index[i],index[j],weight=r[i,j])
                        if(binary):
                            g.add_edge(index[i],index[j])
            d_g[k] = g
    return d_g

def calc_global_eff_worker(k,g):
    try:
        lf = nx.global_efficiency(g)
    except:
        lf = np.nan
    print('Global Eff key: ',k, lf)
    return(k,lf)

def convert_pool_dict(pres_prop):
    d_g_prop_to_dict = {}
    for i in range(len(pres_prop)):
        print(pres_prop[i][0], pres_prop[i][1])
        d_g_prop_to_dict[pres_prop[i][0]] = pres_prop[i][1]
    return d_g_prop_to_dict

def extract_intra_inter_nsp_neurons(c_mat,nsp_ids):
    c = np.copy(c_mat)
    x = list(range(0,c.shape[0]))
    ids_uniq = np.unique(nsp_ids)
    c_mask = np.zeros_like(c)
    for id in ids_uniq:
        m = nsp_ids == id
        x_c = [x[i] for i in range(len(x)) if m[i]]
        for i in x_c:
            for j in x_c:
                c_mask[i,j] = 1
    c_1=c_mask[np.triu_indices(c_mask.shape[0], k = 1)]
    c_2=c[np.triu_indices(c.shape[0], k = 1)]
    c_inter = c_2[c_1 < 1]
    c_intra = c_2[c_1 == 1]
    return c_intra, c_inter

# Callback invoked when the StartEvent happens, sets up our new data.
def start_plot():
    global metric_values, multires_iterations
    
    metric_values = []
    multires_iterations = []
    
# Callback invoked when the EndEvent happens, do cleanup of data and figure.
def end_plot():
    global metric_values, multires_iterations
    
    del metric_values
    del multires_iterations
    # Close figure, we don't want to get a duplicate of the plot latter on.
    plt.close()

# Callback invoked when the IterationEvent happens, update our data and display new figure.    
def plot_values(registration_method):
    global metric_values, multires_iterations
    
    metric_values.append(registration_method.GetMetricValue())                                       
    # Clear the output area (wait=True, to reduce flickering), and plot current data
    clear_output(wait=True)
    # Plot the similarity metric values
    plt.plot(metric_values, 'r')
    plt.plot(multires_iterations, [metric_values[index] for index in multires_iterations], 'b*')
    plt.xlabel('Iteration Number',fontsize=12)
    plt.ylabel('Metric Value',fontsize=12)
    plt.show()
    
# Callback invoked when the sitkMultiResolutionIterationEvent happens, update the index into the 
# metric_values list. 
def update_multires_iterations():
    global metric_values, multires_iterations
    multires_iterations.append(len(metric_values))
    
def align_2d(fixed_image, moving_image, hist_bins=50, sampling_frac=0.8, aff=False, iters=50):
	if (aff):
		initial_transform = sitk.CenteredTransformInitializer(fixed_image,moving_image,sitk.AffineTransform(2),sitk.CenteredTransformInitializerFilter.GEOMETRY)
		print('Using Affine')
	else:
		initial_transform = sitk.CenteredTransformInitializer(fixed_image,moving_image,sitk.Euler2DTransform(),sitk.CenteredTransformInitializerFilter.GEOMETRY)
		print('Using Euler 2D')
	registration_method = sitk.ImageRegistrationMethod()
	# Similarity metric settings.
	registration_method.SetMetricAsMattesMutualInformation(numberOfHistogramBins=hist_bins)
	registration_method.SetMetricSamplingStrategy(registration_method.RANDOM)
	registration_method.SetMetricSamplingPercentage(sampling_frac)
	registration_method.SetInterpolator(sitk.sitkLinear)

	# Optimizer settings.
	registration_method.SetOptimizerAsRegularStepGradientDescent(1.0,.001,iters)
	registration_method.SetOptimizerScalesFromPhysicalShift()

	# Setup for the multi-resolution framework.            
	registration_method.SetShrinkFactorsPerLevel(shrinkFactors = [2,1])
	registration_method.SetSmoothingSigmasPerLevel(smoothingSigmas=[2,1])
	registration_method.SmoothingSigmasAreSpecifiedInPhysicalUnitsOn()

	# Don't optimize in-place, we would possibly like to run this cell multiple times.
	registration_method.SetInitialTransform(initial_transform, inPlace=False)

	# Connect all of the observers so that we can perform plotting during registration.
	registration_method.AddCommand(sitk.sitkStartEvent, start_plot)
	registration_method.AddCommand(sitk.sitkEndEvent, end_plot)
	registration_method.AddCommand(sitk.sitkMultiResolutionIterationEvent, update_multires_iterations) 
	#registration_method.AddCommand(sitk.sitkIterationEvent, lambda: plot_values(registration_method))

	final_transform = registration_method.Execute(sitk.Cast(fixed_image, sitk.sitkFloat32), 
												  sitk.Cast(moving_image, sitk.sitkFloat32))
	print('Final metric value: {0}'.format(registration_method.GetMetricValue()))
	print('Optimizer\'s stopping condition, {0}'.format(registration_method.GetOptimizerStopConditionDescription()))
	return final_transform




def convert_f_to_dff_perc(f_NNet_NSP_mat, perc, win_sz=500):
    #win_sz = 250    
    f_base = np.copy(f_NNet_NSP_mat)
    for j in range(f_base.shape[0]):
        t = np.copy(f_NNet_NSP_mat[j,:])
        f_base[j,:] = ndi.percentile_filter(t, perc, win_sz)
    dff_NNet_NSP_mat = np.divide((f_NNet_NSP_mat - f_base),f_base)
    dff_NNet_NSP_mat[dff_NNet_NSP_mat<0] = 0
    return dff_NNet_NSP_mat

