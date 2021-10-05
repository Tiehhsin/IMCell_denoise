from multiprocessing import Pool
import os
import time
import numpy as np
import pandas as pd
import math
import tifffile as tiff
from skimage import measure
from PIL import Image
from path import *
from utils import matchsize, median_filter

def plot_raw(roi, scale=80, marker_list='all'):
    if img_format == 'panel_dim0_dim1':
        raw_img = tiff.imread(os.path.join(img_path, roi + img_suffix)).transpose(1, 2, 0)
    elif img_format == 'dim0_dim1_panel':
        raw_img = tiff.imread(os.path.join(img_path, roi + img_suffix))

    out_path = os.path.join(output_path, 'raw', roi)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    if marker_list == 'all':
        plot_marker = permutation_marker
    else:
        plot_marker = marker_list
    for marker in plot_marker:
        if marker in permutation_marker:
            marker_num = marker_info.loc[marker]['marker_number']
            img0 = np.squeeze(raw_img[:, :, marker_num])

            img_res_3channel = np.zeros([raw_img.shape[0], raw_img.shape[1], 3])
            for h in range(raw_img.shape[0]):
                for w in range(raw_img.shape[1]):
                    img_res_3channel[h, w, 1] = min(img0[h, w] * scale, 255)
            Image.fromarray(img_res_3channel.astype('uint8')).save(
                os.path.join(out_path, roi + '_' + marker + '_raw.png'))


def plot_denoise(roi, scale=80, marker_list='all', measure_label=True):
    preprocess_img = tiff.imread(os.path.join(preprocess_path, roi + preprocess_suffix))
    if img_format == 'panel_dim0_dim1':
        raw_img = tiff.imread(os.path.join(img_path, roi + img_suffix)).transpose(1, 2, 0)
    elif img_format == 'dim0_dim1_panel':
        raw_img = tiff.imread(os.path.join(img_path, roi + img_suffix))
    raw_mask = tiff.imread(os.path.join(mask_path, roi + mask_suffix))
    mask = matchsize(raw_img, raw_mask)
    if measure_label:
        labels = measure.label(mask, connectivity=2)
    else:
        labels = mask
    labels_index = list(np.unique(labels))
    labels_index.pop(0)

    out_path = os.path.join(output_path, 'denoise', roi)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    if marker_list == 'all':
        plot_marker = permutation_marker
    else:
        plot_marker = marker_list
    for marker in plot_marker:
        if marker in permutation_marker:
            marker_num = marker_info.loc[marker]['marker_number']
            img0 = np.squeeze(preprocess_img[:, :, marker_num])

            background_noise = np.mean(
                pd.read_csv(os.path.join(permutation_path, roi, roi + permutation_suffix), index_col=0)[marker])

            img_res_3channel = np.zeros([preprocess_img.shape[0], preprocess_img.shape[1], 3])
            for h in range(labels.shape[0]):
                for w in range(labels.shape[1]):
                    if labels[h, w] in labels_index:
                        img_res_3channel[h, w, 1] = min(max(img0[h, w] - background_noise, 0) * scale, 255)
                    else:
                        img_res_3channel[h, w, 1] = 0
            Image.fromarray(img_res_3channel.astype('uint8')).save(
                os.path.join(out_path, roi + '_' + marker + 'denoise.png'))


def plot_normalize(roi, scale=80, marker_list='all', measure_label=True, median_post=True):
    preprocess_img = tiff.imread(os.path.join(preprocess_path, roi + preprocess_suffix))
    if img_format == 'panel_dim0_dim1':
        raw_img = tiff.imread(os.path.join(img_path, roi + img_suffix)).transpose(1, 2, 0)
    elif img_format == 'dim0_dim1_panel':
        raw_img = tiff.imread(os.path.join(img_path, roi + img_suffix))
    raw_mask = tiff.imread(os.path.join(mask_path, roi + mask_suffix))
    mask = matchsize(raw_img, raw_mask)
    if measure_label:
        labels = measure.label(mask, connectivity=2)
    else:
        labels = mask
    labels_index = list(np.unique(labels))
    labels_index.pop(0)

    max_mean = pd.read_csv(os.path.join(project_path, 'roi_mean_cell.csv'), index_col='roi')

    out_path = os.path.join(output_path, 'normalize', roi)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    if marker_list == 'all':
        plot_marker = permutation_marker
    else:
        plot_marker = marker_list

    if median_post:
        for marker in plot_marker:
            if marker in permutation_marker:
                max_mean_marker = np.max(max_mean[marker])
                current_mean_marker = max_mean.loc[roi, marker]

                marker_num = marker_info.loc[marker]['marker_number']
                img0 = np.squeeze(preprocess_img[:, :, marker_num])

                background_noise = np.mean(
                    pd.read_csv(os.path.join(permutation_path, roi, roi + permutation_suffix), index_col=0)[marker])

                img_res = np.zeros([preprocess_img.shape[0], preprocess_img.shape[1]])
                img_res_3channel = np.zeros([preprocess_img.shape[0], preprocess_img.shape[1], 3])
                for h in range(labels.shape[0]):
                    for w in range(labels.shape[1]):
                        if labels[h, w] in labels_index:
                            if current_mean_marker == 0:
                                img_res[h, w] = max(img0[h, w] - background_noise, 0)
                            else:
                                img_res[h, w] = max(img0[h, w] - background_noise, 0) * (
                                        max_mean_marker / current_mean_marker)
                        else:
                            img_res[h, w] = 0

                img_res_median = median_filter(img_res)
                for h in range(img_res_median.shape[0]):
                    for w in range(img_res_median.shape[1]):
                        img_res_3channel[h, w, 1] = min(img_res_median[h, w] * scale, 255)
                Image.fromarray(img_res_3channel.astype('uint8')).save(
                    os.path.join(out_path, roi + '_' + marker + '_normalize.png'))

    else:
        for marker in plot_marker:
            if marker in permutation_marker:
                max_mean_marker = np.max(max_mean[marker])
                current_mean_marker = max_mean.loc[roi, marker]

                marker_num = marker_info.loc[marker]['marker_number']
                img0 = np.squeeze(preprocess_img[:, :, marker_num])

                background_noise = np.mean(
                    pd.read_csv(os.path.join(permutation_path, roi, roi + permutation_suffix), index_col=0)[marker])

                img_res_3channel = np.zeros([preprocess_img.shape[0], preprocess_img.shape[1], 3])
                for h in range(labels.shape[0]):
                    for w in range(labels.shape[1]):
                        if labels[h, w] in labels_index:
                            if current_mean_marker == 0:
                                img_res_3channel[h, w, 1] = min(max(img0[h, w] - background_noise, 0) * scale, 255)
                            else:
                                img_res_3channel[h, w, 1] = min(max(img0[h, w] - background_noise, 0) * (
                                        max_mean_marker / current_mean_marker) * scale, 255)
                        else:
                            img_res_3channel[h, w, 1] = 0
                Image.fromarray(img_res_3channel.astype('uint8')).save(
                    os.path.join(out_path, roi + '_' + marker + '_normalize.png'))



if __name__ == '__main__':
    # get marker information
    marker_info = pd.read_csv(marker_info_path, index_col='marker')
    Tmarker = ['CD45', 'CD3', 'CD4', 'CD8a', 'CD27', 'FoxP3', 'CD127', 'CD194', 'CD278']
    Bmarker = ['CD20', 'CD38']
    Mono = ['HLA_DR', 'CD74', 'CD68', 'CD14', 'CD16', 'CD11c', 'CD11b', 'IDO']
    stroma = ['Vimentin', 'Alpha_SMA', 'E_Cadherin', 'EpCAM', 'CA_IX', 'VEGF', 'PDGFRb', 'FAP', 'AXL', 'Collagen_I',
              'Ki67']
    ICB_dic = {'CD134': 'CD134_OX40', 'CD279': 'CD279_PD1', 'CD223': 'CD223_LAG3', 'CD274': 'CD274_PDL1',
               'CD366': 'CD366_TIM3'}
    permutation_marker = Tmarker + Bmarker + Mono + stroma + list(ICB_dic.keys())

    # get roi_list
    images = np.sort(os.listdir(img_path))
    masks = np.sort(os.listdir(mask_path))
    masks = [i for i in masks if 'cell_mask' in i]
    image_name = [i.split(mask_suffix)[0] for i in masks]
    image_name = [i for i in image_name if 'test' not in i]

    # create output folder
    output_path = save_img_path
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    multiprocess = True
    if multiprocess:
        multi_threads = 70
        with Pool(multi_threads) as p:
            p.map(plot_raw, [roi for roi in image_name])
            p.map(plot_denoise, [roi for roi in image_name])
            p.map(plot_normalize, [roi for roi in image_name])
    else:
        for roi in image_name:
            plot_raw(roi)
            plot_denoise(roi)
            plot_normalize(roi)