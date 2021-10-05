from multiprocessing import Pool
import os
import time
import numpy as np
import pandas as pd
import tifffile as tiff
import cv2
from skimage import measure
import argparse
from utils import matchsize
from path import *


def get_positive_mean_cell(fdr_threshold):
    max_mean = np.zeros([len(image_name), len(marker_info.index)])
    for index in range(len(image_name)):
        roi = image_name[index]

        for marker in marker_info.index:
            if marker in permutation_marker:
                marker_num = marker_info.loc[marker]['marker_number']

                fdr_res = pd.read_csv(os.path.join(FDR_path, roi, roi + '_' + marker + '_FDR.csv'), index_col=0)
                fdr_list = list(fdr_res['FDR'])
                thre_list = list(fdr_res['threshold'])
                for fdr_ind in range(len(fdr_list)):
                    if (fdr_list[fdr_ind] > fdr_threshold.loc[marker]['value']) & (fdr_list[fdr_ind] != 0):
                        thre_value = thre_list[fdr_ind]
                        break

                raw_protein = pd.read_csv(os.path.join(protein_path, roi + '.csv'), index_col='ObjectNumber')
                positive_cell = raw_protein[(raw_protein['Area'] > 5) & (raw_protein[marker] >= thre_value)]

                # remove background noise
                background_noise = np.mean(
                    pd.read_csv(os.path.join(permutation_path, roi, roi + permutation_suffix), index_col=0)[marker])
                positive_remove = [max(0, (v - background_noise)) for v in positive_cell[marker]]

                # mean(positive_cell_remove_noise)
                if len(positive_cell) != 0:
                    max_mean[index, marker_num] = np.mean(positive_remove)
                else:
                    max_mean[index, marker_num] = 0
    max_mean_pd = pd.DataFrame(max_mean, columns=list(marker_info.index))
    max_mean_pd['roi'] = image_name
    max_mean_pd = max_mean_pd.set_index("roi")
    max_mean_pd.to_csv(os.path.join(project_path, 'roi_mean_cell.csv'))
    return max_mean_pd


def get_positive_mean_pixel(fdr_threshold):
    max_mean = np.zeros([len(image_name), len(marker_info.index)])
    for index in range(len(image_name)):
        roi = image_name[index]
        preprocess_img = tiff.imread(os.path.join(preprocess_path, roi + preprocess_suffix))
        raw_mask = tiff.imread(os.path.join(mask_path, roi + mask_suffix))
        mask = matchsize(preprocess_img, raw_mask)
        background_noise = pd.read_csv(os.path.join(permutation_path, roi, roi + permutation_suffix), index_col=0)

        # remove background noise
        img_remove = np.zeros([preprocess_img.shape[0], preprocess_img.shape[1], len(marker_info.index)])
        for marker in marker_info.index:
            if marker in permutation_marker:
                marker_num = marker_info.loc[marker]['marker_number']
                background_noise_marker = np.mean(background_noise[marker])
                for h in range(preprocess_img.shape[0]):
                    for w in range(preprocess_img.shape[1]):
                        img_remove[h, w, marker_num] = max(0,
                                                           (preprocess_img[h, w, marker_num] - background_noise_marker))
        # remove background noise protein_expression
        protein_expression_remove = []
        labels = measure.label(mask, connectivity=2)
        labels_index = list(np.unique(labels))
        labels_index.pop(0)
        for i in labels_index:
            imask = mask.copy()
            imask[labels != i] = 0
            imask[labels == i] = 255
            n_pixel = len(imask[imask == 255])
            imasked = cv2.add(img_remove, np.zeros(np.shape(img_remove)), mask=imask.astype(np.uint8))
            protein_avg = np.sum(imasked, axis=(0, 1)) / n_pixel
            protein_expression_remove.append(protein_avg)
        protein_expression_matrix = np.matrix(protein_expression_remove)
        protein_expression = pd.DataFrame(protein_expression_matrix, columns=list(marker_info.index))
        protein_expression['ObjectNumber'] = labels_index

        for marker in marker_info.index:
            if marker in permutation_marker:
                marker_num = marker_info.loc[marker]['marker_number']

                fdr_res = pd.read_csv(os.path.join(FDR_path, roi, roi + '_' + marker + '_FDR.csv'), index_col=0)
                fdr_list = list(fdr_res['FDR'])
                thre_list = list(fdr_res['threshold'])
                for fdr_ind in range(len(fdr_list)):
                    if (fdr_list[fdr_ind] > fdr_threshold.loc[marker]['value']) & (fdr_list[fdr_ind] != 0):
                        thre_value = thre_list[fdr_ind]
                        break

                raw_protein = pd.read_csv(os.path.join(protein_path, roi + '.csv'), index_col='ObjectNumber')
                positive_cell = raw_protein[(raw_protein['Area'] > 5) & (raw_protein[marker] >= thre_value)]
                positive_remove = protein_expression[protein_expression['ObjectNumber'].isin(positive_cell.index)][
                    marker]

                # mean(positive_cell_remove_noise)
                if len(positive_cell) != 0:
                    max_mean[index, marker_num] = np.mean(positive_remove)
                else:
                    max_mean[index, marker_num] = 0
    max_mean_pd = pd.DataFrame(max_mean, columns=list(marker_info.index))
    max_mean_pd['roi'] = image_name
    max_mean_pd = max_mean_pd.set_index("roi")
    max_mean_pd.to_csv(os.path.join(project_path, 'roi_mean_pixel.csv'))
    return max_mean_pd


def get_normalize_result(roi):
    t0 = time.time()

    raw_protein = pd.read_csv(os.path.join(protein_path, roi + '.csv'), index_col='ObjectNumber')
    normalize_protein = pd.DataFrame(list(raw_protein.index), columns=['ObjectNumber'])
    normalize_protein['Position'] = list(raw_protein['Position'])
    normalize_protein['Area'] = list(raw_protein['Area'])

    preprocess_img = tiff.imread(os.path.join(preprocess_path, roi + preprocess_suffix))
    raw_mask = tiff.imread(os.path.join(mask_path, roi + mask_suffix))
    mask = matchsize(preprocess_img, raw_mask)
    background_noise = pd.read_csv(os.path.join(permutation_path, roi, roi + permutation_suffix), index_col=0)
    max_mean = pd.read_csv(os.path.join(project_path, 'roi_mean_cell.csv'), index_col='roi')
    # max_mean = pd.read_csv(os.path.join(project_path, 'roi_mean_pixel.csv'), index_col='roi')

    # remove background noise and scale
    img_scale = np.zeros([preprocess_img.shape[0], preprocess_img.shape[1], len(marker_info.index)])
    for marker in marker_info.index:
        if marker in permutation_marker:
            max_mean_marker = np.max(max_mean[marker])
            current_mean_marker = max_mean.loc[roi][marker]

            marker_num = marker_info.loc[marker]['marker_number']
            background_noise_marker = np.mean(background_noise[marker])
            for h in range(preprocess_img.shape[0]):
                for w in range(preprocess_img.shape[1]):
                    if current_mean_marker == 0:
                        img_scale[h, w, marker_num] = max(0,
                                                          (preprocess_img[h, w, marker_num] - background_noise_marker))
                    else:
                        img_scale[h, w, marker_num] = max(0, (
                                    preprocess_img[h, w, marker_num] - background_noise_marker)) * (
                                                                  max_mean_marker / current_mean_marker)
    protein_expression_scale = []
    labels = measure.label(mask, connectivity=2)
    labels_index = list(np.unique(labels))
    labels_index.pop(0)
    for i in labels_index:
        imask = mask.copy()
        imask[labels != i] = 0
        imask[labels == i] = 255
        n_pixel = len(imask[imask == 255])
        imasked = cv2.add(img_scale, np.zeros(np.shape(img_scale)), mask=imask.astype(np.uint8))
        protein_avg = np.sum(imasked, axis=(0, 1)) / n_pixel
        protein_expression_scale.append(protein_avg)
    protein_expression_matrix = np.matrix(protein_expression_scale)
    for marker in marker_info.index:
        marker_num = marker_info.loc[marker]['marker_number']
        normalize_protein[marker] = protein_expression_matrix[:, marker_num]

    normalize_protein.to_csv(os.path.join(output_path, roi + normalize_suffix))
    print('{} finished: {} seconds process time'.format(roi, str(round(time.time() - t0, 2))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fdr_value', type=float, help='all markers have the same FDR threshold', default=0.1)
    parser.add_argument('--fdr_file', type=str, help='each marker has a different FDR value')
    args = parser.parse_args()
    fdr_value = args.fdr_value
    if args.fdr_file:
        fdr_file = args.fdr_file
    else:
        fdr_file = None

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

    # get fdr threshold
    if args.fdr_file:
        fdr_threshold = pd.read_csv(fdr_file, index_col='marker')
    else:
        fdr_threshold = pd.DataFrame(permutation_marker, columns=['marker'])
        fdr_threshold['value'] = [args.fdr_value] * len(permutation_marker)
        fdr_threshold = fdr_threshold.set_index('marker')

    # get roi_list
    images = np.sort(os.listdir(img_path))
    masks = np.sort(os.listdir(mask_path))
    masks = [i for i in masks if 'cell_mask' in i]
    image_name = [i.split(mask_suffix)[0] for i in masks]
    image_name = [i for i in image_name if 'test' not in i]

    # create output folder
    output_path = normalize_path
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    max_mean = get_positive_mean_cell(fdr_threshold)
    # max_mean = get_positive_mean_pixel(fdr_threshold)
    multiprocess = True
    if multiprocess:
        multi_threads = 70
        with Pool(multi_threads) as p:
            p.map(get_normalize_result, [roi for roi in image_name])
    else:
        for roi in image_name:
            get_normalize_result(roi)