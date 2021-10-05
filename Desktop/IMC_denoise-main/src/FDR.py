from multiprocessing import Pool
import os
import time
import numpy as np
import pandas as pd
import cv2
from scipy import stats
import argparse
from path import *

def get_FDR(roi):
    """
        Calculate the FDR and determine the threshold
        :param roi: sample
        :return: FDR and positive cell threshold result
        """
    # create output folder for each roi
    FDR_out_path = os.path.join(FDR_path, roi)
    if not os.path.exists(FDR_out_path):
        os.makedirs(FDR_out_path)

    for marker in list(marker_info.index):
        if marker in permutation_marker:
            permutation_cell = pd.read_csv(os.path.join(permutation_path, roi, roi + '_permutation_result.csv'),
                                           index_col=0)[marker]
            raw_protein = pd.read_csv(os.path.join(protein_path, roi + '.csv'), index_col='ObjectNumber')
            image_cell = raw_protein[raw_protein['Area'] > 5][marker]

            # Gaussian fitting
            (mu1, sigma1) = stats.norm.fit(permutation_cell)
            (mu2, sigma2) = stats.norm.fit(image_cell)
            x_dummy = np.linspace(mu1, max(permutation_cell), 10000)

            FP_list = []
            TP_list = []
            FDR_list = []
            for thres in x_dummy[::-1]:
                FP = 1 - stats.norm(mu1, sigma1).cdf(thres)
                TP = 1 - stats.norm(mu2, sigma2).cdf(thres)
                if (FP == 0) & (TP == 0):
                    FDR = 0
                else:
                    FDR = FP / (FP + TP)
                FP_list.append(FP)
                TP_list.append(TP)
                FDR_list.append(FDR)
            FDR_pd = pd.DataFrame(list(x_dummy[::-1]), columns=['threshold'])
            FDR_pd['FP'] = FP_list
            FDR_pd['TP'] = TP_list
            FDR_pd['FDR'] = FDR_list
            FDR_pd.to_csv(os.path.join(FDR_out_path, roi + '_' + marker + fdr_suffix))


def get_denoise_result(roi, negative_zero=False):
    """
    Get denoise protein expression
    :param roi: sample
    :param negative_zero: negative cell set zero or not
    :return: denoise protein expression
    """
    raw_protein = pd.read_csv(os.path.join(protein_path, roi + '.csv'), index_col='ObjectNumber')
    denoise_protein = pd.DataFrame(list(raw_protein.index), columns=['ObjectNumber'])
    denoise_protein['Position'] = list(raw_protein['Position'])
    denoise_protein['Area'] = list(raw_protein['Area'])

    for marker in list(marker_info.index):
        if marker in permutation_marker:
            background_noise = np.mean(
                pd.read_csv(os.path.join(permutation_path, roi, roi + permutation_suffix), index_col=0)[marker])

            fdr_res = pd.read_csv(os.path.join(FDR_path, roi, roi + '_' + marker + fdr_suffix), index_col=0)
            fdr_list = list(fdr_res['FDR'])
            thre_list = list(fdr_res['threshold'])
            for fdr_ind in range(len(fdr_list)):
                if (fdr_list[fdr_ind] > fdr_threshold.loc[marker]['value']) & (fdr_list[fdr_ind] != 0):
                    thre_value = thre_list[fdr_ind]
                    break

            if negative_zero:
                positive_cell = raw_protein[(raw_protein['Area'] > 5) & (raw_protein[marker] >= thre_value)]
                positive_cell_num = positive_cell.index
                denoise_protein[marker] = [
                    max(0, raw_protein.loc[ind][marker] - background_noise) if ind in positive_cell_num else 0 for ind
                    in raw_protein.index]
            else:
                denoise_protein[marker] = [max(0, raw_protein.loc[ind][marker] - background_noise) for ind in
                                           raw_protein.index]
    denoise_protein.to_csv(os.path.join(output_path, roi + denoise_suffix))
    print('{} finished'.format(roi))


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
    if fdr_file:
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
    output_path = denoise_path
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    multiprocess = True
    if multiprocess:
        multi_threads = 70
        with Pool(multi_threads) as p:
            p.map(get_FDR, [roi for roi in image_name])
            p.map(get_denoise_result, [roi for roi in image_name])
    else:
        for roi in image_name:
            get_FDR(roi)
            get_denoise_result(roi)
