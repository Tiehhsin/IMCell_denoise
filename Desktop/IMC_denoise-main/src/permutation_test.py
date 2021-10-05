from multiprocessing import Pool
import os
import time
import random
import numpy as np
import pandas as pd
import tifffile as tiff
import cv2
import math
import argparse
from utils import median_filter
from path import *

def extract_highconf_area(img, clip_thres=99):
    """
    Indentify high confidence protein expression region and noise region
    :param img: preprocess image
    :param clip_thres: protein expression threshold
    :param noise_level: noise_level
    :return: high confidence protein expression region image
    """
    clip_img = img.copy()
    highconf = np.zeros(img.shape)
    for c in range(img.shape[2]):
        thmin = np.percentile(np.unique(img[:, :, c]), clip_thres) * noise_level
        clip_img[:, :, c][clip_img[:, :, c] < thmin] = 0
        highconf[:, :, c] = median_filter(clip_img[:, :, c])
    return highconf


def permutation_test(roi):
    """
    :param roi: sample
    :return: permutation test result
    """
    t0 = time.time()
    # create output folder for each roi
    out_path = os.path.join(output_path, roi)
    if not os.path.exists(out_path):
        os.makedirs(out_path)

    random_parameter = np.zeros([permutation_times, len(marker_info), 5])
    permutation_res = np.zeros([permutation_times, len(marker_info)])

    # indentify high confidence protein expression region and noise region
    preprocess_img = tiff.imread(os.path.join(preprocess_path, roi + preprocess_suffix))
    if img_format == 'panel_dim0_dim1':
        raw_img = tiff.imread(os.path.join(img_path, roi + img_suffix)).transpose(1, 2, 0)
    elif img_format == 'dim0_dim1_panel':
        raw_img = tiff.imread(os.path.join(img_path, roi + img_suffix))
    highconf = extract_highconf_area(raw_img)
    # tiff.imsave(os.path.join(highconf_path, roi + highconf_suffix), highconf)

    # get cell info
    raw_protein = pd.read_csv(os.path.join(protein_path, roi + '.csv'), index_col='ObjectNumber')
    cell_area = list(raw_protein['Area'])
    cell_rec_min = list(raw_protein['Rec_Min'])
    cell_rec_max = list(raw_protein['Rec_Max'])
    cell_angle = list(raw_protein['Angle'])

    # permutation test
    for index in range(permutation_times):
        for marker in list(marker_info.index):
            if marker in permutation_marker:
                marker_num = marker_info.loc[marker]['marker_number']

                rec_max_ran = 0
                rec_min_ran = 1
                ellipse_area = 0
                test_position_num = 0
                flag = 1
                while ((rec_max_ran < rec_min_ran) | (ellipse_area < min(cell_area)) | (
                        ellipse_area > max(cell_area)) | (test_position_num <= 5) | (flag == 1)):
                    rec_max_ran = random.randrange(min(cell_rec_max), max(cell_rec_max))
                    rec_min_ran = random.randrange(min(cell_rec_min), max(cell_rec_min))
                    ellipse_area = round(math.pi * (rec_max_ran / 2) * (rec_min_ran / 2))

                    cy_ran = random.randrange(1, preprocess_img.shape[0])
                    cx_ran = random.randrange(1, preprocess_img.shape[1])
                    angle_ran = random.uniform(min(cell_angle), max(cell_angle))
                    orientation_ran = random.randint(0, 1)
                    if orientation_ran == 1:
                        angle_ran = 180 - angle_ran

                    test_mask = np.zeros([preprocess_img.shape[0], preprocess_img.shape[1], 3])
                    test_mask = cv2.ellipse(test_mask, (cx_ran, cy_ran), (rec_min_ran // 2, rec_max_ran // 2),
                                            angle_ran, 0, 360, (0, 255, 0), -1)
                    test_position = np.where(test_mask[:, :, 1] == 255)
                    test_position_num = test_position[0].shape[0]

                    flag = 0
                    for ind in range(test_position_num):
                        if highconf[test_position[0][ind], test_position[1][ind], marker_num] != 0:
                            flag = 1
                            break

                pixel_sum = 0
                for ind in range(test_position_num):
                    pixel_sum += preprocess_img[test_position[0][ind], test_position[1][ind], marker_num]

                permutation_res[index, marker_num] = pixel_sum / test_position_num

                random_parameter[index, marker_num, 0] = cy_ran
                random_parameter[index, marker_num, 1] = cx_ran
                random_parameter[index, marker_num, 2] = rec_max_ran
                random_parameter[index, marker_num, 3] = rec_min_ran
                random_parameter[index, marker_num, 4] = angle_ran

    permutation_result = pd.DataFrame(permutation_res, columns=list(marker_info.index))
    # save permutation test result
    permutation_result.to_csv(os.path.join(out_path, roi + permutation_suffix))
    np.save(os.path.join(out_path, roi + '_random_parameter.npy'), random_parameter)
    print('{} finished: {} seconds process time'.format(roi, str(round(time.time() - t0, 2))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--noise_level', type=float, help='noise level of high confidence area extraction',
                        default=0.05)
    parser.add_argument('--times', type=int, help='permutation times', default=1000)
    args = parser.parse_args()
    noise_level = args.noise_level
    permutation_times = args.times

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
    output_path = permutation_path
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    # if not os.path.exists(highconf_path):
    #     os.makedirs(highconf_path)

    multiprocess = True
    if multiprocess:
        multi_threads = 70
        with Pool(multi_threads) as p:
            p.map(permutation_test, [roi for roi in image_name])
    else:
        for roi in image_name:
            permutation_test(roi)
