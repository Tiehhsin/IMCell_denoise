from multiprocessing import Pool
import os
import time
import numpy as np
import pandas as pd
import tifffile as tiff
import cv2
from skimage import measure
import math
from utils import matchsize
from path import *

"""
   Get cell protein expression and cell information(position, area, major axis, minor axis)
"""

def get_cell_info(roi, measure_label=True):
    """
    :param roi: sample name
    :param measure_label: use skimage.measure function on mask or not
    :return: cell protein expression and cell information
    """
    t0 = time.time()
    preprocess_img = tiff.imread(os.path.join(preprocess_path, roi + preprocess_suffix))
    raw_mask = tiff.imread(os.path.join(mask_path, roi + mask_suffix))
    mask = matchsize(preprocess_img, raw_mask)

    if measure_label:
        labels = measure.label(mask, connectivity=2)
    else:
        labels = mask

    props = measure.regionprops(labels)

    labels_index = list(np.unique(labels))
    labels_index.pop(0)

    cell_protein_expression = []
    cell_pos = []
    cell_area = []
    cell_rec_min = []
    cell_rec_max = []
    cell_angle = []

    for i in labels_index:
        imask = mask.copy()
        imask[labels != i] = 0
        imask[labels == i] = 255
        n_pixel = len(imask[imask == 255])
        imasked = cv2.add(preprocess_img, np.zeros(preprocess_img.shape), mask=imask.astype(np.uint8))
        protein_avg = np.sum(imasked, axis=(0, 1))/n_pixel

        cell_protein_expression.append(protein_avg)
        cell_pos.append(props[i - 1].centroid)
        cell_area.append(props[i - 1].area)
        cell_rec_min.append(
            min((props[i - 1].bbox[2] - props[i - 1].bbox[0]), (props[i - 1].bbox[3] - props[i - 1].bbox[1])))
        cell_rec_max.append(
            max((props[i - 1].bbox[2] - props[i - 1].bbox[0]), (props[i - 1].bbox[3] - props[i - 1].bbox[1])))
        cell_angle.append(180 * math.atan(
            (props[i - 1].bbox[3] - props[i - 1].bbox[1]) / (props[i - 1].bbox[2] - props[i - 1].bbox[0])) / math.pi)

    protein_expression_matrix = np.matrix(cell_protein_expression)
    cell_info = pd.DataFrame(protein_expression_matrix, columns=list(marker_info.index))
    cell_info['ObjectNumber'] = labels_index
    cell_info['Position'] = cell_pos
    cell_info['Area'] = cell_area
    cell_info['Rec_Min'] = cell_rec_min
    cell_info['Rec_Max'] = cell_rec_max
    cell_info['Angle'] = cell_angle

    cell_info = cell_info.set_index('ObjectNumber')
    cell_info.to_csv(os.path.join(output_path, roi + '.csv'))
    print('{} finished: {} seconds process time'.format(roi, str(round(time.time() - t0, 2))))


if __name__ == '__main__':
    # get marker information
    marker_info = pd.read_csv(marker_info_path, index_col='marker')

    # get roi_list
    images = np.sort(os.listdir(img_path))
    masks = np.sort(os.listdir(mask_path))
    masks = [i for i in masks if 'cell_mask' in i]
    image_name = [i.split(mask_suffix)[0] for i in masks]
    image_name = [i for i in image_name if 'test' not in i]

    # create output folder
    output_path = protein_path
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    multiprocess = True
    if multiprocess:
        multi_threads = 70
        with Pool(multi_threads) as p:
            p.map(get_cell_info, [roi for roi in image_name])
    else:
        for roi in image_name:
            get_cell_info(roi)