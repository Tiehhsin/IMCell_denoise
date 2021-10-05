from multiprocessing import Pool
import sys
import os
import time
import numpy as np
import tifffile as tiff
from utils import arc_transform
from path import *

"""
   Preprocess
   Step1, pixel arcsinh transformed
   Step2, hot_pixel removal
"""

def hotpixel_filter(img, hotpixel_thres0=0.98, hotpixel_thres1=4, fliter_win=5):
    """
    :param img: raw image
    :param hotpixel_thre0: "hot pixel" pixel value threshold
    :param hotpixel_thre1: the multiple of "hot pixel" higher than surrounding points
    :param fliter_win: median filter window
    :return image after hot pixel removal
    """

    img = np.squeeze(img)
    height, width = img.shape
    
    # hotpixel_threshold
    threshold = sorted(np.unique(img))[round(len(np.unique(img)) * hotpixel_thres0)-1]
    
    # Zero padding
    pad = fliter_win // 2
    img_pad = np.zeros((height + pad * 2, width + pad * 2), dtype=np.float)
    img_pad[pad: pad + height, pad: pad + width] = img.copy().astype(np.float)
    img_median = img.copy().astype(np.float)

    # median_filter
    for h in range(pad, height + pad):
        for w in range(pad, width + pad):
            if img_pad[h, w] >= threshold:
                median_value = np.median(img_pad[(h - pad):(h + pad + 1), (w - pad):(w + pad + 1)])
                if img_pad[h, w] >= hotpixel_thres1 * median_value:
                    img_median[h-pad, w-pad] = median_value
    return img_median


def preprocess(roi):
    t0 = time.time()
    if img_format == 'panel_dim0_dim1':
        raw_img = tiff.imread(os.path.join(img_path, roi + img_suffix)).transpose(1, 2, 0)
    elif img_format == 'dim0_dim1_panel':
        raw_img = tiff.imread(os.path.join(img_path, roi + img_suffix))

    img_preprocess = np.zeros(raw_img.shape)
    for marker_num in range(raw_img.shape[2]):
        img_arc = arc_transform(raw_img[:, :, marker_num])
        img_preprocess[:, :, marker_num] = hotpixel_filter(img_arc)
    tiff.imsave(os.path.join(output_path, roi + preprocess_suffix), img_preprocess)
    print('{} finished: {} seconds process time'.format(roi, str(round(time.time() - t0, 2))))


if __name__ == '__main__':
    if (img_format != 'panel_dim0_dim1') & (img_format != 'dim0_dim1_panel'):
        print('IMC image format error.')
        sys.exit(0)

    # get roi_list
    images = np.sort(os.listdir(img_path))
    masks = np.sort(os.listdir(mask_path))
    masks = [i for i in masks if 'cell_mask' in i]
    image_name = [i.split(mask_suffix)[0] for i in masks]
    image_name = [i for i in image_name if 'test' not in i]

    # create output folder
    output_path = preprocess_path
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    multiprocess = True
    if multiprocess:
        multi_threads = 70
        with Pool(multi_threads) as p:
            p.map(preprocess, [roi for roi in image_name])
    else:
        for roi in image_name:
            preprocess(roi)