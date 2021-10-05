import numpy as np
import math

def arc_transform(img):
    """
    :param img: raw image
    :return image after arcsinh transform
    """

    img = np.squeeze(img)
    img_arc = np.zeros(img.shape)
    for h in range(img.shape[0]):
        for w in range(img.shape[1]):
            img_arc[h, w] = math.asinh(img[h, w])
    return img_arc


def matchsize(img, mask):
    """
    Match img and mask size
    :param img: raw IMC image
    :param mask: raw mask
    :return: mask of the same size as the image
    """
    x1, y1, z1 = img.shape
    mask = mask[:x1, :y1]
    return mask


def median_filter(img, fliter_win=5):
    """
    :param img: input image
    :param fliter_win: median filter window size
    :return: output image
    """
    img = np.squeeze(img)
    height, width = img.shape
    # Zero padding
    pad = fliter_win // 2
    img_pad = np.zeros((height + pad * 2, width + pad * 2), dtype=np.float)
    img_pad[pad: pad + height, pad: pad + width] = img.copy().astype(np.float)
    img_median = img.copy().astype(np.float)
    # median_filter
    for h in range(pad, height + pad):
        for w in range(pad, width + pad):
            median_value = np.median(img_pad[(h - pad):(h + pad + 1), (w - pad):(w + pad + 1)])
            img_median[h-pad, w-pad] = median_value
    return img_median