import os
import argparse

parser = argparse.ArgumentParser()
# Dataset
marker_info_path = '/mnt/md0/wsi/davinci/temp/snf_data/IMC/denoise_zhang/denoise_data/Melanoma_panel.csv'
img_path = '/mnt/md0/wsi/davinci/temp/snf_data/IMC/denoise_zhang/analysis/fullstacks'
mask_path = '/mnt/md0/wsi/davinci/temp/snf_data/IMC/denoise_zhang/analysis/predictionProbability/BRCA1_threshold-99.7_withAugnoise-0.5/ori_prop'
img_format ='panel_dim0_dim1'
img_suffix = '_fullstacks.tiff'
mask_suffix = '_pred_Probabilities_cell_mask.tiff'

# Project
# parameter
parser.add_argument('--fdr_value', type=float, help='all markers have the same FDR threshold', default=0.1)
parser.add_argument('--fdr_file', type=str, help='each marker has a different FDR value')
# folder
project_path = '/mnt/md0/wsi/davinci/temp/snf_data/IMC/denoise_zhang/denoise_data'
preprocess_path = os.path.join(project_path, 'preprocess')
protein_path = os.path.join(project_path, 'protein')
highconf_path = os.path.join(project_path, 'highconf')
permutation_path = os.path.join(project_path, 'permutation')
FDR_path = os.path.join(project_path, 'FDR')
denoise_path = os.path.join(project_path, 'denoise_result')
normalize_path = os.path.join(project_path, 'normalize_result')
save_img_path = os.path.join(project_path, 'save_img')

preprocess_suffix = '_preprocess.tiff'
highconf_suffix = '_highconf.tiff'
permutation_suffix = '_permutation_result.csv'
fdr_suffix = '_FDR.csv'
denoise_suffix = '_denoise_result.csv'
normalize_suffix = '_normalize_result.csv'