# Denoising of Imaging Mass Cytometry Image for Better Cell Phenotyping

## 1. Introduction
IMC (Imaging Mass Cytometry) is a method to measure different protein expressions in cells. IMC can measure over 100 markers simultaneously while past IF(immunofluorescence) method can only measure several markers. However, in IMC, 1. spillover between different heavy metals and 2. the antibody performance and antigen retrieval condition can be different. These two problem cause the image we get after IMC method can be very noisy.



IMCell in this project is an effective method to denoise the IMC images. The main two parts of the workflow of IMCell is 1. denoise and 2. normalization.

1. denoise

   (1) The raw IMC images be preprocessed by arcsinh function

   The distribution of the raw marker intensities( pixel values in the raw IMC images) can be fairly skew. It's better to transform the skew distribution by arcsinh function to map the distribution in a symmetric range. 

   (2)  The raw IMC images be preprocessed by hotpixel removal

   We can use an easy way to pre-denoise the image at the beginning. If the pixel value is higher than a specific value, replace the pixel value( look our paper for more details).

   (3) Generating decoy cells

   How to quantify the influence of the background noise to the real protein expression? Based on the cell segmentation result by [Dice-XMBD](https://github.com/xmuyulab/Dice-XMBD), we calculated the parameters like major axis, minor axis and orientation of segmented cells.

   Then choose the noise-only region to put decoy cells. The detail of choosing noise-only region is in the paper. We generate ellipses as decoy cells. The parameters of the decoy cells are based on the parameters we got from the real cells in the former step.

   Decoy cells are not supposed to have protein expression, but because of the background noise, they have. We calculated the protein expression of decoy cells in each channel, so we got the quantified influence of background noise.

   (4) FDR method: set the threshold to divide positive cells and negative cells

   - Permutation test to identify positive cells of each channel  
     Some protein channels are very noisy. It seems that the protein is expressed in almost all areas of the tissue image, which may cause many cells to show false positives. The permutation test uses basic statistical methods to find significant protein expression between the cells and the background noise.  
     First, we randomly select an ellipse (circle) area at any position of the ROI, the area of which is within the range of the size of all cells in the ROI, and calculate the protein expression value in this area. This step is repeated a certain number of times to generate a random distribution. Next, we use a one-tailed test to compare the protein expression value of each cell with the random distribution. This test provides us with a p-value which represents whether the cell is positive.  
   - Removal of background noise  
     After identifying the true positive cells through the permutation test, we set the value of the protein expression outside the positive cells to zero to remove background noise. In addition, we subtract the average protein expression value of all negative cells in the ROI from the protein expression value of each positive cell to remove the influence of background noise on positive cells.  

2. Normalization

   To remove the differences between antibody performance and signal-to-noise ratio in the tissue, normalizing the result.

   

## 2. Run the code
1. Download the data you will need:

   (1) [Raw IMC images](https://figshare.com/articles/figure/rawImages/16740379)

   (2) [Cell masks(segmented cells)](https://figshare.com/articles/figure/cellMask/16740664)

   Then put the data in a specific path that you can find.

2. Open terminal and enter into the folder you want to store the project. Then type the command:

   ```
   git clone https://github.com/Tiehhsin/IMCell_denoise.git
   ```

   All the paths used in this project are written in path.py. You need to change the following paths in your own paths.

   - img_path: store the raw IMC images
   - mask_path: store the cell mask (segmented results)
   - project_path: where you store the project

   ```
   #### The paths in my example
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
   ```

   

3. Preprocess the raw IMC images. Including arcsinh transformation and hot pixel removal:

   ```
   python ./preprocess.py
   ```

4. Calculate the parameters of true cells, including major axis, minor axis and orientation angle:

   ```
   python ./cell_info.py
   ```

5. Generate decoy cells and running permutation test:

   ```
   python ./permutation_test.py --noise_level 0.05 --times 10
   ```

6. Positive cells identification by FDR

   ```
   python ./FDR.py --fdr_value 0.1
   ```

7. Normalization

   ```
   python ./batch_normalize.py --fdr_value 0.1
   ```

8. We got the information of each marker and stored them in corresponding .csv files. If the result images is needed, run the following command to get the raw image, denoised image and normalized image. Three functions in plot.py are designed to draw the raw image, denoised image and normalized image.

   ```
   python ./plot_image.py
   ```

   

## References

1. IMCell-XMBD: A statistical approach for robust cell identification and quantification from imaging mass cytometry images. [link]( https://github.com/xmuyulab/IMC_denoise)
2. Dice-XMBD: Deep Learning-Based Cell Segmentation for Imaging Mass Cytometry. [link](https://github.com/xmuyulab/Dice-XMBD)

