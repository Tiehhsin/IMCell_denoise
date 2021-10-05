#!/bin/bash

#### preprocess image data
python /src/preprocess.py

#### protein quantification
python /src/cell_info.py

#### permutation test
python /src/permutation_test.py --noise_level 0.05 --times 10

#### calculate FDR and get denoise result
python /src/FDR.py --fdr_value 0.1

#### get normalize result
python /src/batch_normalize.py --fdr_value 0.1