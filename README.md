# vesicle-quant
Quantify vesicle membranes from negative-stain EM images

Built atop my colleague [Virly Ananda](https://github.com/virlyananda)'s excellent [EM Image Processing scripts](https://github.com/virlyananda/EM-ImageProcessing).  Check out her other cool projects!

This python script will take all `*.tif` files in the current working directory, preprocess them with a gaussian blur and dilation, then use thresholding and background noise reduction to quantify total area of features in the image.  Run with `python preprocess.py` in the working directory containing your tif files.

Example:

![tab_240213_collection_ttyh_002_crop_processing](https://github.com/tribell4310/vesicle-quant/assets/67428134/250c0007-705b-4a09-ae5e-c69d9876a438)
