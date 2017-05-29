Information
===========

nu_corrector is an impementation of our single-image based vignetting or bias correction systems based on the sparsity property of image gradient distribution. We use non-uniformity to denote both vignetting and bias.

"Definition of vignetting and bias":

Vignetting refers to the phenomenon of brightness attenuation away from the image center, and is an artifact that is prevalent in photography. 

Bias of image denotes the spatial variations of intensity/color caused by illumination changes for images taken by a digital camera, by inhomogenious magnetic field for MR images obtained with an MRI machine, or by non-uniform X-ray beam for CT images acquired with a CT scanner.

"Harm of vignetting and bias":

Vignetting and bias can significantly impair computer vision algorithms that rely on precise intensity data. They include photometric methods such as shape from shading, appearance-based techniques such as object recognition and image mosaicing, and many other applications such as image segmentation, image registration, and feature extraction.

"Brief introduction to our work":

In 2006, we proposed probably the first *single-image based* vignetting correction technique [3][4]. In 2008, we designed extremely efficient algorithms for vignetting correction based on sparsity property of image gradient distribution [1]. In 2009, we designed algorithms for estimating bias given a single image. Vignetting is generally assumed to be radially symmetric. In contrast, bias can be a smooth field in any format, which can be represented by for example a bipoly model. nu_corrector is an implementation for [1][2]. The code for [3][4] is not available due to license.

[1]
@InProceedings{Zheng08,
  author = {Y. Zheng, J. Yu, S. B. Kang, S. Lin, C. Kambhamettu},
  title = {Single-image vignetting correction using radial gradient symmetry},
  booktitle = {IEEE Conference on Computer Vision and Pattern Recognition},
  year = {2008},
  month = {June},
  pages = {1--8}
}

[2]
@InProceedings{Zheng09,
  author = {Y. Zheng, M. Grossman, S. Awate, J. Gee },
  title = {Automatic Correction of Intensity Nonuniformity From Sparseness of Gradient Distribution in Medical Images},
  booktitle = {the 12th International Conference on Medical Image Computing and Computer Assisted Intervention},
  year = {2009},
  month = {September},
}

[3]
@InProceedings{Zheng06,
  author = {Y. Zheng, S. Lin, S. Kang},
  title = {Single-Image Vignetting correction},
  booktitle = {IEEE Conference on Computer Vision and Pattern Recognition},
  year = {2006},
  month = {June}
}

[4]
@ARTICLE{,
  AUTHOR =       {Y. Zheng, S. Lin, C. Kambhamettu, J. Yu, S. Kang},
  TITLE =        {Single-Image Vignetting Correction},
  JOURNAL =      {IEEE. Trans. Pattern Analysis and Machine Intelligence},
  YEAR =         {2009},
  volume =       {31},
  pages =        {2243--2256}
}

Basic Usage
===========

The system is implemented in matlab, with one function
written in C++ for efficiency reasons.

Basic steps for running the system:

1. Unpack the code.
2. Start matlab from the directory './code/'.
3. Run the 'demo_vignetting.m' script or 'demo_bias.m' script.

The main functions defined by the code are:

function [im_corrected,im_vignetting]=vignCorrection_nonPara(im_given,varargin)
function [im_corrected,im_bias]=biasCorrection_bipoly(im_given,varargin)

Their usage is demonstrated in the 'demo' script.

By Yuanjie Zheng @ PICSL lab @ UPenn on April 16th, 2010
email: zheng.vision@gmail.com
website: http://sites.google.com/site/zhengvision/