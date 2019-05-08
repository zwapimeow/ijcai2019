# neurips2019

Resources for paper "A Novel Efficient Approach with Data-Adaptive Capability for OMP-based Sparse Subspace Clustering".


EYALEB_DATA : size = [D,N]
              Each column is a face image and N = 38 subjects x 64 images/subject = 2414
              Each image is downsampled from 192x168 to D = 48x42 = 2016
              
EYALEB_LABEL : size = [1,N] 
               Each entry is the label for the corresponding column in EYALEB_DATA

usps_data.mat : size = [(16x16),(10x1100)] 
                imagesize 16x16
                10 numbers(1-9;0) x 1100 images/number
                
usps_lable.mat : size = [(10x1100),1]

Some codes are finetuned from :

【CVPR 2016 paper】C. You, D. Robinson, R. Vidal, Scalable Sparse Subspace Clustering by Orthogonal Matching Pursuit

【ICIP 2018 paper】L. Zhong, Y. S. Zhu, and G. B. Luo. A New Sparse Subspace Clustering by Rotated Orthogonal Matching Pursuit

Download other useful codes according to instructions in paper or provided codes.

Comparison between SSCOMP&oursSSCOMP or SSCROMP&oursSSCROMP can help researchers understand our approach easily.

Run run_EYALEB and run_usps to get results.

:)
