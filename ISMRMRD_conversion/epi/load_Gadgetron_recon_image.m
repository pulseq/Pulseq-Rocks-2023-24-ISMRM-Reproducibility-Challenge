%% To load the Gadgetron-reconstructed .h5 EPI image and convert it to nifti images.

clear; close all; clc;
filename = 'pulseq_epi_out.h5' ;
info = hdf5info(filename) ;
address_data_1 = info.GroupHierarchy.Groups(1).Groups.Datasets(2).Name ;
im = squeeze(double( hdf5read(filename, address_data_1) ) ) ;
im = reshape(im, [80,80,48,30]) ;
figure ;
montage(mat2gray(im(:,:,:,30))) ;
colormap default ;
niftiwrite(im, 'pulseq_epi_gt.nii') ;