%% To convert the Gadgetron-generated .h5 images to nifti .nii images
% Author: Qingping Chen
% Date: 2025.09.23
% Email: qingping.chen@uniklinik-freiburg.de
filename = 'pulseq_mprage_out.h5' ;
info = hdf5info(filename) ;
address_data_1 = info.GroupHierarchy.Groups(1).Groups.Datasets(2).Name ;
im = squeeze(double( hdf5read(filename, address_data_1) ) ) ;
im = permute(im, [3, 2, 1]) ;
im = flip(im, 3) ;
im = flip(im, 1) ;
figure ;
montage(rot90(squeeze(mat2gray(im(:,:,150))))) ;
colormap default ;
niftiwrite(im, 'pulseq_mprage_gt.nii') ;

%% To convert the Gadgetron-reconstructed .h5 EPI image and convert it to nifti images.
clear; close all; clc;
filename = 'pulseq_epi_out.h5' ;
info = hdf5info(filename) ;
address_data_1 = info.GroupHierarchy.Groups(1).Groups.Datasets(2).Name ;
im = squeeze(double( hdf5read(filename, address_data_1) ) ) ;
im = reshape(im, [80,80,48,30]) ;
figure ;
montage(mat2gray(rot90(im(:,:,33,30)))) ;
colormap default ;
niftiwrite(im, 'pulseq_epi_gt.nii') ;