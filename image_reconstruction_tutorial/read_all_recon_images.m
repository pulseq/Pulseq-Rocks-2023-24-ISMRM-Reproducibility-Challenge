%% Load all example recon images
% Author: Qingping Chen
% Date: 23.09.2025
% Email: qingping.chen@uniklinik-freiburg.de

% load all example MPRAGE and EPI images reconstructed by ICE, Gadgetron
% and Open Recon. The images were acquired from a phantom on a Cima.X
% scanner at Freiburg, Germany

%% Load MPRAGE images
pulseq_mprage_ice = niftiread('pulseq_mprage_ice.nii') ;
pulseq_mprage_gt = niftiread('pulseq_mprage_gadgetron.nii') ;
pulseq_mprage_or = niftiread('pulseq_mprage_openrecon.nii') ;
figure;
subplot 131 ;
montage(mat2gray(rot90(pulseq_mprage_gt(:,:,150))) ) ;
title('Gadgetron') ;
subplot 132 ;
montage(mat2gray(rot90(pulseq_mprage_ice(:,:,150))) ) ;
title('ICE') ;
subplot 133 ;
montage(mat2gray(rot90(pulseq_mprage_or(:,:,150))) ) ;
title('Open Recon') ;
exportgraphics(gcf,'mprage.png',...
    'ContentType','vector',...
    'BackgroundColor','none') ;

%% Load EPI images
pulseq_epi_ice = niftiread("pulseq_epi_ice.nii") ;
pulseq_epi_gt = niftiread("pulseq_epi_gadgetron.nii") ;
pulseq_epi_or = niftiread('pulseq_epi_openrecon.nii') ;
figure;
subplot 131 ;
montage(mat2gray(rot90(pulseq_epi_gt(:,:,31,end))) ) ;
title('Gadgetron') ;
subplot 132 ;
montage(mat2gray(rot90(pulseq_epi_ice(:,:,31,end))) ) ;
title('ICE') ;
subplot 133 ;
montage(mat2gray(rot90(pulseq_epi_or(:,:,31,end))) ) ;
title('Open Recon') ;
exportgraphics(gcf,'epi.png',...
    'ContentType','vector',...
    'BackgroundColor','none') ;