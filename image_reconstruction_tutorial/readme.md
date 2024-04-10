# Tutorials for Harmonized Image Reconstruction Using Gadgetron
This directory contains materials and tutorials for image reconstruction using vendor-provided online reconstruction and Gadgetron software on Siemens and GE scanners.
## 1. Prerequisites
### 1.1 Gadgetron Docker image
Download a Gadgetron Docker image (Figure 1) and map a folder in your computer to the Gadgetron Docker image (https://gadgetron.github.io/tutorial/).
![dockerimage](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/assets/26165904/5a80d91d-fec7-4bf7-8dd9-583645902fd7)
**Figure 1** gt1 is the Gadgetron docker image.
### 1.2 ISMRMRD software
Please install the ISMRMRD software (https://github.com/ismrmrd/ismrmrd) in your Matlab/Python software for converting Pulseq raw data to ISMRMRD data.
### 1.3 mapVBVD software
The mapVBVD software is required to read raw data in Siemens TWIX format. Please download the software from https://github.com/pehses/mapVBVD install it to your Matlab/Python software.
### 1.4 DICOM converter
Please download the dicm2nii converter (https://github.com/xiangruili/dicm2nii) to convert the DICOM images to .nii images.
## 2. MPRAGE Image Reconstruction
### 2.1 ICE online reconstruction
Pulseq enables ICE online reconstruction for MPRAGE with GRAPPA acceleration on Siemens scanners. To turn on ICE online recon, please turn on the iPAT card and select `Data handling` in the `Special Card` to be `ICE STD` for Numaris/X scanners (`ICE 3D` for Numaris/4 scanners), as below.
![GRAPPA_ICE_recon_setting](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/assets/26165904/7e07f2d0-dbee-4c3e-9d9d-402acafb8f28)
### 2.2 ISMRMRD conversion
* Please be sure that ISMRMRD software has been added to your Matlab path. Please run the `pulseq2mrd_mprage.m` script with the corresponding `.seq` file (for Siemens users, `mprage.seq`; for GE users, `mprage_noiseScan.seq`) in the same folder to convert the Siemens/GE Pulseq raw data to ISMRMRD data. you might need to adjust the `pulseq2mrd_mprage.m` script a bit for your specific processing.
* For the Siemens product raw data, you can use Gadgetron in the docker image to convert it to ISMRMRD data. Please place the raw data in the mapped folder, and run e.g. `siemens_to_ismrmrd -f meas_MID00562_FID17465_b600.dat -o mprage_data.h5` to produce the ISMRMRD data `mprage_data.h5`. Gadgetron does not support Siemens XA data conversion at the moment - you can use the `siemens2mrd_mprage` script for the product MPRAGE data conversion.
### 2.3 Gadgetron offline reconstruction
* Please place your MPRAGE ISMRMRD data `mprage_data.h5` in the mapped folder
* Open terminal (type “cmd” in the start search bar and enter in Windows)
* Enter Gadgetron docker image: `docker exec -ti gt1 /bin/bash`
* Go to the mapped folder: `cd /opt/data`
* Reconstruct MPRAGE raw data: `gadgetron_ismrmrd_client -f mprage_data.h5 -c Generic_Cartesian_Grappa.xml -o mprage_out.h5`, as below:
![gadgetron_recon](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/assets/26165904/c59577a2-c032-45b7-95a2-a90a979b5dfc)
mprage_out.h5 is the reconstructed *.h5* MPRAGE images.
### 2.4 Read Gadgetron-based .h5 recon images
In Matlab, run the following lines to load the reconstructed mprage_out.h5 MPRAGE images.     
`filename = 'mprage_out.h5' ;`        
`info = hdf5info(filename) ;`      
`address_data_1 = info.GroupHierarchy.Groups(1).Groups.Datasets(2).Name ;`          
`mprage_im = squeeze(double( hdf5read(filename, address_data_1) ) ) ;`   
mprage_im is the MPRAGE images reconstructed by Gadgetron.
## 3. EPI Image Reconstruction
### 3.1 ICE online reconstruction
To turn on ICE online reconstruction for EPI data acquisition, be sure to **turn off** the iPAT card and select `Data handling` in the `Special Card` to be `ICE STD` for Numaris/X scanners (`ICE 2D` for Numaris/4 scanners), as below.
![EPI_ICE_recon_setting](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/assets/26165904/3fc8fad7-9d85-4b26-8682-6331260f41b3)
### 3.2 ISMRMRD conversion
* Please be sure that ISMRMRD software has been added to your Matlab path. Please run the `pulseq2mrd_epi.m` script with the corresponding `epi_rs.seq` file in the same folder to convert the Siemens/GE Pulseq raw data to ISMRMRD data. you might need to adjust the `pulseq2mrd_epi.m` script a bit for your specific processing.
* For the Siemens product raw data, you can use Gadgetron in the docker image to convert it to ISMRMRD data. Please place the raw data in the mapped folder, and run e.g. `siemens_to_ismrmrd -f meas_MID00562_FID17465_b600.dat -z 2 -m IsmrmrdParameterMap_Siemens.xml -x IsmrmrdParameterMap_Siemens_EPI.xsl -o epi_data.h5` to produce the ISMRMRD data `epi_data.h5`. Gadgetron does not support Siemens XA data conversion at the moment - you can use the `siemens2mrd_epi` script for the product EPI data conversion.
### 3.3 Gadgetron offline reconstruction
* Please place your EPI ISMRMRD data `epi_data.h5` in the mapped folder
* Open terminal (type “cmd” in the start search bar and enter in Windows)
* Enter Gadgetron docker image: `docker exec -ti gt1 /bin/bash`
* Go to the mapped folder: `cd /opt/data`
* Reconstruct EPI raw data: `gadgetron_ismrmrd_client -f epi_data.h5 -c epi.xml -o epi_out.h5`, as shown in the Figure in Section 2.3.
epi_out.h5 is the reconstructed *.h5* EPI images.
### 3.4 Read Gadgetron-based .h5 recon images
In Matlab, run the following lines to load the reconstructed epi_out.h5 MPRAGE images.     
`filename = 'epi_out.h5' ;`        
`info = hdf5info(filename) ;`      
`address_data_1 = info.GroupHierarchy.Groups(1).Groups.Datasets(2).Name ;`          
`epi_im = squeeze(double( hdf5read(filename, address_data_1) ) ) ;`    
epi_im is the EPI images reconstructed by Gadgetron.
