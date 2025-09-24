# Pulseq to ISMRMRD Conversion Tutorial
This directory contains materials and tutorials for data format standardization with Pulseq-based MPRAGE and EPI sequences on Siemens, General Electric (GE), Philips, and the United Imaging (UIH) platforms.
## 1. Documents
*  `meas_MID00203_FID24417_pulseq_we_mprage_sag_p2_1mm_iso.dat`: MPRAGE raw data acquired from a phantom on the Siemens Cima.X 3T scanner at Freiburg, Germany.
*  `meas_MID00207_FID24421_pulseq_epirs_iso_2_8mm_slc48_tran.dat`: EPI raw data acquired from a phantom on the Siemens Cima.X 3T scanner at Freiburg, Germany.
* `pulseq2mrd_mprage.m`: MATLAB script to convert the Siemens MPRAGE raw data (`meas_MID00203_FID24417_pulseq_we_mprage_sag_p2_1mm_iso.dat`) to ISMRMRD data (`pulseq_mprage_data.h5`).
* `pulseq2mrd_epi.m`: MATLAB script to convert the Siemens EPI raw data (`meas_MID00207_FID24421_pulseq_epirs_iso_2_8mm_slc48_tran.dat`) to ISMRMRD data (`pulseq_epi_data.h5`).
* `mprage_challenge.seq` and `epi_challenge.seq`: `.seq` files for sequence execution on Siemens platforms. It provides information on LABELs and sequence definitions to support Pulseq to ISMRMRD conversion.
 
## 2. Required toolboxes
Please add the following three toolboxes and their subfolders to MATLAB's path:
### 2.1 [Open-source Pulseq MATLAB toolbox](https://github.com/pulseq/pulseq) to load Pulseq `.seq` files and get LABELs and sequence definitions for ISMRMRD conversion.
### 2.2 [mapVBVD](https://github.com/pehses/mapVBVD#) toolbox to load the Siemens `.dat` raw data.
### 2.3 [ISMRMRD](https://github.com/ismrmrd/ismrmrd#) toolbox for ISMRMRD conversion.

## 3. ISMRMRD conversion for MPRAGE
* Be sure that the Pulseq, mapVBVD, and ISMRMRD toolboxes are in MATLAB's path.
* Run `pulseq2mrd_mprage.m` script. It can convert the Siemens MPRAGE data (`meas_MID00203_FID24417_pulseq_we_mprage_sag_p2_1mm_iso.dat`) to ISMRMRD data (`pulseq_mprage_data.h5`) with the LABELs and sequence definitions loaded from the `mprage_challenge.seq` file.

## 4. ISMRMRD conversion for EPI
* Be sure that the Pulseq, mapVBVD, and ISMRMRD toolboxes are in MATLAB's path.
* Run `pulseq2mrd_epi.m` script. It can convert the Siemens EPI data (`meas_MID00207_FID24421_pulseq_epirs_iso_2_8mm_slc48_tran.dat`) to ISMRMRD data (`pulseq_epi_data.h5`) with the LABELs and sequence definitions loaded from the `epi_challenge.seq` file.
