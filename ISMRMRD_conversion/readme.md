# Tutorials for converting Pulseq-generated data into ISMRMRD format
This directory contains materials and tutorials for data format standardization with Pulseq-based MPRAGE and EPI sequences on Siemens, General Electric (GE), Philips, and the United Imaging (UIH) platforms.
## 1. Required toolboxes
Please add the following three toolboxes and their subfolders to MATLAB's path:
### 1.1 [Open-source Pulseq MATLAB toolbox](https://github.com/pulseq/pulseq) to load Pulseq `.seq` files and get LABELs and sequence definitions.
### 1.2 [mapVBVD](https://github.com/pehses/mapVBVD#) toolbox to load the Siemens `.dat` raw data.
### 1.3 [ISMRMRD](https://github.com/ismrmrd/ismrmrd#) toolbox for ISMRMRD conversion.

## 2. MPRAGE
This folder contains materials to convert a Siemens MPRAGE `.dat` raw data acquired from Cima.X 3T to ISMRMRD data.
* Be sure that the Pulseq, mapVBVD, and ISMRMRD toolboxes are in MATLAB's path.
* Run `pulseq2mrd_mprage.m` script. It can convert the Siemens MPRAGE data (`meas_MID00203_FID24417_pulseq_we_mprage_sag_p2_1mm_iso.dat`) to ISMRMRD data (`pulseq_mprage_data.h5`) with the LABELs and sequence definitions loaded from the `mprage_challenge.seq` file.

## 3. EPI
This folder contains materials to convert a Siemens EPI `.dat` raw data from Cima.X 3T to ISMRMRD data.
* Be sure that the Pulseq, mapVBVD, and ISMRMRD toolboxes are in MATLAB's path.
* Run `pulseq2mrd_epi.m` script. It can convert the Siemens EPI data (`meas_MID00207_FID24421_pulseq_epirs_iso_2_8mm_slc48_tran.dat`) to ISMRMRD data (`pulseq_epi_data.h5`) with the LABELs and sequence definitions loaded from the `epi_challenge.seq` file.
