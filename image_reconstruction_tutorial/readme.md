# Tutorials for Harmonized Image Reconstruction Using Gadgetron
This directory contains materials and tutorials for 
* image reconstruction using Gadgetron for Pulseq-generated ISMRMRD data from Siemens, General Electric (GE), Philips, and the United Imaging (UIH) platforms, and
* vendor-based online image reconstruction for Pulseq-based MPRAGE and EPI sequences on the Siemens platform.
## 1. Instructions for Gadgetron offline reconstruction
### 1.1 Gadgetron installation (on a Windows 11 system)
* Download and install the [Docker](https://www.docker.com/) software in your computer. You may need to install/update the Windows Sub Linux (WSL) system for the Docker installation.
* Open your terminal (PowerShell with an administrative account in Windows) and navigate to the folder you would like to map to the Gadgetron Docker container.
* Run `docker run -t --name gt_latest --detach --volume ${pwd}:/opt/data ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_nocuda:latest`. If Docker is not recognized, set Docker to connect to `C:\Program Files\Docker\Docker\resources\bin` in the Environment Path in Windows. This will download and then launch the [latest Gadgetron version](https://gadgetron.readthedocs.io/en/latest/building.html) in a Docker container. It will also mount your current folder as a data folder inside the container. The installed Gadgetron software within the Docker container is as follows:
<img width="1906" height="414" alt="image" src="https://github.com/user-attachments/assets/2c8db889-7739-45ae-b97a-71100d6bd158" />     

* Start the Gadgetron software within the container by clicking the Triangle below "Actions" in the image above.      

* Run this command `docker exec -ti gt_latest /bin/bash` in the PowerShell (as below). This will execute your Gadgetron container.    
 <img width="658" height="60" alt="image" src="https://github.com/user-attachments/assets/86b6d7e1-da37-4cda-a430-474d70b83b33" />

### 1.2 Data preparation
* Place your Pulseq-generated MPRAGE (`pulseq_mprage_data.h5`) and EPI (`pulseq_epi_data.h5`) ISMRMRD data in the mounted folder.
* Run the command in Terminal: `cd /opt/data` to enter the mounted folder.
* <img width="1830" height="388" alt="image" src="https://github.com/user-attachments/assets/c5644b5f-67b6-4c8a-8c40-1ecdba6dd5a3" />
### 1.3 Gadgetron reconstruction
* For MPRAGE reconstruction, run `gadgetron_ismrmrd_client -f pulseq_mprage_data.h5 -c Generic_Cartesian_Grappa.xml -o pulseq_mprage_out.h5` in the PowerShell to generate the Gadgetron-reconstructed MPRAGE image `pulseq_mprage_out.h5`. `pulseq_mprage_data.h5` is the Pulseq-generated MPRAGE ISMRMRD data.
<img width="1640" height="271" alt="image" src="https://github.com/user-attachments/assets/32706c0e-9f27-4aeb-b7ce-644086755dd3" />
* For EPI reconstruction, run `gadgetron_ismrmrd_client -f pulseq_epi_data.h5 -c epi.xml -o pulseq_epi_out.h5` in the PowerShell to generate the Gadgetron-reconstructed EPI image: `pulseq_epi_out.h5`. `pulseq_epi_data.h5` is the Pulseq-generated EPI ISMRMRD data.     
<img width="1333" height="276" alt="image" src="https://github.com/user-attachments/assets/b6286d76-a03f-47b8-ad45-a40706fca4dc" />

### 1.4 Load Gadgetron-reconstructed images
* Run `load_mprage_gt_recon.m` script to load Gadgetron-reconstructed MPRAGE images (`pulseq_mprage_out.h5`) and save as nifti images (`pulseq_mprage_gt.nii`).
* Run `load_epi_gt_recon.m` script to load Gadgetron-reconstructed EPI images (`pulseq_epi_out.h5`) and save as nifti images (`pulseq_epi_gt.nii`).

## 2. Instructions for vendor-based online reconstruction on Siemens
### 2.1 ICE reconstruction
Before executing the Pulseq-based sequences on a Siemens platform, you can enable ICE online Reconstruction following the instructions below:
* Navigate to the Special Card, set Data handling to ICE STD for NUMARIS/X (e.g. XA60A and XA61A), and ICE 3D/2D for MPRAGE/EPI for NUMARIS/4 (e.g. VB, VD, and VE), as below:

<img width="945" height="407" alt="image" src="https://github.com/user-attachments/assets/a6b0a2b8-63e7-480c-acdc-3b63c5972799" />

* Select `Sum-of-Square` for coil combination.
* Be sure that the maximal pixel intensity does not violate the intensity threshold of **4096**.
### 2.2 Open Recon reconstruction
* Download the Open Recon SDK from [this link](https://www.magnetom.net/t/open-recon-sdk-downloads/7223), and follow the instructions to install Open Recon on your scanner. Note that currently, Open Recon is only available on XA50 or higher versions. For more information, please visit the [Open Recon online forum](https://www.magnetom.net/c/openrecon/35).
* You can also refer to our ISMRM abstract: *Chen et al., Open-source, Flexible, and Reproducible Workflow for Data Acquisition, Reconstruction, and Post-processing Based on Pulseq and Open Recon, ISMRM, Hawaii, 2025* for more information for the integration of Pulseq with Open Recon.

## 3. Reconstructed images
* The MPRAGE images reconstructed using Gadgetron, ICE, and Open Recon are shown below:

<img width="677" height="254" alt="mprage" src="https://github.com/user-attachments/assets/36839e0a-6109-4535-9805-3a40e1046df0" />

* The EPI images reconstructed using Gadgetron, ICE, and Open Recon are shown below:

<img width="666" height="202" alt="epi" src="https://github.com/user-attachments/assets/9e1cdfda-f27a-44f0-817b-24ec08f405e2" />



 






