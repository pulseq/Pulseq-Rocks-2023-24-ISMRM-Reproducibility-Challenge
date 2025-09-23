# Image Reconstruction Tutorial
This directory contains materials and tutorials for 
* image reconstruction using Gadgetron for Pulseq-generated ISMRMRD data from Siemens, General Electric (GE), Philips, and the United Imaging (UIH) platforms, and
* vendor-based online image reconstruction for Pulseq-based MPRAGE and EPI sequences on the Siemens platform.
## 1. Documents
The documents contained in this tutorial are listed below:
* `pulseq_mprage_out.h5` and `pulseq_epi_out.h5`: Gadgetron-reconstructed Pulseq-based MPRAGE and EPI `.h5` images. The data were acquired from a phantom on the Cima.X 3T at Freiburg, Germany.
* `convert_h5_to_nifti.m`: MATLAB script to convert `pulseq_mprage_out.h5` and `pulseq_epi_out.h5` to NIFTI images (`.nii`).
* `pulseq_mprage_gadgetron.nii` and `pulseq_epi_gadgetron.nii`: the Pulseq-based MPRAGE and EPI `.nii` images converted by the script `convert_h5_to_nifti.m`.
* `pulseq_mprage_ice.nii` and `pulseq_epi_ice.nii`: the Pulseq-based MRPAGE and EPI NIFTI images reconstructed by Siemens' ICE software.
* `pulseq_mprage_openrecon.nii` and `pulseq_epi_openrecon.nii`: the Pulseq-based MRPAGE and EPI NIFTI images reconstructed by Siemens' Open Recon software.
* `read_all_recon_images.m`: MATLAB script to read all NIFTI images and produce `.png` figures.
## 2. Instructions for Gadgetron offline reconstruction
### 2.1 Gadgetron installation (on a Windows 11 system)
* Download and install the [Docker](https://www.docker.com/) software in your computer. You may need to install/update the Windows Sub Linux (WSL) system for the Docker installation.
* Open your terminal (PowerShell with an administrative account in Windows) and navigate to the folder you would like to map to the Gadgetron Docker container, e.g. by running this command `cd C:\Users\chenq\Downloads\Docker_folder`, as below:

<img width="1154" height="176" alt="image" src="https://github.com/user-attachments/assets/aad83ff4-b6ec-4b47-852a-ef14d0feac0b" />

* Run `docker run -t --name gt_latest --detach --volume ${pwd}:/opt/data ghcr.io/gadgetron/gadgetron/gadgetron_ubuntu_rt_nocuda:latest`. If Docker is not recognized, set Docker to connect to `C:\Program Files\Docker\Docker\resources\bin` in the Environment Path in Windows. This will download and then launch the [latest Gadgetron version](https://gadgetron.readthedocs.io/en/latest/building.html) in a Docker container. It will also mount your current folder as a data folder inside the container. The installed Gadgetron software within the Docker container is as follows:

<img width="1901" height="385" alt="docker" src="https://github.com/user-attachments/assets/91d67445-be01-4eef-8d25-dd0cdedd9932" />

* Open your Docker software with the administrative account. Start the Gadgetron software within the Docker container by clicking the triangle button in the red rectangle in the figure above.      

* Run `docker exec -ti gt_latest /bin/bash` in the PowerShell (as below). This executes your Gadgetron container in PowerShell.
  
 <img width="658" height="60" alt="image" src="https://github.com/user-attachments/assets/86b6d7e1-da37-4cda-a430-474d70b83b33" />

### 2.2 Data preparation
* Place your Pulseq-generated MPRAGE (`pulseq_mprage_data.h5`) and EPI (`pulseq_epi_data.h5`) ISMRMRD data in the mounted folder, for example, in the folder: `C:\Users\chenq\Downloads\Docker_folder`.
* Run the command in Terminal: `cd /opt/data` to enter the mounted folder.
* <img width="1830" height="388" alt="image" src="https://github.com/user-attachments/assets/c5644b5f-67b6-4c8a-8c40-1ecdba6dd5a3" />
### 2.3 Gadgetron reconstruction
* For MPRAGE reconstruction, run `gadgetron_ismrmrd_client -f pulseq_mprage_data.h5 -c Generic_Cartesian_Grappa.xml -o pulseq_mprage_out.h5` in the PowerShell to call Gadgetron to reconstruct the converted Pulseq-generated MPRAGE ISMRMRD data (`pulseq_mprage_data.h5`) into an image (`pulseq_mprage_out.h5`), as shown in the figure below: 
<img width="1640" height="271" alt="image" src="https://github.com/user-attachments/assets/32706c0e-9f27-4aeb-b7ce-644086755dd3" />

* For EPI reconstruction, run `gadgetron_ismrmrd_client -f pulseq_epi_data.h5 -c epi.xml -o pulseq_epi_out.h5` in the PowerShell to call Gadgetron to reconstruct the converted Pulseq-generated EPI ISMRMRD data (`pulseq_epi_out.h5`) into an image (`pulseq_epi_out.h5`), as shown in the figure below:     
<img width="1333" height="276" alt="image" src="https://github.com/user-attachments/assets/b6286d76-a03f-47b8-ad45-a40706fca4dc" />

### 2.4 Convert Gadgetron-reconstructed images from `.h5` to NIFTI `.nii`
* Run `convert_h5_to_nifti.m` script to convert the Gadgetron-reconstructed MPRAGE (`pulseq_mprage_out.h5`) and EPI (`pulseq_mprage_gt.nii`) images from `.h5` to NIFTI `.nii` images.

## 3. Instructions for vendor-provided online reconstruction on Siemens
Vendor-provided online reconstructions, ICE and OpenRecon, are enabled for Pulseq-based MPRAGE and EPI sequences on Siemens platforms.
### 3.1 ICE reconstruction
Before executing the Pulseq-based sequences on a Siemens platform, you can enable ICE online Reconstruction following the instructions below:
* Navigate to the Special Card, set Data handling to `ICE STD` for NUMARIS/X (e.g. XA60A and XA61A), or `ICE 3D`/`ICE 2D` for MPRAGE/EPI for NUMARIS/4 (e.g. VB, VD, and VE), as below:

<img width="472.5" height="203.5" alt="image" src="https://github.com/user-attachments/assets/a6b0a2b8-63e7-480c-acdc-3b63c5972799" />

* Select `Sum-of-Square` for coil combination.
* Be sure that the maximal pixel intensity does not violate the intensity threshold of **4096**.
### 3.2 Open Recon reconstruction
* Download the Open Recon SDK from [this link](https://www.magnetom.net/t/open-recon-sdk-downloads/7223).
* Follow the `python-modules\OpenReconGettingStartedPython.pdf` document inside the SDK package to build your Open Recon container.
* You may want to configure your own Open Recon container to specify how Open Recon interfaces with the product reconstruction and the end user interface presented on the scanner. In this case, please refer to the document `OpenReconJsonConfig.pdf` located within the SDK package for instructions on configuring the Open Recon container using a JSON-formatted descriptor.
* Follow the `README.pdf` inside the SDK package to install Open Recon on your scanner. Note that Open Recon is now only partially/completely available on XA50 or higher versions. For more information, please visit the [Open Recon online forum](https://www.magnetom.net/c/openrecon/35).
* The open-source, flexible, and reproducible workflow for data acquisition, reconstruction, and post-processing based on Pulseq and Open Recon is shown below:

<img width="329.375" height="610" alt="workflow" src="https://github.com/user-attachments/assets/4ce158d2-8d91-470b-b1fc-3bcf91f45459" />

Workflow overview. **(A)** TSE sequence designed in Pulseq. **(B)** The Pulseq interpreter loads the .seq file and streams events to scanners. **(C)** Data acquisition occurs on various scanners. **(D)** Acquired data are streamed into the ICE pipeline and emitted to the OR container for custom reconstruction or post-processing using the ISMRMRD6 format. The processed data are injected into the pipeline and **(E)** converted to DICOM images and then sent to the host computer for online display. If OR is not available, raw data can be exported for offline reconstruction in the same OR container.

For more information, please refer to our ISMRM abstract: *Chen et al., Open-source, Flexible, and Reproducible Workflow for Data Acquisition, Reconstruction, and Post-processing Based on Pulseq and Open Recon, ISMRM, Hawaii, 2025*.

## 4. Example reconstructed images
This tutorial provides example Pulseq-based images reconstructed by Gadgetron, ICE, and Open Recon, acquired from a phantomn on a Cima.X 3T Siemens scanner at Freiburg, Germany.
* The MPRAGE phantom images from Siemens Cima.X 3T:

<img width="677" height="254" alt="mprage" src="https://github.com/user-attachments/assets/36839e0a-6109-4535-9805-3a40e1046df0" />

* The EPI phantom images from Siemens Cima.X 3T:

<img width="666" height="202" alt="epi" src="https://github.com/user-attachments/assets/9e1cdfda-f27a-44f0-817b-24ec08f405e2" />



 













