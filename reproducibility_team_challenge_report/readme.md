# 2023-24 ISMRM Reproducibility Challenge Report
This directory contains the report of the *Pulseq Rocks* team for the 2023-24 ISMRMRD Reproducibility Team Challenge event.

## 1. Measurements
### 1.1 Scanners
* University Medical Center Freiburg, Germany: Siemens 3T Cima.X (XA61A), Prisma (XA60A), and Trio (VB19A) scanners equipped with 20-channel, 20-channel, and 12-channel receiver-only coils, respectively.
* Massachusetts General Hospital, USA: Siemens Prisma (XA30A) scanner equipped with a 20-channel receiver-only coil.
* University of Michigan, USA: GE SIGNA UHP 3T scanner equipped with a 32-channel Nova coil.
### 1.2 Objects to be scanned
* A phantom
* The brain of a healthy human volunteer
### 1.3 Data acquisition
* Vendor-based and Pulseq-based 3D MPRAGE sequences with GRAPPA acceleration and noise scan
* Vendor-based and Pulseq-based 2D multi-slice EPI sequences with ramp sampling and a three-echo navigator

**Note**: please visit the "data_acquisition_tutorial" directory in this repository for more information.        
### 1.4 Image reconstruction
* Offline Gadgetron reconstruction for GE Pulseq-based data and Siemens vendor-based and Pulseq-based data
* Online “Image Calculation Environment” (ICE) reconstruction on Siemens scanners for Siemens-based and Pulseq-based sequences
* Online reconstruction on the GE scanner at the University of Michigan for GE-based product sequences

**Note**: please visit the "image_reconstruction_tutorial" directory in this repository for more information.    
## 2. Results and Discussions
Please read the `Pulseq_Rocks_2023_24_ISMRM_reproducibility_challenge_report.pdf` document for more information.
## 3. Conclusion
* We successfully establish an efficient, open-source, easy-to-learn MRI data acquisition and image reconstruction workflow based on Pulseq and Gadgetron and validate it for MPRAGE and EPI protocols.
* The reproducibility challenge study demonstrates that it is feasible to use this workflow to harmonize data acquisition and image reconstruction across scanner system versions, centers, and vendors.
* The preliminary results indicate that this workflow has excellent potential to enhance efficiency, transparency, and reproducibility for data acquisition and reconstruction in large-scale MRI studies.

**Note**: This work is not published yet. We plan to publish it as soon as possible.
