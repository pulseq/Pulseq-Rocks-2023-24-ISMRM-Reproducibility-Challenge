# Introduction
This repository was originally built for hosting materials, tutorials, and reports of the ***Pulseq Rocks*** team for the [2023-24 ISMRM Reproducibility Challenge](https://challenge.ismrm.org/2023-24-reproducibility-challenge/) event. The report of the Reproducibility Challenge is now in the subfolder [reproducibility_tema_challenge_report](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/blob/main/reproducibility_team_challenge_report/PulseqRocks_Chen_v2.pdf).
![reproducibility_qr_graphic_sticker](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/assets/26165904/adef1b4e-0391-4506-ab49-15db699dae20)
## Good News
Our *Pulseq Rocks* team won first place in this highly competitive Challenge because PULSEQ ROCKS!!!
## Reproducibility Challenge Team
### Team name
* *Pulseq Rocks*
### Original author sub-team members
* Qingping Chen, Frank Zijlstra, Patrick Hucker, Sebastian Littin, and Maxim Zaitsev, from University Medical Center Freiburg, Germany
### Replicator sub-team members
* Amaya Murguia, Andrea Jacobson, David Frey, Scott Peltier, and Jon-Fredrik Nielsen, from the University of Michigan (UoM), USA
* Pengcheng Xu and Berkin Bilgic, from Massachusetts General Hospital (MGH), USA
## Reproducibility Tasks
Our task is to replicate the abstract titled *"Open-Source, Cross-Platform Workflow for MRI Data Acquisition and Image Reconstruction Based on the Pulseq Framework"* (program number: 0948, ISMRM 2024) on Siemens and General Electric (GE) magnetic resonance scanners in different research centers. This abstract was awarded the **ISMRM Magna cum Laude Merit Award**.

# Update on 24.09.2025      
The repository is now updated to contain materials and tutorials of our work titled *An open-source, reproducible workflow for MRI data acquisition and reconstruction using advanced Pulseq: multi-site and cross-vendor validation*. This workflow was successfully validated on eight scanners, encompassing two magnetic fields, five sites, and four vendors.
It contains four directories:
* `data_acquisition_tutorial`: tutorial for data acquisition with MPRAGE and EPI sequences on Siemens, General Electric (GE), Philips, and the United Imaging (UIH) platforms.
* `ISMRMRD_conversion_tutorial`: tutorial for converting Pulseq-generated data to ISMRMRD format.
* `image_reconstruction_tutorial`: tutorial for image reconstruction using Gadgetron on all platforms and vendor-based online reconstruction on the Siemens platform.
* `required_software`: host all required software for this workflow.

## Motivation
Reproducibility in MRI is limited by variability in data acquisition, formatting, and reconstruction across sites and vendors. To mitigate the reproducibility challenges in MRI, we established a cross-vendor, open-source workflow for transparent and reproducible MRI acquisition and reconstruction. This workflow was evaluated with MPRAGE and EPI protocols on a phantom and a healthy human brain on eight scanners, encompassing five sites and four vendors. Vendor-native data acquisition and reconstruction protocols were also performed for comparison.
### Participating sites
* Freiburg, Germany
* Boston, USA
* Michigan, USA
* Utrecht, the Netherlands
* Shanghai, China
### Participating vendors
* Siemens
* General Electric (GE)
* Philips
* the United Imaging (UIH)

## Open-source, cross-vendor, reproducible workflow for MRI data acquisition and reconstruction using advanced Pulseq
We established a cross-vendor, open-source workflow for transparent and reproducible MRI acquisition and reconstruction. The extended Pulseq framework with advanced features (e.g., LABEL) was used to harmonize data acquisition. Pulseq-generated raw k-space data were standardized into the ISMRMRD format using LABELs and sequence definitions embedded in the Pulseq files. Gadgetron provided open-source, vendor-independent image reconstruction and post-processing. Additionally, vendor-based online reconstructions (ICE and OpenRecon) were enabled for Pulseq-based sequences on Siemens platforms. The whole workflow is shown below:

<img width="437.8333" height="676" alt="workflow" src="https://github.com/user-attachments/assets/06cc3497-7523-4c1c-abc5-289b2dd43a2e" />

Overview of the complete workflow. **(A)** EPI sequence designed using the Pulseq MATLAB toolbox. **(B)** The Pulseq interpreter loads the .seq file and streams sequence events in real time. **(C)** These events are executed on scanners from multiple sites and vendors. **(D)** Acquired raw k-space data are either converted offline to ISMRMRD format and reconstructed using Gadgetron, or streamed directly into ICE or OpenRecon for online reconstruction on Siemens platforms. **(E)** Final images are generated across sites or vendors using harmonized Gadgetron pipelines. Alternatively, ICE or OpenRecon sends DICOM images directly to the Siemens scanner console.

## Acknowledgment

This work was supported by research grants NIH R01 EB032378, NIH U24 NS120056, EU MRITwins 101078393, EURAMET 22HLT02 A4IM, and DFG INST 39 1365-1. 


## References
* Layton KJ, Kroboth S, Jia F, et al. Pulseq: a rapid and hardware-independent pulse sequence prototyping framework. Magn Reson Med. 2017;77(4):1544-1552. doi:10.1002/mrm.26235
* Nielsen JF, Noll DC. TOPPE: a framework for rapid prototyping of MR pulse sequences. Magn Reson Med. 2018;79(6):3128-3134. doi:10.1002/mrm.26990
* Roos THM, Versteeg E, Gosselink M, et al. pTx-Pulseq in hybrid sequences: accessible and advanced hybrid open-source MRI sequences on Philips scanners. Magn Reson Med. 2025;94(5):1946-1962. doi:10.1002/mrm.30601
* United Imaging. United Imaging MR ADEPT Platform | Pulseq on UIH MR. https://adept-forge.github.io/. Accessed September 17, 2025.
* Hansen MS, Sørensen TS. Gadgetron: an open source framework for medical image reconstruction. Magn Reson Med. 2013;69:1768-1776. doi:10.1002/mrm.24389
* Xue H, Inati S, Sørensen TS, Kellman P, Hansen MS. Distributed MRI reconstruction using Gadgetron-based cloud computing. Magn Reson Med. 2015;73(3):1015-1025. doi:10.1002/mrm.25213
* Inati SJ, Naegele JD, Zwart NR, et al. ISMRM raw data format: A proposed standard for MRI raw datasets. Magn Reson Med. 2017;77(1):411-421. doi:10.1002/mrm.26089
* Chen Q, Wehkamp N, Wan C, et al. Automated, open-source, vendor-independent quality assurance protocol based on the Pulseq framework. Magn Reson Mater Physics, Biol Med. 2025;38(3):533-546. doi:10.1007/s10334-025-01247-1
* Chen Q, Hucker P, Shafiekhani M, Zaitsev M. Open-source, flexible, and reproducible workflow for data acquisition, reconstruction, and post-processing based on Pulseq and Open Recon. In: Proceedings of the 33rd Annual Meeting of ISMRM, Honolulu, HI, 10-15 May, 2025. ; 2025.
* Chen Q, Zijlstra F, Hucker P, Littin S, Zaitsev M. Open-source, cross-platform workflow for MRI data acquisition and image reconstruction based on the Pulseq framework. In: Proceedings of the 32nd Annual Meeting of ISMRM, Singapore, 04-09 May, 2024. ; 2024.
