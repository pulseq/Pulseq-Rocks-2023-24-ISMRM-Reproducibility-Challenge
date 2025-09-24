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
We established a cross-vendor, open-source workflow for transparent and reproducible MRI acquisition and reconstruction. The extended Pulseq framework with advanced features (e.g., LABEL) was used to harmonize data acquisition. Pulseq-generated raw k-space data were standardized into the ISMRMRD format using LABELs and sequence definitions embedded in the Pulseq files. Gadgetron provided open-source, vendor-independent image reconstruction and post-processing. Additionally, vendor-based online reconstructions (ICE and OpenRecon) were enabled for Pulseq-based sequences on Siemens platforms.



**Figure 1** Overview of the whole workflow. **(A)** MPRAGE sequence diagram and its GRAPPA pattern designed in the Pulseq Matlab software. **(B)** The Pulseq interpreter loads the .seq file and streams events to scanners. **(C)** Three different Siemens scanners for data acquisition. **(D)** The acquired data are streamed into ICE/Gadgetron for online reconstruction. If Gadgetron is not installed on the scanner, raw data can be exported to perform offline reconstruction. **(E)** ICE/Gadgetron sends the DICOM images to the MRI host computer within seconds/minutes after measurement.
## Acknowledgment
This work is supported by research grants NIH R01 EB032378 and NIH U24 NS120056. 
## References
* Van Horn JD, Toga AW. Multisite neuroimaging trials. Curr Opin Neurol. 2009;22(4):370-378.
* Layton KJ, Kroboth S, Jia F, et al. Pulseq: A rapid and hardware-independent pulse sequence prototyping framework. Magn Reson Med. 2017;77(4):1544-1552.
* Hansen MS, Sørensen TS. Gadgetron: An open source framework for medical image reconstruction. Magn Reson Med. 2013;69(6):1768-1776.
* Xue H, Inati S, Sørensen TS, Kellman P, Hansen MS. Distributed MRI reconstruction using Gadgetron-based cloud computing. Magn Reson Med. 2015;73(3):1015-1025.
* Inati SJ, Naegele JD, Zwart NR, et al. ISMRM Raw data format: A proposed standard for MRI raw datasets. Magn Reson Med. 2017;77(1):411-421. doi:10.1002/mrm.26089
* Mugler JP. Rapid Three-dimential T1-weighted MR Imaging with the MP-RAGE sequence. J Magn Reson Imaging. 1991;1(561-567).
* Stehling MK, Turner R, Mansfield P. Echo-planar imaging: Magnetic resonance imaging in a fraction of a second. Science (80- ). 1991;254(5028):43-50.
* Griswold MA, Jakob PM, Heidemann RM, et al. Generalized Autocalibrating Partially Parallel Acquisitions (GRAPPA). Magn Reson Med. 2002;47(6):1202-1210.
