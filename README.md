# Introduction
This repository contains materials, tutorials, and reports of the *Pulseq Rocks* team for the [2023-24 ISMRM Reproducibility Challenge](https://challenge.ismrm.org/2023-24-reproducibility-challenge/) event.
![reproducibility_qr_graphic_sticker](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/assets/26165904/adef1b4e-0391-4506-ab49-15db699dae20)
# Good News
### Our *Pulseq Rocks* team won first place in this highly competitive Challenge because PULSEQ ROCKS!!!
## Reproducibility Challenge Team
### Team name
* *Pulseq Rocks*
### Original author sub-team members
* Qingping Chen, Frank Zijlstra, Patrick Hucker, Sebastian Littin, and Maxim Zaitsev, from University Medical Center Freiburg, Germany
### Replicator sub-team members
* Amaya Murguia, Andrea Jacobson, David Frey, Scott Peltier, and Jon-Fredrik Nielsen, from the University of Michigan (UoM), USA
* Pengcheng Xu and Berkin Bilgic, from Massachusetts General Hospital (MGH), USA
## Reproducibility Tasks
Our task is to replicate the abstract titled *"Open-Source, Cross-Platform Workflow for MRI Data Acquisition and Image Reconstruction Based on the Pulseq Framework"* (program number: 0948, ISMRM 2024) on Siemens and General Electric (GE) magnetic resonance scanners in different research centers.   
### Motivation
There is currently a worldwide effort to discover potential research and promote clinical use of magnetic resonance imaging (MRI). However, a significant practical barrier is the effort required to harmonize sequences and reconstruction algorithms on vendor-specific development platforms to ensure that data are acquired and reconstructed in a consistent and reproducible manner across different platforms. In order to directly compare results between platforms or to pool data from multiple platforms for increasing statistical power, open-source science and software tools are desirable for efficiently harmonizing data acquisition and image reconstruction in multi-site MRI studies.
### Data Acquisition and Image Reconstruction Workflow
In the abstract, to alleviate the above-mentioned challenge, we develop an open-source, cross-platform, easy-to-learn workflow based on the Pulseq framework to enable efficient, transparent, reproducible data acquisition (Figure 1). We extend Pulseq to integrate with Siemens’ “Image Calculation Environment” (ICE) platform and Gadgetron to establish a complete data acquisition and reconstruction workflow. ICE is integrated into the Siemens magnetic resonance system, while Gadgetron can be used to reconstruct data from various vendors by employing vendor-independent ISMRMRD data format. Two example sequences, Magnetization Prepared RApid Gradient Echo (MPRAGE) and Echo-Planar Imaging (EPI), were developed based on the extended Pulseq and executed on three Siemens scanners to validate the workflow.
![workflow](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/assets/26165904/71345df4-6293-4298-8950-404a543cc111)
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
