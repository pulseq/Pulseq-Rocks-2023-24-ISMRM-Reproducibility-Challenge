# Tutorials for converting Pulseq-generated data into ISMRMRD format
This directory contains materials and tutorials for data format standardization with Pulseq-based MPRAGE and EPI sequences on Siemens, General Electric (GE), Philips, and the United Imaging (UIH) platforms.
## 1. Required software
Please add the following three toolboxes and their subfolders to MATLAB's path:
### 1.1 [Open-source Pulseq MATLAB toolbox](https://github.com/pulseq/pulseq) to load Pulseq `.seq` files and get LABELs and sequence definitions.
### 1.2 [mapVBVD](https://github.com/pehses/mapVBVD#) software to load the Siemens `.dat` raw data. 
### 1.3 [ISMRMRD](https://github.com/ismrmrd/ismrmrd#) software for ISMRMRD conversion.

## 2. MPRAGE
We developed a 3D MPRAGE sequence using Pulseq in Matlab. The MPRAGE sequence diagram is featured with (1) a two-fold GRAPPA acceleration, (2) a noise scan, (3) and a water-only RF excitation for fat suppression. The sequence diagram and the acceleration pattern are shown below:
<img width="2563" height="902" alt="mprage_diagram" src="https://github.com/user-attachments/assets/fc26137e-333b-4d10-af44-16bae71b236f" />
* For Siemens users, please install the most updated Siemens Pulseq interpreter on your Siemens scanner. Please run the `writeMPRAGE_grappa_WE.m` script with the option `vendor='siemens'` to generate the `mprage_challenge.seq` Pulseq file. Please place the `mprage_challenge.seq` to the subfolder `%CustomerSeq%\pulseq` on your scanner. For more information about sequence installation and execution, please refer to the Pulseq C2P manual. For data acquisition, please refer to the instruction `data_acquisition_reconstruction_instruction_siemens.pdf`.
* For users from other vendors (GE, Philips, and UIH), please run the `writeMPRAGE_grappa_WE.m` script with the corresponding vendor option (e.g., `vendor='ge'`). For installation of the vendor-specific Pulseq interpreter and the execution of the Pulseq .seq files, please refer to the websites or contact the developers mentioned above in `1.3 Closed-source, vendor-dependent Pulseq interpreter`. For data acquisition, please refer to the instruction `data_acquisition_instruction_GE_Philips_UIH.pdf`.
* We also run product MPRAGE protocols for comparison. For product sequence configuration on Siemens, please refer to the document `siemens_protocol.pdf`. For product sequence configuration on other vendors, please refer to `data_acquisition_instruction_GE_Philips_UIH.pdf` for the information on the key sequence parameters.
## 3. EPI
We developed a 2D multi-slice EPI sequence with fat saturation, ramp sampling, and a three-echo navigator using Pulseq. The EPI sequence diagram is shown below.
<img width="1909" height="791" alt="epi_diagram" src="https://github.com/user-attachments/assets/3227e030-35c6-4d55-8706-a30d76c0eab1" />
* Please run the `writeEpiRS_label.m` script with the correct vendor option (e.g., `vendor = 'siemens';`) to generate the `epi_challenge.seq` Pulseq file for sequence execution. For data acquisition on Siemens platforms, please refer to the instruction `data_acquisition_reconstruction_instruction_siemens.pdf`. For data acquisition on non-Siemens platforms, please refer to the instruction `data_acquisition_instruction_GE_Philips_UIH.pdf`.
* We also run product EPI protocols for comparing the performance between the vendor-based and Pulseq-based EPI sequences. For product sequence configuration on Siemens, please refer to the document `siemens_protocol.pdf`. For product sequence configuration on other vendors, please refer to `data_acquisition_instruction_GE_Philips_UIH.pdf` for the information on the key sequence parameters.


