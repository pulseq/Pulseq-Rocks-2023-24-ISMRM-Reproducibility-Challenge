# Example Data, Images, and Scripts for Data Acquisition and Image Reconstruction Tutorials
The example data, images, and scripts for the [data acquisition](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/tree/main/data_acquisition_tutorial) and [image reconstruction](https://github.com/pulseq/Pulseq-Rocks-2023-24-ISMRM-Reproducibility-Challenge/tree/main/image_reconstruction_tutorial) tutorials can be downloaded via this link: 
https://www.dropbox.com/scl/fi/o39yenr3dz5hf8csjs32m/example_rawdata_reconImage_script.zip?rlkey=rprm958feu2y1r870ci0928oe&dl=0.
The example data were obtained from a liquid-filled phantom on a Siemens Prisma 3T scanner (Version: XA60A).
The package contains the following data, images, and scripts.
## `epi` Sub-Folder
* `epi_out.h5`: the Gadgetron-reconstructed EPI images.
* `epi_rs.seq`: the `.seq` file for the Pulseq-based EPI sequences.
* `meas_MID00090_FID10635_pulseq_epirs_tran_slc48_iso_2_8mm.dat`: the Siemens `.dat` EPI raw data.
* `pulseq2mrd_epi.m`: the script to convert the `.dat` EPI raw data to ISMRMRD data.
* `readdata.m`: the script to load the `epi_out.h5` file.
* `writeEpiRS_label.m`: the script to generate the `epi_rs.seq` file.
## `mprage` Sub-Folder
* `mprage_out.h5`: the Gadgetron-reconstructed MPRAGE images.
* `mprage.seq`: the `.seq` file for the Pulseq-based MPRAGE sequences.
* `meas_MID00077_FID10622_pulseq_t1_mprage_sag_p2.dat`: the Siemens `.dat` MPRAGE raw data.
* `pulseq2mrd_mprage.m`: the script to convert the `.dat` MPRAGE raw data to ISMRMRD data.
* `readdata.m`: the script to load the `mprage_out.h5` file.
* `writeMPRAGE_grappa_4siemens.m`: the script to generate the `mprage.seq` file.
