%% Convert the Pulseq-generated EPI Siemens .dat data to ISMRMRD data
% Author: Qingping Chen
% Email: qingping.chen@uniklinik-freiburg.de
% Date: 2025.09.16
% Prerequisites: 
% 1. Install Pulseq Matlab software from https://github.com/pulseq/pulseq
% 2. Install the ISMRMRD software from https://github.com/ismrmrd/ismrmrd
% 3. Install the mapVBVD software to load Siemens .dat raw data from https://github.com/pehses/mapVBVD
% 4. and add them to Matlab's path

%% load the .seq file and read label and sequence definitions
clear; close all; clc ;
seq = mr.Sequence() ;                           % Create a new sequence object
seqfilename = 'epi_challenge.seq' ;
seq.read(seqfilename) ;                         % load .seq file
adc_lbl = seq.evalLabels('evolution','adc') ;   % exact ADC labels

REPlbl = adc_lbl.REP ; % label index for noise repetition
SLClbl = adc_lbl.SLC ; % label index for slice
LINlbl = adc_lbl.LIN ; % label index for phase encoding
NAVlbl = adc_lbl.NAV ; % flag for navigator
AVGlbl = adc_lbl.AVG ; % flag for current acquisition
SEGlbl = adc_lbl.SEG ; % flag for segment
REVlbl = adc_lbl.REV ; % flag for reverse
LINlbl_min = min(LINlbl) ; LINlbl_max = max(LINlbl) ; % min and max LIN, C++ numbering style
SLClbl_min = min(SLClbl) ; SLClbl_max = max(SLClbl) ;
AVGlbl_min = min(AVGlbl) ; AVGlbl_max = max(AVGlbl) ;
SEGlbl_min = min(SEGlbl) ; SEGlbl_max = max(SEGlbl) ;
REPlbl_min = min(REPlbl) ; REPlbl_max = max(REPlbl) ;
Nnav = sum(NAVlbl(1:(LINlbl_max+1))) ;                % number of navigator echoes
Nrep = REPlbl_max + 1 ;                               % number of repetitions
Nscan = size(LINlbl, 2) ;                             % total scan number

% other sequence parameters
% encoded space fov
fov = 1e3 * seq.getDefinition('FOV') ;
% x: readout; y: phase encoding; z: partition/slice

e_fov_x = fov(1) ; e_fov_y = fov(2) ; e_fov_z = fov(3) ;
readout_os = seq.getDefinition('ReadoutOversamplingFactor') ;
e_fov_x = e_fov_x * readout_os ;
sliceThickness = seq.getDefinition('SliceThickness') ;
sliceGap = seq.getDefinition('SliceGap') ;
slicePositions = seq.getDefinition('SlicePositions') ; % mm
sliceOrder = slicePositions/(sliceThickness+sliceGap) + (size(slicePositions,1)-1) / 2 ;
sliceOrder = (size(slicePositions,1)-1) - sliceOrder ;
SLClbl = repmat(sliceOrder', [LINlbl_max + 1 + Nnav, Nrep]) ;
SLClbl = SLClbl(:) ;
% encoded space matrix size
e_matrixSize_x = LINlbl_max + 1 ;
e_matrixSize_y = LINlbl_max + 1 ;
e_matrixSize_z = SLClbl + 1 ;
% reconspace fov and matrix size
r_fov_x = e_fov_x / readout_os ; r_fov_y = e_fov_y ; r_fov_z = e_fov_z ;
r_matrixSize_x = e_matrixSize_x ;
r_matrixSize_y = e_matrixSize_y ;
r_matrixSize_z = e_matrixSize_z ;

%% load Siemens .dat data
datafile = dir('meas_MID00207_FID24421_pulseq_epirs_iso_2_8mm_slc48_tran.dat') ;
twix_obj = mapVBVD(datafile.name) ;                % load Siemens .dat data
phasecorData = twix_obj{end}.phasecor.unsorted() ; % navigator data
imageData = twix_obj{end}.image.unsorted() ;
Nsample = size(imageData, 1) ;  % number of ADC
Ncoil = size(imageData, 2) ;    % number of coils

% combine noise data, ref data, and image data together (typical for Siemens data)
totalData = zeros(Nsample, Ncoil, Nscan) ;  % initiliaze combined data
phasecor_count = 1 ;                        % initialize navigator counter
image_count = 1 ;                           % initialize image scan counter
for i = 1:Nscan
    if NAVlbl(i) == 1 % navigator scan
        totalData(:,:,i) = phasecorData(:,:,phasecor_count) ;
        phasecor_count = phasecor_count + 1 ;
    else % image scan
        totalData(:,:,i) = imageData(:,:,image_count) ;
        image_count = image_count + 1 ;
    end
end

%% initialize ISMRMRD data and set header

% Output file Name
filename = 'pulseq_epi_data.h5' ;
dset = ismrmrd.Dataset(filename) ;
acqblock = ismrmrd.Acquisition(Nscan) ;

% Set the header elements that don't change
acqblock.head.version(:) = 1 ;
acqblock.head.number_of_samples(:) = Nsample ;
acqblock.head.center_sample(:) = floor(Nsample/2) ;
acqblock.head.active_channels(:) = Ncoil ;
acqblock.head.read_dir  = repmat([0 0 1]', [1 Nscan]) ;
acqblock.head.phase_dir = repmat([0 1 0]', [1 Nscan]) ;
acqblock.head.slice_dir = repmat([1 0 0]', [1 Nscan]) ;

%% Loop over the acquisitions, set the header, set the data and append
for i = 1:Nscan
    % Set the header elements that change from acquisition to the next
    % c-style counting
    acqblock.head.scan_counter(i) = i-1 ;
    % Note next entry is k-space encoded line number (not i which
    % is just the sequential acquisition number)
    acqblock.head.idx.kspace_encode_step_1(i) = LINlbl(i) ;
    acqblock.head.idx.slice(i) = SLClbl(i) ;
    acqblock.head.idx.repetition(i) = REPlbl(i) ;
    acqblock.head.idx.average(i) = AVGlbl(i) ;
    acqblock.head.idx.segment(i) = SEGlbl(i) ;

    % Set the flags
    acqblock.head.flagClearAll(i) ;
    if NAVlbl(i) == 1
        acqblock.head.flagSet('ACQ_IS_PHASECORR_DATA', i) ;
    end
    if REVlbl(i) == 1
        acqblock.head.flagSet('ACQ_IS_REVERSE', i) ;
        % flip the data 
        acqblock.data{i} = flip(squeeze(totalData(:,:,i))) ; %QC: For Siemens data, flip it back! because mapVBVD flip it when loading the data
    else
        % fill the data
        acqblock.data{i} = squeeze(totalData(:,:,i)) ;
    end
end

% Append the acquisition block
dset.appendAcquisition(acqblock) ;

%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the xml header %
%%%%%%%%%%%%%%%%%%%%%%%%
% We create a matlab struct and then serialize it to xml.
% Look at the xml schema to see what the field names should be

header = [] ;

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 123194105 ; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.systemVendor = 'ISMRMRD Labs' ;
header.acquisitionSystemInformation.systemModel = 'Virtual Scanner' ;
header.acquisitionSystemInformation.receiverChannels = Ncoil ;

% The Encoding (Required)
header.encoding.trajectory = 'cartesian' ;
header.encoding.encodedSpace.fieldOfView_mm.x = fov(1) ;
header.encoding.encodedSpace.fieldOfView_mm.y = fov(1) ;
header.encoding.encodedSpace.fieldOfView_mm.z = sliceThickness * 1e3 ;
header.encoding.encodedSpace.matrixSize.x = e_matrixSize_x ;
header.encoding.encodedSpace.matrixSize.y = e_matrixSize_y ;
header.encoding.encodedSpace.matrixSize.z = 1 ;
% Recon Space
header.encoding.reconSpace.fieldOfView_mm.x = r_fov_x ;
header.encoding.reconSpace.fieldOfView_mm.y = r_fov_y ;
header.encoding.reconSpace.fieldOfView_mm.z = sliceThickness * 1e3 ;
header.encoding.reconSpace.matrixSize.x = r_matrixSize_x ;
header.encoding.reconSpace.matrixSize.y = r_matrixSize_y ;
header.encoding.reconSpace.matrixSize.z = 1 ;
% Encoding Limits
header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0 ;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = Nsample-1 ;
header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(Nsample/2) ;
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = LINlbl_min ;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = LINlbl_max ;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor((LINlbl_max+1)/2) ;
header.encoding.encodingLimits.kspace_encoding_step_2.minimum = 0 ;
header.encoding.encodingLimits.kspace_encoding_step_2.maximum = 0 ;
header.encoding.encodingLimits.kspace_encoding_step_2.center = 0 ;
header.encoding.encodingLimits.average.minimum = 0 ;
header.encoding.encodingLimits.average.maximum = 0 ;
header.encoding.encodingLimits.average.center = 0 ;
header.encoding.encodingLimits.slice.minimum = SLClbl_min ;
header.encoding.encodingLimits.slice.maximum = SLClbl_max ;
header.encoding.encodingLimits.slice.center = 0 ;
header.encoding.encodingLimits.repetition.minimum = REPlbl_min ;
header.encoding.encodingLimits.repetition.maximum = REPlbl_max ;
header.encoding.encodingLimits.repetition.center = 0 ;

gridPara = 1e6*seq.getDefinition('TrapezoidGriddingParameters') ;
rampUpTime = gridPara(1) ;
flatTopTime = gridPara(2) ;
rampDownTime = gridPara(3) ;
acqDelayTime = gridPara(4) ;
numSamples = Nsample ;
numberOfNavigator = 3 ;
dwellTime = gridPara(5)/numSamples ;
etl = LINlbl_max+1 ;
header.encoding.trajectoryDescription.identifier = 'ConventionalEPI' ;
header.encoding.trajectoryDescription.userParameterLong   = ...
    struct('value', {etl,numberOfNavigator,rampUpTime,rampDownTime,flatTopTime,acqDelayTime,numSamples},...
    'name',{'etl','numberOfNavigators','rampUpTime','rampDownTime','flatTopTime','acqDelayTime','numSamples'}) ;
header.encoding.trajectoryDescription.userParameterDouble.value   = dwellTime ;
header.encoding.trajectoryDescription.userParameterDouble.name   = 'dwellTime' ;
header.encoding.trajectoryDescription.comment = 'Conventional EPI sequence' ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_0 = 1 ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 = 1 ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2 = 1 ;
header.encoding.parallelImaging.calibrationMode = 'other' ;

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header) ;
dset.writexml(xmlstring) ;

%% Close the dataset
dset.close() ;