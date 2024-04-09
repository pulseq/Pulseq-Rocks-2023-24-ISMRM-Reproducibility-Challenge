ccc ;
%% load Siemens .dat data
datafile = dir('*meas*.dat') ; % C:\Users\chenq\Downloads\Experiments\Exp_20231102_cimax_invivo_phantom_QA\invivo\MPRAGE\pulseq_ice
twix_obj = mapVBVD(datafile.name) ; % load Siemens .dat data
phasecorData = twix_obj{end}.phasecor.unsorted() ; % navigator data
imageData = twix_obj{end}.image.unsorted() ;
Nsample = size(imageData, 1) ; % number of ADC
Ncoil = size(imageData, 2) ; % number of coils
hdr = twix_obj{end}.hdr ;
Npe = twix_obj{end}.image.NLin ;
Nslice = twix_obj{end}.image.NSli ;
Nnav = 3 ;
REP_image = reshape(twix_obj{end}.image.Rep, [Npe, Nslice]) ;
REP_phasecor = reshape(twix_obj{end}.phasecor.Rep, [Nnav, Nslice]) ;
SLC_image = reshape(twix_obj{end}.image.Sli, [Npe, Nslice]) ;
SLC_phasecor = reshape(twix_obj{end}.phasecor.Sli, [Nnav, Nslice]) ;
LIN_image = reshape(twix_obj{end}.image.Lin, [Npe, Nslice]) ;
LIN_phasecor = reshape(twix_obj{end}.phasecor.Lin, [Nnav, Nslice]) ;
NAV_image = reshape(zeros(1,size(imageData,3)), [Npe, Nslice]) ;
NAV_phasecor = reshape(ones(1,size(phasecorData,Nnav)), [Nnav, Nslice]) ;
AVG_image = reshape(twix_obj{end}.image.Ave, [Npe, Nslice]) ;
AVG_phasecor = reshape(twix_obj{end}.phasecor.Ave, [Nnav, Nslice]) ;
SEG_image = reshape(twix_obj{end}.image.Seg, [Npe, Nslice]) ;
SEG_phasecor = reshape(twix_obj{end}.phasecor.Seg, [Nnav, Nslice]) ;
REV_image = reshape(twix_obj{end}.image.IsReflected, [Npe, Nslice]) ;
REV_phasecor = reshape(twix_obj{end}.phasecor.IsReflected, [Nnav, Nslice]) ;

REPlbl = [REP_phasecor; REP_image] ; REPlbl = REPlbl(:)-1 ;
SLClbl = [SLC_phasecor; SLC_image] ; SLClbl = SLClbl(:)-1 ;
LINlbl = [LIN_phasecor; LIN_image] ; LINlbl = LINlbl(:)-1 ;
AVGlbl = [AVG_phasecor; AVG_image] ; AVGlbl = AVGlbl(:)-1 ;
SEGlbl = [SEG_phasecor; SEG_image] ; SEGlbl = SEGlbl(:)-1 ;
REVlbl = [REV_phasecor; REV_image] ; REVlbl = REVlbl(:) ;
NAVlbl = [NAV_phasecor; NAV_image] ; NAVlbl = NAVlbl(:) ;
LINlbl_min = min(LINlbl) ; LINlbl_max = max(LINlbl) ; % min and max LIN
SLClbl_min = min(SLClbl) ; SLClbl_max = max(SLClbl) ;
AVGlbl_min = min(AVGlbl) ; AVGlbl_max = max(AVGlbl) ;
SEGlbl_min = min(SEGlbl) ; SEGlbl_max = max(SEGlbl) ;
REPlbl_min = min(REPlbl) ; REPlbl_max = max(REPlbl) ;
Nscan = size(LINlbl, 1) ; % total scan number
% other sequence parameters
% encoded space fov
fov = 1e3 * [220, 220, 144] *1e-3 ; %seq.getDefinition('FOV') ;
% x: readout; y: phase encoding; z: partition/slice

e_fov_x = fov(1) ; e_fov_y = fov(2) ; e_fov_z = 3;%fov(3) ;
readout_os = 2;%seq.getDefinition('ReadoutOversamplingFactor') ;
% e_fov_x = e_fov_x * readout_os ;
% phaseResolution = 1;%seq.getDefinition('PhaseResolution') ;
% e_res_x = e_fov_x / Nsample ; % resolution in readout direction
% e_res_y = e_fov_y / (LINlbl_max+1) ; % resolution in phase encoding direction
sliceThickness = 3e-3 ;%seq.getDefinition('SliceThickness') ;
% sliceGap = 0; %seq.getDefinition('SliceGap') ;
% slicePositions = seq.getDefinition('SlicePositions') ; % mm
% sliceOrder = slicePositions/(sliceThickness+sliceGap) + (size(slicePositions,1)-1) / 2 ;
% sliceOrder = (size(slicePositions,1)-1) - sliceOrder ;
% SLClbl_first = sliceOrder(1) ;
% SLClbl_last = sliceOrder(end) ;
% SLClbl = repmat(sliceOrder', [Nscan/(SLClbl_max+1),1]) ;

% SLClbl = SLClbl(:) ;
% encoded space matrix size
e_matrixSize_x = Nsample/readout_os ;
e_matrixSize_y = Npe ;
e_matrixSize_z = 1;%SLClbl_max + 1 ;
% reconspace fov and matrix size
r_fov_x = e_fov_x / readout_os ; r_fov_y = e_fov_y ; r_fov_z = e_fov_z ;
r_matrixSize_x = e_matrixSize_x / readout_os ;
r_matrixSize_y = e_matrixSize_y ;
r_matrixSize_z = e_matrixSize_z ;

%% combine noise data, ref data, and image data together (typical for Siemens data)
totalData = zeros(Nsample, Ncoil, Nscan) ; % initiliaze combined data
phasecor_count = 1 ; % initialize navigator counter
image_count = 1 ; % initialize image scan counter
for i = 1:Nscan
    if NAVlbl(i) == 1 % navigator scan
        totalData(:,:,i) = phasecorData(:,:,phasecor_count) ;
        phasecor_count = phasecor_count + 1 ;
    else % image scan
        totalData(:,:,i) = imageData(:,:,image_count) ;
        image_count = image_count + 1 ;
    end
end

%%
% Output file Name
filename = 'epi_data.h5';
dset = ismrmrd.Dataset(filename);
acqblock = ismrmrd.Acquisition(Nscan) ;
% Set the header elements that don't change
acqblock.head.version(:) = 1 ;
acqblock.head.number_of_samples(:) = Nsample ;
acqblock.head.center_sample(:) = floor(Nsample/2) ;
acqblock.head.active_channels(:) = Ncoil ;
acqblock.head.read_dir  = repmat([0 0 1]', [1 Nscan]) ;
acqblock.head.phase_dir = repmat([0 1 0]', [1 Nscan]) ;
acqblock.head.slice_dir = repmat([1 0 0]', [1 Nscan]) ;
%%
% Loop over the acquisitions, set the header, set the data and append
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
%     if LINlbl(i) == LINlbl_min
%         acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', i) ;
%         acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', i) ;
%     end
%     if LINlbl(i) == LINlbl_max
%         acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', i) ;
%         acqblock.head.flagSet('ACQ_LAST_IN_SLICE', i) ;
%     end
%     if LINlbl(i) == LINlbl_min && SLClbl(i) == SLClbl_first
%         acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', i) ;
%     end
%     if LINlbl(i) == LINlbl_max && SLClbl(i) == SLClbl_last
%         acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', i) ;
%     end
    if NAVlbl(i) == 1
        acqblock.head.flagSet('ACQ_IS_PHASECORR_DATA', i) ;
%         acqblock.head.flagSet('ACQ_IS_RTFEEDBACK_DATA', i) ;
    end
    if REVlbl(i) == 1
        acqblock.head.flagSet('ACQ_IS_REVERSE', i) ;
        % flip and fill the data
        acqblock.data{i} = flip(squeeze(totalData(:,:,i))) ;
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

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = 123194105 ; % 3T

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.systemVendor = 'ISMRMRD Labs' ;
header.acquisitionSystemInformation.systemModel = 'Virtual Scanner' ;
header.acquisitionSystemInformation.receiverChannels = Ncoil ;

% The Encoding (Required)
header.encoding.trajectory = 'cartesian' ;
header.encoding.encodedSpace.fieldOfView_mm.x = 220 ;
header.encoding.encodedSpace.fieldOfView_mm.y = 220 ;
header.encoding.encodedSpace.fieldOfView_mm.z = sliceThickness * 1e3 ;
header.encoding.encodedSpace.matrixSize.x = 80 ;
header.encoding.encodedSpace.matrixSize.y = 80 ;
header.encoding.encodedSpace.matrixSize.z = 1 ;
% Recon Space
header.encoding.reconSpace.fieldOfView_mm.x = 220 ;
header.encoding.reconSpace.fieldOfView_mm.y = 220 ;
header.encoding.reconSpace.fieldOfView_mm.z = sliceThickness * 1e3 ;
header.encoding.reconSpace.matrixSize.x = 80 ;
header.encoding.reconSpace.matrixSize.y = 80 ;
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

rampUpTime = hdr.Meas.alRegridRampupTime(1) ;
flatTopTime = hdr.Meas.alRegridFlattopTime(1) ;
rampDownTime = hdr.Meas.alRegridRampdownTime(1) ;
acqDelayTime = hdr.Meas.alRegridDelaySamplesTime(1) ;
numSamples = hdr.Meas.alRegridDestSamples(1) ;
numberOfNavigator = 3 ;
dwellTime = 1e-3 * hdr.Meas.alDwellTime(1);
etl = numSamples ;
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


% Commented code below appears not necessary - saw this parameter after converting
% a scanner file using siemens_to_ismrmrd
% header.userParameters.userParameterLong.name = 'EmbeddedRefLinesE1' ;
% header.userParameters.userParameterLong.value = ACShw *2  ;

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header) ;
dset.writexml(xmlstring) ;

%% Write the dataset
dset.close() ;

%%
% h5_file_name = 'epi_data.h5';
% dataset = ismrmrd.Dataset(h5_file_name,'dataset');
% hdr1 = ismrmrd.xml.deserialize(dataset.readxml);
% 
% data_struct = dataset.readAcquisition();    % the acquired lines are stored in data
% mdh = data_struct.head;                  % mdh header
% total_scan_num = data_struct.getNumber;

%%
filename = 'epi_out.h5' ;
info = hdf5info(filename) ;
address_data_1 = info.GroupHierarchy.Groups(1).Groups.Datasets(2).Name ;
im = squeeze(double( hdf5read(filename, address_data_1) ) ) ;
figure ;
montage(mat2gray(im)) ;
colormap gray ;
