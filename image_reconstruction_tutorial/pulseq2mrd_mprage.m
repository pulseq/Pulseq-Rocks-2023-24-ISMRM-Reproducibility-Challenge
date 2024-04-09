ccc ;
%% load .seq file and read label
seq = mr.Sequence() ;              % Create a new sequence object
seqfilename = 'mprage.seq' ;
seq.read(seqfilename,'detectRFuse') ; % load .seq file
adc_lbl=seq.evalLabels('evolution','adc') ; % exact ADC labels
figure; plot(adc_lbl.REF) ; % plot labels
hold on; plot(adc_lbl.LIN);plot(adc_lbl.NOISE) ;
plot(adc_lbl.IMA) ;plot(adc_lbl.PAR) ;
legend('REF','LIN', 'NOISE','IMA','PAR') ;
title('evolution of labels/counters') ;

NOISElbl = adc_lbl.NOISE ; % label for noise measurement
LINlbl = adc_lbl.LIN ; % label for Phase encoding
IMAlbl = adc_lbl.IMA ; % label for ACS lines for both calibration and undersampling
REFlbl = adc_lbl.REF ; % label for ACS lines for calibration
PARlbl = adc_lbl.PAR ; % label for partition encoding
LINlbl_min = min(LINlbl) ; LINlbl_max = max(LINlbl) ; % min and max LIN
PARlbl_min = min(PARlbl) ; PARlbl_max = max(PARlbl) ; % min and max PAR
Nscan = size(LINlbl, 2) ; % total scan number

%% load Siemens .dat data
datafile = dir('meas*.dat') ; % C:\Users\chenq\Downloads\Experiments\Exp_20231102_cimax_invivo_phantom_QA\invivo\MPRAGE\pulseq_ice
twix_obj = mapVBVD(datafile.name) ; % load Siemens .dat data
refData = twix_obj{end}.refscan.unsorted() ; % ACS calibration raw data
noiseData = twix_obj{end}.noise.unsorted() ; % noise measurement data
imageData = twix_obj{end}.image.unsorted() ; % regular undersampling data for imaging (a few also belongs to ACS data) 
Nsample = size(imageData, 1) ; % number of ADC
Ncoil = size(imageData, 2) ; % number of coils

% other sequence parameters
% encoded space fov
e_fov = 1e3 * seq.getDefinition('FOV') ;
% x: readout; y: phase encoding; z: partition/slice
switch seq.getDefinition('OrientationMapping')
    case 'COR'
        e_fov_x = e_fov(3) ; e_fov_y = e_fov(1) ; e_fov_z = e_fov(2) ;
    case 'SAG'
        e_fov_x = e_fov(3) ; e_fov_y = e_fov(2) ; e_fov_z = e_fov(1) ;
    otherwise
        e_fov_x = e_fov(1) ; e_fov_y = e_fov(2) ; e_fov_z = e_fov(3) ;
end
readout_os = seq.getDefinition('ReadoutOversamplingFactor') ;
e_fov_x = e_fov_x * readout_os ;
phaseResolution = seq.getDefinition('PhaseResolution') ;
e_res_x = e_fov_x / Nsample ; % resolution in readout direction
e_res_y = phaseResolution * e_res_x ; % resolution in phase encoding direction
e_res_z  = e_fov_z / (PARlbl_max+1) ; % resolution in partition encoding direciton
% encoded space matrix size
e_matrixSize_x = Nsample ;
e_matrixSize_y = e_fov_y / e_res_y ;
e_matrixSize_z = PARlbl_max + 1 ;
% reconspace fov and matrix size
r_fov_x = e_fov_x / readout_os ; r_fov_y = e_fov_y ; r_fov_z = e_fov_z ;
r_matrixSize_x = e_matrixSize_x / readout_os ;
r_matrixSize_y = e_matrixSize_y ;
r_matrixSize_z = e_matrixSize_z ;

kSpaceCenterLine = seq.getDefinition('kSpaceCenterLine') ;

%% combine noise data, ref data, and image data together (typical for Siemens data)
totalData = zeros(Nsample, Ncoil, Nscan) ; % initiliaze combined data
noise_count = 1 ; % initialize noise counter
image_count = 1 ; % initialize regular undersmapling counter
ref_count = 1 ; % initialize ACS counter (REF and IMA)
for i = 1:Nscan
    if NOISElbl(i) == 1 % noise measurement
        totalData(:,:,i) = noiseData(:,:,noise_count) ;
        noise_count = noise_count + 1 ;
    else % not noise measurement
        if REFlbl(i) == 1 % is ACS region
            if IMAlbl(i) == 1
                % both calibration and part of the undersampled pattern
                totalData(:,:,i) = refData(:,:,ref_count) ;
                ref_count = ref_count + 1 ;
                image_count = image_count + 1 ; % this is also part of the regular undersampling in imageData
            else
                % in ACS block but not part of the regular undersampling
                totalData(:,:,i) = refData(:,:,ref_count) ;
                ref_count = ref_count + 1 ;
            end
        else
            totalData(:,:,i) = imageData(:,:,image_count) ;
            image_count = image_count + 1 ;
        end
    end
end
%%
figure ;
plot(abs(totalData(:,:,1000))) ;
%%
% Output file Name
filename = 'mprage_data.h5';
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
    acqblock.head.idx.kspace_encode_step_2(i) = PARlbl(i) ;
%     acqblock.head.idx.repetition(acqno) = rep - 1 ;

    % Set the flags
    acqblock.head.flagClearAll(i) ;
    if LINlbl(i) == LINlbl_min
        acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP1', i) ;
    elseif LINlbl(i) == LINlbl_max
        acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP1', i) ;
    end
    if PARlbl(i) == PARlbl_min
        acqblock.head.flagSet('ACQ_FIRST_IN_ENCODE_STEP2', i) ;
    elseif PARlbl(i) == PARlbl_max
        acqblock.head.flagSet('ACQ_LAST_IN_ENCODE_STEP2', i) ;
    end
    %         if pe_idx(acqno) == 1 && par_idx(acqno) == 1
    %             acqblock.head.flagSet('ACQ_FIRST_IN_SLICE', acqno);
    %             acqblock.head.flagSet('ACQ_FIRST_IN_REPETITION', acqno);
    %         elseif pe_idx(acqno) == J(end) && par_idx(acqno) == Npar
    %             acqblock.head.flagSet('ACQ_LAST_IN_SLICE', acqno);
    %             acqblock.head.flagSet('ACQ_LAST_IN_REPETITION', acqno);
    %         end
    if NOISElbl(i) == 1
        acqblock.head.flagSet('ACQ_IS_NOISE_MEASUREMENT', i) ;
    end
    if REFlbl(i) == 1
        if IMAlbl(i) == 1
            % both calibration and part of the undersampled pattern
            acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION_AND_IMAGING', i) ;
        else
            % in ACS block but not part of the regular undersampling pattern Ysamp_u
            acqblock.head.flagSet('ACQ_IS_PARALLEL_CALIBRATION', i) ;
        end
    end
    % fill the data
    acqblock.data{i} = squeeze(totalData(:,:,i)) ;
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
header.encoding.encodedSpace.fieldOfView_mm.x = e_fov_x ;
header.encoding.encodedSpace.fieldOfView_mm.y = e_fov_y ;
header.encoding.encodedSpace.fieldOfView_mm.z = e_fov_z ;
header.encoding.encodedSpace.matrixSize.x = e_matrixSize_x ;
header.encoding.encodedSpace.matrixSize.y = e_matrixSize_y ;
header.encoding.encodedSpace.matrixSize.z = e_matrixSize_z ;
% Recon Space
header.encoding.reconSpace.fieldOfView_mm.x = r_fov_x ;
header.encoding.reconSpace.fieldOfView_mm.y = r_fov_y ;
header.encoding.reconSpace.fieldOfView_mm.z = r_fov_z ;
header.encoding.reconSpace.matrixSize.x = r_matrixSize_x ;
header.encoding.reconSpace.matrixSize.y = r_matrixSize_y ;
header.encoding.reconSpace.matrixSize.z = r_matrixSize_z ;
% Encoding Limits
header.encoding.encodingLimits.kspace_encoding_step_0.minimum = 0 ;
header.encoding.encodingLimits.kspace_encoding_step_0.maximum = Nsample-1 ;
header.encoding.encodingLimits.kspace_encoding_step_0.center = floor(Nsample/2) ;
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = LINlbl_min ;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = LINlbl_max ;
header.encoding.encodingLimits.kspace_encoding_step_1.center = kSpaceCenterLine ;
header.encoding.encodingLimits.kspace_encoding_step_2.minimum = PARlbl_min ;
header.encoding.encodingLimits.kspace_encoding_step_2.maximum = PARlbl_max ;
header.encoding.encodingLimits.kspace_encoding_step_2.center = floor((PARlbl_max+1)/2) ;
header.encoding.encodingLimits.repetition.minimum = 0 ;
header.encoding.encodingLimits.repetition.maximum = 0 ;
header.encoding.encodingLimits.repetition.center = 0 ;

header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_0 = 1 ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_1 = 2 ;
header.encoding.parallelImaging.accelerationFactor.kspace_encoding_step_2 = 1 ;
header.encoding.parallelImaging.calibrationMode = 'embedded' ;

% header.measurementInformation.measurementID = '000002' ;
% header.measurementInformation.patientPosition = 'HFS' ;
% header.measurementInformation.measurementDependency = struct('dependencyType', {'SenMap','Noise'},'measurementID',{'000001','000001'}) ;

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
% h5_file_name = 'mprage_data.h5';
% dataset = ismrmrd.Dataset(h5_file_name,'dataset');
% hdr = ismrmrd.xml.deserialize(dataset.readxml);
% 
% data_struct = dataset.readAcquisition();    % the acquired lines are stored in data
% mdh = data_struct.head;                  % mdh header
% total_scan_num = data_struct.getNumber;

filename = 'mprage_out.h5' ;
info = hdf5info(filename) ;
address_data_1 = info.GroupHierarchy.Groups(1).Groups.Datasets(2).Name ;
im = squeeze(double( hdf5read(filename, address_data_1) ) ) ;
% im = permute(im, [3 2 1]) ;
figure ;
montage(mat2gray(im(:,:,:))) ;
colormap gray ;
