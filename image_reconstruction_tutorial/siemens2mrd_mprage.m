ccc ;
% convert only the image data, excluding noise data and coilsens data
datafile = dir('meas*.dat') ;
dataset = ismrmrd.Dataset('triodata.h5') ;
hdr = ismrmrd.xml.deserialize(dataset.readxml) ;
data_acq = dataset.readAcquisition() ; % data_acq.data: 192*16 3564=36*(96+3)
mdh = data_acq.head ;
twix_obj = mapVBVD(datafile.name) ;
noiseData = twix_obj{end}.noise.unsorted() ; % noise measurement data
refscandata = twix_obj{2}.refscan.unsorted() ;
imagedata = twix_obj{2}.image.unsorted() ;
nCoil = size(imagedata, 2) ;
% undersample region: 1:2:101, size=51; 135:2:237, size=52;
% 103(REF&IMA):1:134(REF), size = 32
ImaRefData = zeros(512, nCoil, 25920) ;
ImaRefData(:,:,1:51*192) = imagedata(:,:,1:51*192) ;
ImaRefData(:,:,51*192+1:(51+32)*192) = refscandata ;
ImaRefData(:,:,(51+32)*192+1:end) = imagedata(:,:,(51+16)*192+1:end) ;
% check:
error = sum(sum(abs(imagedata(:,:,51*192+1))-abs(refscandata(:,:,1))))
error1 = sum(sum(abs(imagedata(:,:,(51+16)*192))-abs(refscandata(:,:,31*192))))

totalData = zeros(512, nCoil, size(ImaRefData,3)+1) ;
totalData(:,:,1) = noiseData ;
totalData(:,:,2:end) = ImaRefData ;
nAcq = size(totalData, 3) ;
%%
filename = 'mprage_data.h5';

dset = ismrmrd.Dataset(filename);

% It is very slow to append one acquisition at a time, so we're going
% to append a block of acquisitions at a time.
% In this case, we'll do it one repetition at a time to show off this
% feature.  Each block has nYsamp acquisitions
acqblock = ismrmrd.Acquisition(nAcq);

% Set the header elements that don't change
acqblock.head = data_acq.head ;
acqblock.head.active_channels(:) = nCoil ;
%% add image scan
for i=1:nAcq
    acqblock.data{i} = squeeze(totalData(:, :, i )) ;
end
% Append the acquisition block
dset.appendAcquisition(acqblock) ;

%%%%%%%%%%%%%%%%%%%%%%%%
%% Fill the xml header %
%%%%%%%%%%%%%%%%%%%%%%%%
% We create a matlab struct and then serialize it to xml.
% Look at the xml schema to see what the field names should be
% 
% Commented code below appears not necessary - saw this parameter after converting
% a scanner file using siemens_to_ismrmrd
% header.userParameters.userParameterLong.name = 'EmbeddedRefLinesE1' ;
% header.userParameters.userParameterLong.value = ACShw *2  ;

%% Serialize and write to the data set
hdr.acquisitionSystemInformation.receiverChannels = nCoil ;
xmlstring = ismrmrd.xml.serialize(hdr);
dset.writexml(xmlstring);

%% Write the dataset
dset.close();
