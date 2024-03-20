%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%    ADAPTIVE AR BCI Experiment
%    ver.0
%
%    2024.03.11
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change for individual recorder host
recorderip = '127.0.0.1';

% Establish connection to BrainVision Recorder Software 32Bit RDA-Port
% (use 51234 to connect with 16Bit Port)
con = pnet('tcpconnect', recorderip, 51244);

% Check established connection and display a message
stat = pnet(con,'status');
if stat > 0
    disp('connection established');
end

% %%%%%%%%%%%%%% Only Testing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'IP address:','Port number:','Subject number:'};
dlgtitle = 'Input';
fieldsize = [1 45; 1 45; 1 45];
definput = {'127.0.0.1','1668','Sub01'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput);

adip = answer{1};
port = str2double(answer{2});
SubName = answer{3};

[file,path] = uigetfile({'*.mat'});

load([path,file]);

param.calibrate=true;
param.device = 3; % (All:0/Doorlock:1/AirConditioner:2/Lamp:3/Bluetoothspeaker:4) ');
socket_sender(adip, port, param.device+200);
param.Numtrial = 0;

% param.H                         = figure(1);
% set(param.H, 'color', 'w','position',[2000, 50, 900, 500]);
% 
% for i = 1:length(param.Ch)
%     param.SH(i)                       = axes;
%     hold(param.SH(i),'on');
% 
%     subplot(6,6,i,param.SH(i));
% 
%     param.h(i,1)                        = plot(nan,nan, 'parent',param.SH(i));
%     param.h(i,2)                        = plot(nan,nan, 'parent',param.SH(i));
%     param.h(i,3)                        = plot(nan,nan, 'parent',param.SH(i));
%     param.h(i,4)                        = plot(nan,nan, 'parent',param.SH(i));
% 
% end
%     set(gcf, 'Position', [200, 50, 1800, 1000])
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 10]);

param.trD.mode = 'testing';

param.dir =path; % directory where data saved
files = split(path,'\');
logpath =  cell2mat(join(files(1:end-2),'/'));

% -- parameters for adaptation
param.trD.threshold = 0.1;
param.trD.adaptmode = 'margin';
param.DSP.lambda = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% shared('shareVar');

%%
Delay = 0.7*param.Fs;

% --- Main reading loop ---
header_size = 24;
finish = false;
while ~finish
    try
        % check for existing data in socket buffer
        tryheader = pnet(con, 'read', header_size, 'byte', 'network', 'view', 'noblock');
        while ~isempty(tryheader)
            
            % Read header of RDA message
            hdr = ReadHeader(con);
            
            % Perform some action depending of the type of the data package
            switch hdr.type
                case 1       % Start, Setup information like EEG properties
                    disp('Start');
                    % Read and display EEG properties
                    props = ReadStartMessage(con, hdr);
                    disp(props);
                    
                    save('setting.mat','props');
                    % Reset block counter to check overflows
                    lastBlock = -1;
                    
                    % set data buffer to empty
                    data1s = [];
                    trig = [];
                    
                    %
                case 4       % 32Bit Data block
                    % Read data and markers from message
                    [datahdr, data, markers] = ReadDataMessage(con, hdr, props);
                    
                    % check tcpip buffer overflow
                    if lastBlock ~= -1 && datahdr.block > lastBlock + 1
                        disp(['******* Overflow with ' int2str(datahdr.block - lastBlock) ' blocks ******']);
                        
                        logger('Overflow','',logpath);
                    end
                    lastBlock = datahdr.block;
                    
                    % print marker info to MATLAB console
                    tmp = zeros(1,datahdr.points);
                    
                    if datahdr.markerCount > 0
                        disp(['markercount : ', num2str(datahdr.markerCount)])
                        for m = 1:datahdr.markerCount
                            tmp(markers(m).position+1) = str2num(markers(m).description(2:end));
                            fprintf('Trigger: %d, Location: %d .. ', tmp(tmp~=0), find(tmp~=0) );
                            fprintf('\n')
                        end
                    end
                    
                    %--- Process EEG data,
                    EEGData = reshape(data, props.channelCount, length(data) / props.channelCount);
                    data1s = [data1s EEGData];
                    trig = [trig tmp];
                    dims = size(data1s,2);
                    
                    %-- get clean data
                    if param.calibrate == 0 && dims > param.caltime
                        param = Calibration(data1s, param);
                        disp('Calibration is done');
                        param.calibrate = true;
                        data1s = [];
                        trig = [];
                    end
                    
                    %-- If trigger type 13 is detected, finish the block
                    if ~isempty(find(tmp == param.Sys(2)))
                        param.switch_on                 = true;
                        param.Numtrial                  = param.Numtrial + 1;
                        disp(['Trial:',num2str(param.Numtrial)])
                        logger(['Trial:',num2str(param.Numtrial)],'test',logpath);
                    end
                    
                    %-- check if one block end
                    if param.switch_on %
                        param.switch_on = false;
                        sig_vec                = data1s;
                        trigger                = trig;
                        
                        %%%% Give a 700ms delay on trigger : for trigger reconstruction %%%
                        
                        trigger_re = trigger;
                        trigger_re(1:Delay) = [];
                        trigger_re = [trigger_re zeros(1,Delay)];
                        
                        %-- check mode: training or testing
                        switch param.trD.mode
                            case 'testing'
                                [C, param] = P300_processing_adaptive(sig_vec(1:param.NumChIni,:),trigger_re,param);

                                save([param.dir,'/', param.SubName,'_Testing',num2str(param.Numtrial)],'sig_vec','trigger','trigger_re');
                               
                                socket_sender(adip,port,C); % send result to controller
                                
                                fprintf('Selected : %d \n',C);
                                param.prediction(param.Numtrial) = C;
                                logger(['Predicted: ',num2str(C)],'test',logpath);
                                shareVar.value = C;

                                %-- for camera control
                                if param.device == 1 && C == 3
                                    runCamera();
                                elseif param.device == 1 && C == 4
                                    closeCamera();
                                end
                                
                                
                                save([param.dir,'/param'],'param'); % save parameter
                                
                                data1s = [];
                                trig = [];
                                
                        end
                    end
                    
                case 3       % Stop message
                    save([param.dir,'/param'],'param'); % save parameter
                    
                    disp('Stop');
                    data = pnet(con, 'read', hdr.size - header_size);
                    finish = true;
                    
                otherwise    % ignore all unknown types, but read the package from buffer
                    data = pnet(con, 'read', hdr.size - header_size);
            end
            tryheader = pnet(con, 'read', header_size, 'byte', 'network', 'view', 'noblock');
        end
    catch
        er = lasterror;
        disp(er.message);
    end
end % Main loop

% Close all open socket connections
pnet('closeall');

% Display a message
disp('connection closed');

% end



%% ***********************************************************************
% Read the message header
function hdr = ReadHeader(con)
% con    tcpip connection object

% define a struct for the header
hdr = struct('uid',[],'size',[],'type',[]);

% read id, size and type of the message
% swapbytes is important for correct byte order of MATLAB variables
% pnet behaves somehow strange with byte order option
hdr.uid = pnet(con,'read', 16);
hdr.size = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
hdr.type = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
end

%% ***********************************************************************
% Read the start message
function props = ReadStartMessage(con, hdr)
% con    tcpip connection object
% hdr    message header
% props  returned eeg properties

% define a struct for the EEG properties
props = struct('channelCount',[],'samplingInterval',[],'resolutions',[],'channelNames',[]);

% read EEG properties
props.channelCount = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
props.samplingInterval = swapbytes(pnet(con,'read', 1, 'double', 'network'));
props.resolutions = swapbytes(pnet(con,'read', props.channelCount, 'double', 'network'));
allChannelNames = pnet(con,'read', hdr.size - 36 - props.channelCount * 8);
props.channelNames = SplitChannelNames(allChannelNames);
end

%% ***********************************************************************
% Read a data message
function [datahdr, data, markers] = ReadDataMessage(con, hdr, props)
% con       tcpip connection object
% hdr       message header
% props     eeg properties
% datahdr   data header with information on datalength and number of markers
% data      data as one dimensional arry
% markers   markers as array of marker structs

% Define data header struct and read data header
datahdr = struct('block',[],'points',[],'markerCount',[]);

datahdr.block = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
datahdr.points = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
datahdr.markerCount = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));

% Read data in float format
data = swapbytes(pnet(con,'read', props.channelCount * datahdr.points, 'single', 'network'));

% Define markers struct and read markers
markers = struct('size',[],'position',[],'points',[],'channel',[],'type',[],'description',[]);
for m = 1:datahdr.markerCount
    marker = struct('size',[],'position',[],'points',[],'channel',[],'type',[],'description',[]);
    
    % Read integer information of markers
    marker.size = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
    marker.position = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
    marker.points = swapbytes(pnet(con,'read', 1, 'uint32', 'network'));
    marker.channel = swapbytes(pnet(con,'read', 1, 'int32', 'network'));
    
    % type and description of markers are zero-terminated char arrays
    % of unknown length
    c = pnet(con,'read', 1);
    while c ~= 0
        marker.type = [marker.type c];
        c = pnet(con,'read', 1);
    end
    
    c = pnet(con,'read', 1);
    while c ~= 0
        marker.description = [marker.description c];
        c = pnet(con,'read', 1);
    end
    
    % Add marker to array
    markers(m) = marker;
end
end

%% ***********************************************************************
% Helper function for channel name splitting, used by function
% ReadStartMessage for extraction of channel names
function channelNames = SplitChannelNames(allChannelNames)
% allChannelNames   all channel names together in an array of char
% channelNames      channel names splitted in a cell array of strings

% cell array to return
channelNames = {};

% helper for actual name in loop
name = [];

% loop over all chars in array
for i = 1:length(allChannelNames)
    if allChannelNames(i) ~= 0
        % if not a terminating zero, add char to actual name
        name = [name allChannelNames(i)];
    else
        % add name to cell array and clear helper for reading next name
        channelNames = [channelNames {name}];
        name = [];
    end
end
end