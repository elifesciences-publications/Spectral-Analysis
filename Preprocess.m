% Preprocess - Conditions multifile signals loaded by LoadFiles
% Looks at the multifile data structure created by LoadFiles and preprocesses
% the timeseries in preparation for subsequent wavelet analysis

% Author: Avinash Pujala


%% Declaring a Few Global Variables
global lightChannelIndex globData globTime;
lightChannelIndex = -1;

%% Extracting Basic File Information
fNames = fieldnames(dataStruct); % Converting to 'char' type variable for ease of size indexing
nFiles = length(fNames);

%% Minor Options
stopband =[58 62];
denoise = 'n'; % 'y' = removes 60Hz by stopbanding b/w values specified in "stopband"; 'n' = no fitering;
threshdenoising ='n'; %%%% 'y' = gives the option of setting threshold for artifact truncation; 'n' =  no option;

peakStringency = 1;
yShift = 25;

%% Checking Sampling Interval Consistency
sichck = diff(samplingInts);
if any(sichck ~= 0) % Given the code in createdatastructure.m, I am not sure that this will ever be true.
    errordlg('Signals Sampled at Different Intervals!')
    break
else
    samplingInt = samplingInts(1);
end

%% Comparing Signal Lengths to Determine the Smallest Permissible Common Signal Length
minMat=[]; maxMat =[]; loopCounter = 0;
clear data
data(nFiles,1) = struct;
time = cell(nFiles,1);
for fileNum = 1:nFiles
    data(fileNum).raw = dataStruct.(fNames{fileNum,:});
    time{fileNum} = timeAxisStruct.(fNames{fileNum,:});
    if loopCounter < 1 % This is to prevent asking for artifact detection mode repeatedly by applying the first detection mode to all the files
        [data(fileNum).raw,time{fileNum},~,selection] = artifactalign(data(fileNum).raw,time{fileNum});
    else
        [data(fileNum).raw,time{fileNum}] = artifactalign(data(fileNum).raw,time{fileNum},selection);
    end
    minMat =[minMat; time{fileNum}(1)];
    maxMat =[maxMat; time{fileNum}(end)];    
    loopCounter = loopCounter+1;
end

%% Setting the Default Values for Processing Parameters to Last Stored Values

if ~exist('ch')
    ch = ['1 2 3 4'];
end

if ~exist('hpf')
    hpf = 50;
end

if ~exist('lpf')
    lpf = 2;
end
if ~exist('firstCommonTime')
    firstCommonTime = max(minMat); % Starting point for common time vector
end

if ~exist('lastCommonTime')
    lastCommonTime = min(maxMat); % End point for common time vector
end

if ~exist('freqRange')
    freqRange = [0.1 4];
end

if ~exist('stringency')
    stringency = 1;
end

if ~exist('isoThresh')
    isoThresh = 0.5;
end

if ~exist('phaseType')
    phaseType = 'All';
end

if ~exist('traceType')
    traceType = 'Smooth';
end

if ~exist('icc')
    icc = -1;
end

if ~exist('hpf_icc')
    hpf_icc = 0.0001;
end

if ~exist('lpf_icc')
    lpf_icc = 2;
end

if ~exist('lightChannelIndex')
    lightChannelIndex = -1;
end

%% Inputting Processing Paramters
maxFr = floor(0.5*(1/samplingInt));
minFr = 2*(1/(min(maxMat)-max(minMat)));
trPrompt = ['Time Range (in sec): Permissible Range = [' num2str(max(minMat)) '    '  num2str(min(maxMat)) ' ]'];
frPrompt =  ['Freq Range (in Hz): Permissible Range = [' num2str(minFr) '    '  num2str(maxFr) ' ]'];
prompts = {'Channel Pairs to Analyze (Include IC channel)', 'Highpass Before Rectification (If < 20, no rectification)',...
    'Lowpass after Rectification', trPrompt, frPrompt, 'Statistical Stringency (1 = 95% CI, 2 ~ 98% CI )',...
    'Isoline Threshold (Threshold for Contour Line Plotted on the Avg XW Spectrum Such That Half the Power is Below this Line)',...
    'Phase Filtering ("Alt" = Alternation, "Synch" = Only Synchronous, "All" = No Phase Filtering)'...
    'Trace type for plotting ("Raw" or "Smooth")','Specify Intracellular Channel(s), if any (''-1'' = no channel)','IC Highpass','IC Lowpass'...
    'Specify Light Channel Index so as to Overlay Trace on Avg XWTs ("-1") = no channel'};
dlgTitle = 'Processing Parameters';
numLines = 1;
defaults = {num2str(ch),num2str(hpf),num2str(lpf),num2str([max(minMat) min(maxMat)]),...
    num2str(freqRange),num2str(stringency),num2str(isoThresh), phaseType,traceType,num2str(icc),...
    num2str(hpf_icc),num2str(lpf_icc),num2str(lightChannelIndex)};
answers = inputdlg(prompts,dlgTitle,numLines,defaults);

[ch,hpf,lpf] = deal(str2num(answers{1}),str2num(answers{2}),str2num(answers{3}));
[timeRange,freqRange,stringency,isoThresh] = deal(str2num(answers{4}),str2num(answers{5}),str2num(answers{6}),str2num(answers{7}));
[phaseType,traceType] = deal(answers{8},answers{9});
[icc,hpf_icc,lpf_icc,lightChannelIndex] = deal(str2num(answers{10}),str2num(answers{11}),str2num(answers{12}),str2num(answers{13}));

%%% Channel Number Check
if numel(ch)<2
    errordlg('Select at least two channels to compare! Enter channel number twice for autowavelet')
    break
end

if (max(minMat)-timeRange(1))> samplingInt
    error_msg = {'Time range out of bounds! First value of the variable "timeRange" must at least equal' num2str(max(minMat))};
    errordlg(error_msg,'Time Range Out of Bounds!')
    break
elseif (timeRange(2)-min(maxMat))> samplingInt;
    error_msg = {'Time range out of bounds! Second value of the variable "timeRange" cannot exceed' num2str(min(maxMat))};
    errordlg(error_msg,'Time Range Out of Bounds!')
    break
else
    firstCommonTime = timeRange(1);
    lastCommonTime = timeRange(end);
end

commonTime = firstCommonTime:samplingInt:lastCommonTime; % Creates common time vector
lenTime =length(commonTime);

%% File Loop
icc_pos = find(ch==icc);
extra_pos = find(ch~=icc);
ch_extra = ch;
ch_extra(icc_pos)=[];
% artChk = questdlg('Auto-remove Stimulus Artifacts?');
if selection ~=4
    artChk = 'yes';
else
    artChk = 'no';
end
% temp = cell(nFiles,1);
% signal = cell(nFiles,1);
for fileNum = 1:nFiles
    data(fileNum).fName = fNames{fileNum};
    %%%% Truncating Signals to Common Time Portion
    [fpt,lpt] = deal([]);
    fstr = num2str(fileNum);
    dTime = diff(commonTime);
    if any(dTime<=0)
        errordlg('Time Axis Vector Inconsistent: Please Make Sure Time Axis in Each File Has Evenly Spaced Ascending Values;')
    end    
    fpt = find(time{fileNum} >= firstCommonTime,1);
    lpt = fpt + lenTime;
    commonPts = fpt:lpt;
    lenDiff = length(commonTime) - length(commonPts);
    lpt = lpt+lenDiff;
    commonPts = fpt:lpt;
    if lightChannelIndex ~=-1
        data(fileNum).raw = data(fileNum).raw(commonPts,[ch lightChannelIndex]);% Inserting light channel into 'data'
    else
        data(fileNum).raw = data(fileNum).raw(commonPts,ch);
    end
    % Highpass Filtering Signals
    data(fileNum).hp(:,extra_pos) = chebfilt(data(fileNum).raw(:,extra_pos),samplingInt,hpf,'high');
%     temp{fileNum}(:,extra_pos) = chebfilt(data{1}(:,extra_pos),samplingInt,hpf,'high');
    if ~isempty(icc_pos)
%         temp{fileNum}(:,icc_pos) = chebfilt(data{fileNum}(:,icc_pos),samplingInt,hpf_icc,'high');
        data(fileNum).hp(:,icc_pos) = chebfilt(data(fileNum).raw(:,icc_pos),samplingInt,hpf_icc,'high');
    end
    
    if strcmpi(artChk,'yes')
        data(fileNum).hp = autoartremove(data(fileNum).hp,commonTime);
%         temp{fileNum} = autoartremove(temp{fileNum},time);
    end
    
    %%%% Filtfilt.m Error Message
    if any(isnan(data(fileNum).hp(:)))
        errordlg('Signal Filtering Error! Please re-specify the time range')
    end
    
    if strcmpi(threshdenoising,'y')
%         temp{fileNum} = chebfilt(overthreshremove(temp{fileNum},time),samplingInt,hpf,'high');
        data(fileNum).hp = chebfilt(overthreshremove(data(fileNum).hp,time),samplingInt,hpf,'high');
    end
    %%%% Filtfilt.m Error Message
    if any(isnan(data(fileNum).hp(:)))
        errordlg('Signal Filtering Error! Please re-specify the time range')
    end
    %%%% Removing 60Hz noise
    if strcmpi(denoise,'y')
        data(fileNum).hp = chebfilt(double(data(fileNum).hp),samplingInt,stopband,'stop');
%         temp{fileNum} = chebfilt(temp{fileNum},samplingInt,stopband,'stop');
        %%%% Filtfilt.m Error Message
        if any(isnan(data(fileNum).hp(:)))
            errordlg('Signal Filtering Error! Please re-specify the time range')
        end
    end
    data(fileNum).smooth = data(fileNum).hp;
    %%%% Filtfilt.m Error Message
    if any(isnan(data(fileNum).hp(:)))
        errordlg('Signal Filtering Error! Please re-specify the time range')
    end
    
    %%%% Signal Rectification or absolute values
    if hpf >= 20 % Rectification only when signal is highpassed over 20Hz
        %%%%%% Altneratively, use this for absolute values %%%%
%         signal{fileNum} = abs(signal{fileNum});
        data(fileNum).smooth = abs(data(fileNum).smooth);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
%     signal{fileNum}(:,extra_pos) = chebfilt(signal{fileNum}(:,extra_pos),samplingInt,lpf,'low');
    data(fileNum).smooth(:,extra_pos) = chebfilt(data(fileNum).smooth(:,extra_pos),samplingInt,lpf,'low');
    if ~isempty(icc_pos)
        data(fileNum).smooth(:,icc_pos) = chebfilt(data(fileNum).smooth(:,icc_pos),samplingInt,lpf_icc,'low');
    end
    %%%% Filtfilt.m Error Message
    if any(isnan(data(fileNum).smooth(:)))
        errordlg('Signal Filtering Error! Please re-specify the time range')
    end
end

data(1).time = commonTime;
globData = data(1).raw;
globTime = data(1).time;


%% Choice of Removing Light Artifacts
% ques = questdlg('Automatically remove light artifacts?','LIGHT ARTIFACT REMOVAL','Yes','No','No');
% if strcmpi(ques,'yes')
%     dataStruct_mod = slowartifactremove(dataStruct,samplingInt);
% else
%     dataStruct_mod = dataStruct;
% end

%% Running Follow Up Code
% xwplotmd