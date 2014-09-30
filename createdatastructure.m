function [ dataStruct,timeAxisStruct,samplingInts ] = createdatastructure( files,paths )
% CREATEDATASTRUCTURE Creates structure variables for data and time axis
% from files entered as character input
% [dataStruct,timeAxisStruct,samplingInts] =
% createdatastructure(files,paths);
%
% *****  Author: AP ******
if iscell(files) % When only 1 file is loaded MATLAB loads as type 'char' instead of 'cell'.
    files = files'; % Columnizing the cells.
else
    files = {files}; % 'char' type to 'cell' type conversion.
end
files_char = char(files);
loopCounter = 0;
for ff = 1:size(files_char,1)
    loopCounter = loopCounter + 1;
    f = fullfile(paths,files{ff,:});
    dotind = strfind(files{ff,:},'.'); % Index of '.' in the file name
    extension = files_char(ff,dotind+1:dotind+ 3); % Type of file extension
    
    if strcmpi(extension,'abf')
        [data, samplingInts(ff), fileInfo]=abfload(f);
        samplingInts(ff) = samplingInts(ff)*1e-6; % Expresses sampling interval in seconds.
        timeAxis = (0:length(data)-1)*samplingInts(ff); % Creates a time axis vector starting with time zero
    elseif strcmpi(extension,'atf')
        [header, labels, comments, data] = import_atf(f);
        timeAxis = data(:,1);
        timeAxis = timeAxis*1e-3; % Expresses samplingInt in seconds
        samplingInts(ff) = timeAxis(2)-timeAxis(1);
        data(:,1)=[];
        break
    elseif strcmpi(extension,'mat')
        ds =  load(f);
        fn = fieldnames(ds);
        data = ds.(fn{end});
        firstCol = data(:,1);
        if any(diff(firstCol)<0)
            errordlg('First column in the selected .mat file must be time')
            break
        end
        timeAxis = data(:,1);
%         samplingInts(ff) = timeAxis(2)-timeAxis(1);
        samplingInts(ff) = mode(diff(timeAxis));
        data(:,1) = [];
    else
        errordlg('File Loading Error!!! Try changing the format of the file to ".atf" or ".mat"');
    end
    fName = files_char(ff,1:dotind-1);
    fName = num2str(fName);
    fName=['T' fName];
    dataStruct.(fName)= data;
    timeAxisStruct.(fName)= timeAxis;
end

maxSi = max(samplingInts);
nsf  = floor(1/maxSi);
if any((diff(samplingInts)~=0)) % Checking to see if all files sampled at the same interval
    answer = questdlg('FILES SAMPLED AT DIFFERENT INTERVALS, SUBSAMPLE TO EQUALIZE?','RESAMPLING','YES','NO','YES');
    if strcmpi(answer,'Yes')        
        fNames = fieldnames(dataStruct);
        for ff = 1:size(files_char,1)            
            fn = fNames{ff,:};
            [dataStruct.(fn),timeAxisStruct.(fn)] = reducedata(dataStruct.(fn),timeAxisStruct.(fn),nsf);
        end
    end    
end
samplingInts = maxSi*ones(size(samplingInts));

end


