global lightChannelIndex
    Light_Channel=4;

 lightChannelIndex=4;
     percent=99;

align_opto=0;

isoline=0;
% if ~exist('end_light')
%     end_light=110;
% end

%%%%%%%%%
%% Optogenetic?
prompts = { 'Light_Channel plot (-1 if none)', 'Light channel global', 'Fraction of time to keep before and after the light artefact','align_opto (0 for start of light, 1 for end -1 for beginning)', 'isoline (0 no, 1 yes)'}; % 'time light starts', 'time light ends'};
dlgTitle='optogenetic';
numLines=1;
defaults= {num2str(Light_Channel), num2str(lightChannelIndex), num2str(percent), num2str(align_opto), num2str(isoline)};%, num2str(start_light), num2str(end_light)};
answers = inputdlg(prompts,dlgTitle,numLines,defaults);
[Light_Channel, lightChannelIndex, percent, align_opto, isoline]=deal(str2num(answers{1}),str2num(answers{2}), str2num(answers{3}), str2num(answers{4}),str2num(answers{5}));%
%[Light_Channel, start_light, end_light]=deal(str2num(answers{1}),str2num(answers{2}),str2num(answers{3}));
if Light_Channel~=-1
    [start_light, end_light,dataStruct, timeAxisStruct]=resize3(Light_Channel, dataStruct, samplingInts, timeAxisStruct, percent);
end