function varargout = ComputeXWT(varargin)

% My custom written XWT based on xwt.m by Grinsted et al
% [Wxy,freq,coi,sig95] = ComputeXWT(x,y,time,freqRange,dj,stringency,phaseType,sigmaXY);

%% Fixed and variable parameters
fourier_factor      = 1.0330; % Conversion factor for changing wavelet scales to periods (true for wavenumber = 6)
Pad                 = 1; % Zero padding of signals (Pad = 1: Pad with zeroes; Pad = 0: Don't zero pad)
motherWavelet       = 'Morlet'; %%%('Morlet', 'Paul','DOG') - For now use only Morlet, other wavelets will give erroneous results
figdisp             = 'n'; %%% ('n' = does not display figures; [] = displays figs );
plotfig             = 'no';
yShifter            = 1.5; %(Factor by which the max value of a trace is multiplied...
%  and shifted along the y-axis to prevent overlap; default = 1.5)
powerSpectrumType   = 'both';  %('linear', 'log', 'both', 'none')

% peakDetectionThresh = 0.5; % Determines the peak detection for global wavelet spectrum plotted to the right of XWT
% time_axis_xticks    = 'regular'; %('train' - displays the stimulus train, 'regular' - displays time regularly; default:'train')


%% Input check
if nargin < 4
    errordlg('At least 4 inputs required')
    return;
elseif nargin == 4
    dj = 1/64;
elseif nargin == 5
    dj = varargin{5};
    stringency = 1;
    phaseType = 'synch';
elseif nargin == 6
    dj = varargin{5};
    stringency = varargin{6};
    phaseType = 'synch';
elseif nargin == 7
    dj = varargin{5};
    stringency = varargin{6};
    phaseType = varargin{7};
elseif nargin == 8
     dj = varargin{5};
    stringency = varargin{6};
    phaseType = varargin{7};
    sigmaXY = varargin{8};
else
    errordlg('Too many inputs')
    return;
end

x = varargin{1}; x = x(:);
y = varargin{2}; y = y(:);
time = varargin{3}; time = time(:);
samplingInt = mode(diff(time));

freqRange = varargin{4}; % Power for frequencies only within this range will be calculated and displayed
if freqRange(2) >= floor(1/samplingInt)
    errordlg('High frequency value above Nyquist limit, please respecify')
end

scaleRange = 1./(freqRange*fourier_factor); % Scale range corresponding to frequency range.
S0 = min(scaleRange);
MaxScale = max(scaleRange);

%% PHASE DISPLAY
if ~(strcmpi(phaseType,'alt')|| strcmpi(phaseType,'synch') || strcmpi(phaseType,'all'))
    errordlg('Input for the variable "phaseType" not specified properly')
end

%% POWER SPECTRUM DISPLAY
if strcmpi(powerSpectrumType,'linear')|| strcmpi(powerSpectrumType,'log') || strcmpi(powerSpectrumType,'both') || strcmpi(powerSpectrumType,'none')
else
    errordlg('Input for the variable "powerSpectrumType" not specified properly')
end

%% XWT
% clear statMat chLabelMat fNamesMat
% nChannelPairs = size(sigMat,3)-1;
% nFiles = size(sigMat,2);
% endCol = nFiles*nChannelPairs+1;
% statMat = cell(15, nFiles*nChannelPairs + 1 + 2*nChannelPairs);
% % statMat = cell(15, size(sigMat,2)*(size(sigMat,3)-1)+ 1 + 2*(size(ch,2)-1));
% statMat(:,1)= deal({'Peak f','Mean f', 'Std f','Peak Ph', 'Mean Ph', 'Std Ph',...
%     'Mean Pow','Std Pow','Synch Pow','Alt Pow','Alt Pow/Synch Pow' ,'Tot Pow','Tot Pow Ratio','Time Range','Freq Range'});
% % chLabelMat =cell(1, size(sigMat,2)*(size(sigMat,3)-1)+1+2*(size(ch,2)-1));
% chLabelMat =cell(1, nFiles*nChannelPairs + 1 + 2*nChannelPairs);
% chLabelMat{1} = 'Channels';
% % fNamesMat =cell(1, size(sigMat,2)*(size(sigMat,3)-1)+1+2*(size(ch,2)-1));
% fNamesMat = cell(1,nFiles*nChannelPairs + 1 + 2*nChannelPairs);
% fNamesMat{1}='File Names';
% cellNum = 1;
% time_varying_power_mat =[];
% fileCounter = 0;
% clear Wxy3d % 3D array with stacking Wxy matrices generated for each pair of channels from each file along the z-dimension
% clear sigxy


x = [time(:) x(:)];
y = [time(:) y(:)];
[Wxy,period,~,coi,sig95]= xwt(x,y,Pad, dj,'S0',S0, 'ms',MaxScale, 'Mother', motherWavelet);
freq =1./period;
if nargin ==8
    Wxy = Wxy/sigmaXY;
else
    [~,~,sigmaX] = ZscoreByHist(x(:,2));
[~,~,sigmaY] = ZscoreByHist(y(:,2));
Wxy = Wxy./(sigmaX*sigmaY);
end

%% Removing regions outside COI
ftmat = repmat(freq(:), 1, length(time));
coimat = repmat(1./coi(:)',length(freq), 1);
% Wxy(ftmat < coimat) = 0; % Removing regions outside of COI

%% Removing insignificant points
% Wxy(sig95 < stringency) = 0;
Wxy(abs(Wxy) < stringency)= 0;

%% Phase-based filtering
a  = angle(Wxy);
if strcmpi(phaseType,'alt')
    Wxy(abs(a)<=(0.5*pi))=0; % Keeps only those matrix elements with angles >= (0.75*pi) = 135 deg
elseif strcmpi(phaseType,'synch')
    Wxy(abs(a)>pi/2)=0; % Keeps only those matrix elements with angles <= pi/2
elseif strcmpi(phaseType,'all')
    % Do nothing
else
    errordlg('Please input proper phase type')
end

varargout{1} = Wxy;
varargout{2} = freq;
varargout{3} = coi;
varargout{4} = sig95;


end
%% Global Power Spectrum

% globalPowSpec = sum(abs(Wxy),2);
% globalPowSpec_maxnorm = globalPowSpec/ max(globalPowSpec); % Power units indicate cumulative probability
% normlog_globalPowSpec = log2(globalPowSpec);
% blah = normlog_globalPowSpec;
% blah(blah==-inf) = NaN;
% blah = (blah-min(blah))/(max(blah)-min(blah)); % Normalization by subtracting min and dividing by amplitude
% normlog_globalPowSpec = blah;
% [maxtab,~] = peakdet(globalPowSpec_maxnorm, peakDetectionThresh);
% if isempty(maxtab)
%     [maxtab(:,2), maxtab(:,1)] = max(globalPowSpec_maxnorm(:));
%     maxtab = maxtab(1,:);
% end
% 
% pv = find(maxtab(:,2)== max(maxtab(:,2)));
% pf = round(freq(maxtab(pv,1))*100)/100;

%
%
%         %% PLOTTING FIGURES
%
%         if lower(figdisp) =='y'
%
%             figure('Name', ['XW Power, File ' num2str(fileNum) ', Channels '...
%                 num2str(ch(chNum)) ' vs ' num2str(ch(chNum+1))],'color','w','renderer','painter')
%             figPos = get(gcf,'position');
%             if strcmpi(powerSpectrumType,'none')
%                 set(gcf,'position',[figPos(1)*0.6 figPos(2)-figPos(4)/2 figPos(3)*1.3...
%                     figPos(4)*2]);
%             else
%                 set(gcf,'position',[figPos(1)*0.6 figPos(2)-figPos(4)/2 figPos(3)* 1.7...
%                     figPos(4)*1.5]);
%             end
%             %% XWT AXES
%             ax1 = axes; box off
%             aPos = get(ax1,'position');
%             if strcmpi(powerSpectrumType,'none')
%                 aPos = [aPos(1)*0.8 aPos(2)+ aPos(4)*(1/3) aPos(3) aPos(4)*0.75];
%             else
%                 aPos = [aPos(1)*0.8 aPos(2)+ aPos(4)*(1/3) aPos(3)*(0.95) aPos(4)*0.75];
%             end
%             set(ax1,'position', aPos,'drawmode','fast')
%         end
%
%         if lower(figdisp)=='y'
% %             plotwave(Wxy,time_reduced,period,coi,sig95,...
% %                 sigmas(fileNum,chNum), sigmas(fileNum,chNum+1)
% plotwave(Wxy,time,period,coi,sig95,std(x),std(y))
%             set(ax1,'color','k','xtick',[], 'xcolor','w','ycolor','k','drawmode','fast')
%             xlim([time_reduced(1) time_reduced(end)]) %%%% This line is NECESSARY to ensure that x-axis is aligned with traces below
%             xlabel('')
%             ylims1 = get(ax1,'ylim');
%             yticklabels_ax1 = str2num(get(ax1,'yticklabel'));
%             dyticklabels_ax1 = diff(yticklabels_ax1);
%             %         if numel(dyticklabels_ax1)==1,yScale = 'log';
%             %         elseif dyticklabels_ax1(1)== dyticklabels_ax1(2), yScale = 'linear';
%             %         else yScale = 'log';
%             %         end
%             yScale = 'log';
%             yl = ylabel('Frequency (Hz)','fontsize',14);
%             ylpos = get(yl,'pos');
%             aPos = get(ax1,'position');
%
%             %% XW POWER SPECTRUM AXES
%             if strcmpi(powerSpectrumType,'none')
%             else
%                 ax2 = axes; hold on, box off
%                 aPos2 = get(ax2,'pos');
%                 aPos2 = [aPos(1) + aPos(3) aPos(2) aPos2(3)*(0.15) aPos(4)];
%                 set(ax2,'pos', aPos2, 'color','none','tickdir','out','fontsize',11,'drawmode','fast');
%                 xlabel([{'Normalized'}; {'Power'}])
%                 if imag(maxtab(:,2))
%                     maxtab(:,2)= 1+i;
%                     peakFreqs = 1+i;
%                 else
%                     peakFreqs = freq(maxtab(:,1));
%                 end
%                 if strcmpi(yScale,'linear')
%                     logPeakFreqs = peakFreqs;
%                     freq2 = freq;
%                 else
%                     logPeakFreqs = log2(peakFreqs);
%                     freq2 = log2(freq);
%                 end
%                 switch powerSpectrumType
%                     case 'linear'
%                         plot(norm_globalPowSpecow,freq2,'k','linewidth',2)
%                         leg = {'Linear Spectrum'}
%                     case 'log'
%                         plot(normlog_globalPowSpecow,freq2,'k','linewidth',2)
%                         set(ax2,'xtick',[],'xcolor','w','drawmode','fast')
%                         leg = {'Log Spectrum'};
%                         ax2b = axes('position',aPos2,'xaxislocation','bottom',...
%                             'yaxislocation','right','color','none','xscale','log','ytick',[],'ycolor','k','fontsize',11);
%                         xt = [0.25 0.5 1];
%                         set(ax2b,'xtick',xt);
%                         xlabel([{'Normalized'}; {'Power'}])
%                     case 'both'
%                         plot(globalPowSpec_maxnorm, freq2,'k','linewidth',2)
%                         plot(normlog_globalPowSpec, freq2,'k:','linewidth',2,'parent',ax2)
%                         leg ={'Linear'; 'Log'};
%                         ax2b = axes('position',aPos2,'xaxislocation','top',...
%                             'yaxislocation','right','tickdir','out','color','none','xscale','log','ytick',[],'ycolor','w','fontsize',11);
%                         xt = [0.25 0.5 1];
%                         set(ax2b,'xtick',xt,'xlim',[0 1]);
%                         hold off
%                 end
%
%                 set(ax2,'ylim',ylims1,'xlim',[0 1],'xtick',[0.5 1],'ytick',[],'drawmode','fast')
%                 xvals = 0.35*ones(size(peakFreqs));
%
%                 clear txt
%                 for pf = 1:length(peakFreqs)
%                     txt{pf} = [num2str(round(peakFreqs(pf)*100)/100) ' Hz'];
%                     % txt = {round(peakFreqs*100)/100};
%                 end
%
%                 text(xvals,logPeakFreqs,txt,'fontsize',11,'color','r','parent',ax2);
%                 legend(ax2,leg,'fontsize',10) % This line needs to be here to legend can be moved by hand after fig is generated
%
%             end
%
%             %% TIME SERIES AXES
%             ax3=axes; hold on, box off
%             aPos3 = get(ax3,'position');
%             aPos3 = [aPos(1) aPos3(2) aPos(3) aPos3(4)*(1/3)];
%             set(ax3,'position',aPos3,'tickdir','out','color','w','ycolor','w','drawmode','fast');
%
%
%             if strcmpi(traceType,'raw')
%                 fstr = num2str(fileNum);
%                 tempSig = eval(['temp' num2str(fileNum) '(:,chNum);']);
%                 tempSig = truncatedata(tempSig,time,[firstTime lastTime]);
%
%                 tempTime = time; % Adding firstTime to the time vector here
%                 % will set the time of the first stimulus in the stimulus
%                 % train to a value of zero
%
%                 plot(tempTime,tempSig + (yShifter/2)*max(tempSig),'k','linewidth',1.5)
%                 tempSig = eval(['temp' num2str(fileNum) '(:,chNum+1);']);
%                 tempTime = linspace(firstTime,lastTime,length(tempSig));
%                 %                 tempSig = zscore(truncatedata(tempSig,time,[firstTime lastTime]));
%                 tempSig = truncatedata(tempSig,time,[firstTime lastTime]);
%                 plot(tempTime,tempSig-(yShifter/2)*max(tempSig),...
%                     'k','linewidth',1.5)
%
%             elseif strcmpi(traceType,'smooth')
%                 tempSig(:,chNum) = zscore(sigMat(:,fileNum,chNum));
%                 tempTime = linspace(firstTime,lastTime,length(tempSig)); % Adding firstTime to the time vector here
%                 % will set the time of the first stimulus in the stimulus
%                 % train to a value of zero
%                 %                 plot(tempTime,tempSig(:,chNum),colors(chNum),'linewidth',1.5)
%                 plot(tempTime,tempSig(:,chNum)+2.5,'k','linewidth',1.5)
%                 tempSig(:,chNum+1) = zscore(sigMat(:,fileNum,chNum+1));
%                 plot(tempTime,tempSig(:,chNum+1)-2.5,'k','linewidth',1.5)
%             end
%             yl2 = ylabel([{'Normalized'};{'Amplitude'}],'fontsize',14,'color','k');
%             ylpos2 = get(yl2,'pos');
%             ylpos2_mod = ylpos2;
%             ylpos2_mod(1) = tempTime(1)-abs(ylpos(1))+ 0.65; ylpos2_mod(2)=0;
%             set(yl2,'pos',ylpos2_mod);
%             set(ax3,'ytick',[])
%             axis([tempTime(1) tempTime(end) -inf inf])
%             set(gca,'fontsize',14)
%             switch lower(time_axis_xticks)
%                 case 'regular'
%                     %                 xtick = fix(linspace(time(1),time(end),10));
%                     %                 dxt = xtick(2)-xtick(1);
%                     %                 xtick = xtick(1):dxt:xtick(end);
%                     %                 set(ax3,'ytick',[],'xtick',xtick)
%                     set(ax3,'ytick',[])
%                 case 'train'
%                     stimtrain % Running this program generates 'tStimArts' a vector of stimulus artifact times
%                     xlabel(['Stim Train (' num2str(stimDur) 'sec) @ ' num2str(stimFreq) ' Hz'],'fontsize',14)
%                     set(ax3,'ytick',[],'xticklabel',[],'xtick',[tStimArts - tStimArts(1)])
%                     hold off;
%             end
%         end
%
%
%         %% TIME-VARYING MEAN FREQUENCIES AND XW POWERS
%
%         [mfvec,pfvec] = instantaneouswavefreq(Wxy,freq);
%
%         eval(['time_varying_meanfreq_f' num2str(fileNum) 'ch' num2str(chNum)...
%             num2str(chNum+1) ' = mfvec;']);
%         eval(['time_varying_pfreq_f' num2str(fileNum) 'ch' num2str(chNum)...
%             num2str(chNum+1) ' = pfvec;']);
%         tvpower = instantaneouswavepow(Wxy);
%
%         eval(['time_varying_power_f' num2str(fileNum) 'ch' num2str(chNum)...
%             num2str(chNum+1) ' = tvpower;']);
%         %         eval(['tvpf_comb' num2str(fileNum) 'ch = tvpower;']);
%         eval(['time_varying_power_mat = [time_varying_power_mat; time_varying_power_f' num2str(fileNum) 'ch' num2str(chNum)...
%             num2str(chNum+1) '];']);
%
%         %% OPTION FOR DYNAMIC FREQ AND PHASE PLOTS
%
%         switch plotfig
%             case 'Yes'
%                 dynamicfreqpowplot
%             case 'No'
%         end
%
%
%         Wxy3d.(chStr)(:,:,fileNum) = [Wxy];
%         sigxy.(chStr)(fileNum,:) = sigmas(fileNum,chNum)* sigmas(fileNum,chNum+1);
%
%     end
% end
%
%
% %% Averaged XW Plots
%  Wxy = W_coi_sig_alt;
% [Wxy_avg,Wxy_iso,masterVar,W_gm,S,S_prelog,Wxy_avg_channels]  = avgxwt(Wxy,freq,time_reduced,coi,sigMat,isoThresh,ch);
% avgCheck = 1;
%
%
% %% Updating Statistics (statMat) and Storing in an Excel File
% clear intraIsoVars
% intraIsoVars = cell(size(statMat,1)+2,2*(size(ch,2)-1));
% shift = 1;
%
% for cp = 1:nChannelPairs
%
%     scalars =[statMat{1,cp+1:nChannelPairs:endCol}]; %%% Peak Frequencies
%     weights = [statMat{12,cp+1:nChannelPairs:endCol}]; %%% Weights = Total Power
%     mu = circ_mean(scalars(:),weights(:)); %%% Weighted Mean
%     sig = circ_std(scalars(:),weights(:)); %%% Weighted Std
%     statMat{1,endCol+shift} = mu; %%% Mean of Peak Frequencies
%     statMat{1,endCol+shift+1} = sig; % Std of Peak Frequencies
%     intraIsoVars{2,shift} = masterVar.meanPeakFreq(cp);
%     intraIsoVars{2,shift+1} = masterVar.stdPeakFreq(cp);
%
%     scalars =[statMat{2,cp+1:nChannelPairs:endCol}]; %%% Mean Frequencies
%     mu = circ_mean(scalars(:),weights(:));
%     sig = circ_std(scalars(:),weights(:));
%     statMat{2,endCol+shift} = mu; %%%% Mean of Mean Frequencies
%     statMat{2,endCol+shift+1} = sig; %%% Std of Mean Frequencies
%     intraIsoVars{3,shift} = masterVar.meanFreq(cp);
%     intraIsoVars{3,shift+1} = masterVar.stdMeanFreq(cp);
%
%     scalars =[statMat{4,cp+1:nChannelPairs:endCol}]; %%% Peak Phases
%     scalars = scalars*(pi/180); %%% Converting to radians
%     mu = circ_mean(scalars(:),weights(:)); mu =  mu*(180/pi);
%     sig = circ_std(scalars(:),weights(:)); sig = sig*(180/pi);
%     statMat{4,endCol+shift} = mu; %%%% Mean of Peak Phases
%     statMat{4,endCol+shift+1} = sig; %%% Std of Peak Phases
%     intraIsoVars{5,shift} = masterVar.meanPeakPhase(cp);
%     intraIsoVars{5,shift+1} = masterVar.stdPeakPhase(cp);
%
%     scalars =[statMat{5,cp+1:nChannelPairs:endCol}]; %%% Mean Phases
%     scalars = scalars*(pi/180); %%% Converting to radians
%     mu = circ_mean(scalars(:),weights(:)); mu = mu*(180/pi);
%     sig = circ_std(scalars(:),weights(:)); sig = sig*(180/pi);
%     statMat{5,endCol+shift} = mu; %%%% Mean of Mean Phases
%     statMat{5,endCol+shift+1} = sig; %%% Std of Mean Phases
%     intraIsoVars{6,shift} = masterVar.meanPhase(cp);
%     intraIsoVars{6,shift+1} = masterVar.stdMeanPhase(cp);
%
%     scalars =[statMat{7,cp+1:nChannelPairs:endCol}]; %%% Mean Powers
%     mu = mean(scalars(:));
%     sig = std(scalars(:));
%     statMat{7,endCol+shift} = mu; %%%% Mean of Mean Powers
%     statMat{7,endCol+shift+1} = sig; %%% Std of Mean Powers
%     intraIsoVars{8,shift} = masterVar.meanPow(cp);
%
%     scalars =[statMat{12,cp+1:nChannelPairs:endCol}]; %%% Total Powers
%     mu = mean(scalars(:));
%     sig = std(scalars(:));
%     statMat{12,endCol+shift} = mu; %%%% Mean of Total Powers
%     statMat{12,endCol+shift+1} = sig; %%% Std of Mean Powers
%     intraIsoVars{13,shift} = masterVar.totPow(cp);
%
%     chLabelMat{1,endCol+shift} = ['Cross-file Avg ' num2str(ch(cp)) ' & ' num2str(ch(cp+1))];
%     chLabelMat{1,endCol+shift+1} = ['Cross-file STD ' num2str(ch(cp)) ' & ' num2str(ch(cp+1))];
%     fNamesMat{1,endCol+shift} = ['Cross-file Avg ' num2str(ch(cp)) ' & ' num2str(ch(cp+1))];
%     fNamesMat{1,endCol+shift+1} = ['Cross-file STD ' num2str(ch(cp)) ' & ' num2str(ch(cp+1))];
%     intraIsoVars{1,shift} = ['Intra Iso Cross-file Avg ' num2str(ch(cp)) ' & ' num2str(ch(cp+1))];
%     shift = shift + 2; %%% Since two cell cols are added in each loop
% end
% statMat =[chLabelMat; statMat; fNamesMat];
% statMat = [statMat intraIsoVars];
%
%
%
% % labCompSid = 'S-1-5-21-12604286-656692736-1848903544';
% % % myLaptopSid = 'S-1-5-21-2395549063-1654931228-1756539298'; % Dell E1505
% % myLaptopSid = 'S-1-5-21-3197369867-541179473-1092110829'; % Dell XPS
% % sid =getsid;
% % strangeCompCheck =0;
% % switch sid
% %     case labCompSid
% %         [success, message] = ...
% %             xlswrite('C:\Documents and Settings\pujalaav.nih\My Documents/temp.xls',statMat);
% %         if success
% %             'Data has been written to temp in My Documents'
% %         else
% %             errordlg('Data writing failed! Excel file must be closed for writing data')
% %         end
% %     case myLaptopSid
% %         [success,message]= ...
% %             xlswrite('C:\Users\Avi\Documents\temp.xls',statMat);
% %
% %         if success
% %             'Data has been written to temp in My Documents'
% %         else
% %             errordlg('Data writing failed! Excel file must be closed for writing data')
% %         end
% %     otherwise
% %         'This is not your lab or home computer, so data has not been written'
% %         strangeCompCheck = 1;
% % end
%
% %% OPTION FOR PHASE PLOTS
% % ButtonName = questdlg('Plot phase distributions?', ...
% %     'Phase distribution', ...
% %     'Yes', 'No', 'No');
% % switch ButtonName,
% %     case 'Yes',
% %         plotphase
% %     case 'No',
% % end %
%
%
% %% NEED FIXING --> The following lines of code need to be fixed to take multidimensional Wxy into acct
% % [mfvec,pfvec] = instantaneouswavefreq(Wxy_iso(:,:,1,1),freq);
% % tvpower = instantaneouswavepow(Wxy_iso(:,:,1,1));
% dynamicfreqpowplot
% plotphase
%
%
%
% %% Creating a Master Structure Variable that Saves All the Important Variables
%
% % answer = questdlg('Save Data?','Saving the Master Variable','No','Yes','No');
% answer = 'no';
% if strcmpi(answer,'Yes')
%     clear mName
% %     [master.Data, master.Time,master.Time_reduced, master.statMat] =...
% %         deal(dataStruct,timeAxisStruct,time_reduced,statMat);
% %
% %     [master.Wxy, master.Wxy_avg, master.meanPhaseVec, master.maxPhaseVec] =...
% %         deal(Wxy, Wxy_avg, meanPhaseVec, maxPhaseVec);
% %
% %     [master.samplingInt,master.meanPowVec, master.meanFreqVec, master.maxFreqVec] = ...
% %         deal(samplingInt, tvpower, mfvec, pfvec);
%
%     tm = repmat(time_reduced(:),1,nFiles);
%     tmat = sigMat;
%     tmat(:,:,3) = tm;
%     masterVar.Signals = tmat;
%     [masterVar.Wxy, masterVar.statMat] = deal(Wxy,statMat);
%     p = paths;
%     p(union(strfind(p,':'),strfind(p,'\')))='_';
%     mName = ['masterVariable_' p date];
%     save(mName,'masterVar')
% end
%
%
%
%
