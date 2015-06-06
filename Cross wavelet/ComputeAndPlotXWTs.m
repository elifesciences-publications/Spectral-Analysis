% ComputeAndPlotXWTs  Generates and plots crosswavelet spectra for pairs of signals
% This is a variant of XWPLOT in which data is stored as multidimensional
% arrays instead of as structures. File numbers change across the 3rd
% dimension and channel pair numbers vary across the 4th dimension.
% For e.g., Wxy(:,:,2,3) is the matrix of XW coefficients from the 2nd file
% and computed for the 3rd pair of channels being compared.

% Author: AP

%% Variables Glossary
% Wxy ==> XWT coefficients. Values are complex numbers whose magnitude
%         cross wavelet power and angle gives the relative phase between
%         the two signals.
% sig95 ==> Matrix the size of Wxy containing numbers indicating how many
%           variances away the corresponding point in Wxy is from a 'red'
%           or 'white' noise background spectrum. For instance, if
%           sig95(4,5)= 3, that implies Wxy(4,5)is 3 variances
%           (95% confidence level) away from noise with the same mean as the
%           time series.
%

%% Wavelet Parameters
Wxy = struct;
Wxy.dj              = 1/2^6; % (Resolution of wavelet scales. Must be less than 1/10, e.g. 1/24)
Wxy.pad             = 1; % Zero padding of signals (wave.pad = 1: pad with zeroes; wave.pad = 0: don't zero pad)
nPhaseBins          = 90; % Number of bins for phase histograms
Wxy.motherWavelet   = 'Morlet'; %%%('Morlet', 'Paul','DOG') - For now use only Morlet, other wavelets will give erroneous results
avgCheck            = 0;
pkDetThr            = 0.25;
time_axis_xticks    = 'regular'; %('train' - displays the stimulus train, 'regular' - displays time regularly; default:'train')
figdisp             = 'y'; %%% ('n' = does not display figures; [] = displays figs );

%% Calculating fixed & derived wavelet parameters
Wxy.fourierFactor = 1.0330; % This is the factor into which wavelet scales
% divide to yield frequency values. DO NOT CHANGE THIS VALUE (as long as
% you are using Morlet with wave# = 6)
Wxy.scaleRange = 1./(freqRange*Wxy.fourierFactor); % Scale range corresponding to frequency range.
Wxy.S0 = min(Wxy.scaleRange);
Wxy.maxScale = max(Wxy.scaleRange);

if freqRange(2) >= floor(1/samplingInt); errordlg('High frequency value above Nyquist limit'), end
newSamplingFreq = max(freqRange(2)*2.5,20); % This ensures that the new sampling rate is well over twice the largest frequency.

% plotfig  = questdlg('Would you like to plot time-varying frequencies for individual files?','To Plot or Not to Plot?','No','Yes','No');
plotfig = 'no';

%% Group & normalize signals
sigMat = zeros(size(data(1).smooth,1),size(data(1).smooth,2),size(data,1)); % T-by-C-by-F matrix,
% where T = number of time points, C = num of channels, F = num
% of files
sigmaxy = zeros(size(sigMat,2),size(sigMat,3));
for file = 1:size(data,1)
    for chNum = 1:size(sigMat,2)
        sigMat(:,chNum,file) = data(file).smooth(:,chNum);
        sigmaxy(chNum,file) = std(sigMat(:,chNum,file));
    end
end
%# Put the signals from each channel end-to-end and compute global std for each channel
sigmaXY = std(reshape(permute(sigMat,[1 3 2]),size(sigMat,1)*size(sigMat,3),size(sigMat,2)),[],1);
scaleMat = permute(repmat(repmat(sigmaXY',1,size(sigMat,3))./sigmaxy,[1,1,size(sigMat,1)]),[3 1 2]);
%# Normalize signals such that global and local stds match
sigMat = sigMat.*scaleMat;

firstTime = data(1).time(1);
lastTime = data(1).time(end);

%% Subsample the Signals to Speed Computation
dt = floor((1/samplingInt)/newSamplingFreq);
sigMat = sigMat(1:dt:end,:,:);
Wxy.time = linspace(data(1).time(1),data(1).time(end),size(sigMat,1));
dtr = Wxy.time(2)-Wxy.time(1);
endPt = Wxy.time(1)+ (dtr*length(Wxy.time));
Wxy.time = Wxy.time(1):dtr:endPt-dtr; %%% This step is necessary to ensure equal timesteps within Wxy.time!!!
lenTime  = length(Wxy.time);

%% Slow Artifact Removal (optogenetic light signal arifact)
% sigMat = slowartifactremove(sigMat) %# In case, I want to remove light
%# artifacts from optogenetic
%# trials.

%% Statistical Parameters
level = 2^2; % (default: level = 2)
threshType = 'sigma'; %('sigma' - thresholds based on std; 'level' -
% thresholds based on level);
if isempty(stringency)
    stringency = 1; % I made a modification to the code in xwt.m which allows the user to select a region of the time series for which to calculate standard deviation, and therefore this "stringency" value actually reflects the # of standard deviations the activity must be to be considered significant.
end
if strcmpi(threshType,'level')|| strcmpi(threshType,'sigma')
else
    errordlg('Input for the variable "threshType" not specified properly')
end

%% Trace Display
if strcmpi(traceType,'smooth')|| strcmpi(traceType,'raw')
else
    errordlg('Input for the variable "traceType" not specified properly')
end
yShifter = 1.5; %(Factor by which the max value of a trace is multiplied...
% and shifted along the y-axis to prevent overlap; default = 1.5)

%% Phase Display
if strcmpi(phaseType,'alt')|| strcmpi(phaseType,'synch') || strcmpi(phaseType,'all')
else
    errordlg('Input for the variable "phaseType" not specified properly')
end

%% Power Spectrum Display
powerSpectrumType = 'both';  %('linear', 'log', 'both', 'none')
if strcmpi(powerSpectrumType,'linear')|| strcmpi(powerSpectrumType,'log') || strcmpi(powerSpectrumType,'both') || strcmpi(powerSpectrumType,'none')
else
    errordlg('Input for the variable "powerSpectrumType" not specified properly')
end

%% XWT
clear statMat chLabelMat fNamesMat
nChannelPairs = size(sigMat,3)-1;
nFiles = size(sigMat,3);
endCol = nFiles*nChannelPairs+1;
statMat = cell(15, nFiles*nChannelPairs + 1 + 2*nChannelPairs);
statMat(:,1)= deal({'Peak f','Mean f', 'Std f','Peak Ph', 'Mean Ph', 'Std Ph',...
    'Mean Pow','Std Pow','Synch Pow','Alt Pow','Alt Pow/Synch Pow' ,'Tot Pow','Tot Pow Ratio','Time Range','Freq Range'});
chLabelMat =cell(1, nFiles*nChannelPairs + 1 + 2*nChannelPairs);
chLabelMat{1} = 'Channels';
fNamesMat = cell(1,nFiles*nChannelPairs + 1 + 2*nChannelPairs);
fNamesMat{1}='File Names';
cellNum = 1;
time_varying_power_mat =[];
fileCounter = 0;
lf = round(log2(Wxy.maxScale/Wxy.S0)/Wxy.dj)+10; % Aribitrarily chose 10

W = [];
for file = 1:nFiles % File Number Loop # 1
    fStr =['f' num2str(file)];
    fileCounter = fileCounter + 1;
    for chNum = 1:nChannelPairs % Channel Number Loop # 1
        chStr = ['ch' num2str(chList(chNum)) num2str(chList(chNum+1))];
        cellNum = cellNum+1;
        
        [Wxy.raw(:,:,file,chNum),period,scale, coi, sig95]= xwt([Wxy.time(:) sigMat(:,chNum,file)],[Wxy.time(:) sigMat(:,chNum+1,file)],...
            Wxy.pad, Wxy.dj,'S0',Wxy.S0, 'ms', Wxy.maxScale, 'Mother', Wxy.motherWavelet);
       
        
        if file == 1
            Wxy.freq =1./period;
            ftmat = repmat(Wxy.freq(:), 1, lenTime);
            coimat = repmat(1./coi(:)',length(Wxy.freq), 1);            
        end
        
        temp = Wxy.raw(:,:,file,chNum);
        temp(ftmat<coimat) = 0; % Removing regions outside of COI
        Wxy.coi(:,:,file,chNum) = temp;
        
        %# Sigma or level based thresholding
        if strcmpi(threshType,'sigma')
            temp  = Wxy.raw(:,:,file,chNum);
            temp(sig95 < stringency) = 0;
            Wxy.sig(:,:,file,chNum) = temp;
            
            temp = Wxy.coi(:,:,file,chNum);
            temp(sig95 < stringency) = 0;
            Wxy.coi(:,:,file,chNum) = temp;
        elseif strcmpi(threshType,'level')
            temp = Wxy.raw(:,:,file,chNum);
            temp(temp < level) = 0;
            Wxy.sig(:,:,file,chNum) = temp;
            
            temp = Wxy.coi(:,:,file,chNum);
            temp(temp < level) = 0;
            Wxy.coi(:,:,file,ch) = temp;
        end
        
        %# Phase filtering       
        Axy  = angle(Wxy.raw(:,:,file,chNum));
        if strcmpi(phaseType,'alt')
            disp('Filtering out synchronous phases')
            temp = Wxy.sig(:,:,file,chNum);
            temp(abs(Axy)<=(0.5*pi)) = 0;
            Wxy.sig(:,:,file,chNum) = 0;
            
            temp = Wxy.coi(:,:,file,chNum);
            temp(abs(Axy)<=(0.5*pi)) = 0;
            Wxy.coi(:,:,file,chNum) = temp;
            
        elseif strcmpi(phaseType,'synch')
            disp('Filtering out alternating phases')
            temp = Wxy.sig(:,:,file,chNum);
            temp(abs(Axy)>pi/2) = 0;
            Wxy.sig(:,:,file,chNum) = temp;
            
            temp = Wxy.coi(:,:,file,chNum);
            temp(abs(Axy)>pi/2) = 0;
            Wxy.coi(:,:,file,chNum) = temp;
        elseif strcmpi(phaseType,'all')
            disp('No phase filtering')
        else
            errordlg('Phase filtering improperly specified!')
        end
        
        %# Estimating freq
        W = abs(Wxy.sig(:,:,file,chNum));
        sumFP = sum(ftmat(:).*W(:));
        sumP = sum(W(:));
        
        Wxy.meanFreq(file,chNum) = round((sum(ftmat(:).*W(:))/sumP)*100)/100;
        Wxy.stdFreq(file,chNum) = circ_std(ftmat(W~=0),W(W~=0));
        %         Wxy.freqWithMostPow(file,chNum) = round(ftmat(find(W==max(W(:))))*100)/100;
        Wxy.freqWithMostPow(file,chNum) = 1/period(find(sum(W,2)==max(sum(W,2))));
        
        %# Computing pow spectrum
        Wxy.powSpec{file,chNum} = mean(abs(W),2);
        blah = Wxy.powSpec{file,chNum};
        blah(blah==0) = nan;
        Wxy.powSpec_log{file,chNum} = log2(blah);
        [maxtab,~] = peakdet(Wxy.powSpec{file,chNum},pkDetThr*max(Wxy.powSpec{file,chNum}));
        
        if ~isempty(maxtab)
            pv = find(maxtab(:,2)== max(maxtab(:,2)));
            pf = round(Wxy.freq(maxtab(pv,1))*100)/100;
        else
            [maxtab(:,1),maxtab(:,2), pv, pf] = deal(nan);
        end
        
        %# Estimating phase
        W = Wxy.sig(:,:,file,chNum);
        W(W==0)=[];
        nPhaseBins = min([numel(angle(W)), 90]);
        [Wxy.phaseDist{file,chNum},theta] = hist(angle(W(:)),nPhaseBins); % Unweighted phase histogram
        [ph_dist3,vals] = hist3([angle(W(:)) abs(W(:))],[nPhaseBins,nPhaseBins]);
        powmat = repmat(vals{2},size(ph_dist3,1),1);
        Wxy.phaseDistWt{file,chNum} = ph_dist3.*powmat;
        Wxy.phaseDistWt{file,chNum} = sum(Wxy.phaseDistWt{file,chNum},2)'; % Power-weighted phase histogram
        mphase = angle(sum(W));
        mphase(mphase<0) = mphase(mphase<0)+ 2*pi; % 0 to 360 instead of -180 to +180 deg
        Wxy.meanPhase(file,chNum) = round(mphase*(180/pi)*100)/100;
        Wxy.stdPhase(file,chNum) = round(100*circ_std(angle(W(:)),abs(W(:)))*(180/pi))/100;
        
        if ~isempty(Wxy.phaseDist{file,chNum})
            Wxy.phaseDist{file,chNum} = [Wxy.phaseDist{file,chNum}(:); Wxy.phaseDist{file,chNum}(1)]; % Ligates the plot
            Wxy.phaseDist{file,chNum} = Wxy.phaseDist{file,chNum}/max(Wxy.phaseDist{file,chNum});
            
            Wxy.phaseDistWt{file,chNum} = [Wxy.phaseDistWt{file,chNum}(:); Wxy.phaseDistWt{file,chNum}(1)];
            Wxy.phaseDistWt{file,chNum} = Wxy.phaseDistWt{file,chNum}/max(Wxy.phaseDistWt{file,chNum});
            
            Wxy.theta(file,chNum) = {[theta(:); theta(1)]};
            phf = find(Wxy.phaseDistWt{file,chNum} == max(Wxy.phaseDistWt{file,chNum}));
            if ~isempty(phf)
                phf = phf(1);
                Wxy.phaseWithMostPow(file,chNum) = round(Wxy.theta{file,chNum}(phf)*180/pi);
            else
                Wxy.phaseWithMostPow(file,chNum) = nan;
            end
        end
        
        if Wxy.phaseWithMostPow(file,chNum) < 0,
            Wxy.phaseWithMostPow(file,chNum) = Wxy.phaseWithMostPow(file,chNum)+360;
        end
        
        %# Estimating power
        Wxy.totPow(file,chNum)  = round(sum(abs(W)));
        Wxy.totPow_alt(file,chNum) = round(sum(abs(W(angle(W) > pi/2 | angle(W) < -pi/2))));
        Wxy.totPow_synch(file,chNum) = round(sum(abs(W(angle(W) < pi/2 & angle(W) > -pi/2))));
        Wxy.meanPow_alt(file,chNum) = round(mean(abs(W(angle(W) > pi/2 | angle(W) < -pi/2))));
        Wxy.meanPow_synch(file,chNum) = round(mean(abs(W(angle(W) < pi/2 & angle(W) > -pi/2))));
        Wxy.stdPow_alt(file,chNum) = round(std(abs(W(angle(W) > pi/2 | angle(W) < -pi/2))));
        Wxy.stdPow_synch(file,chNum) = round(std(abs(W(angle(W) < pi/2 & angle(W) > -pi/2))));
        
        
        %# Plotting figures
        if strcmpi(figdisp,'y')
            figure('Name', ['XW Power, File ' num2str(file) ', Channels '...
                num2str(ch(chNum)) ' vs ' num2str(ch(chNum+1))],'color','w','renderer','painter')
            figPos = get(gcf,'position');
            if strcmpi(powerSpectrumType,'none')
                set(gcf,'position',[figPos(1)*0.6 figPos(2)-figPos(4)/2 figPos(3)*1.3...
                    figPos(4)*2]);
            else
                set(gcf,'position',[figPos(1)*0.6 figPos(2)-figPos(4)/2 figPos(3)* 1.7...
                    figPos(4)*1.5]);
            end
            ax1 = axes; box off
            aPos = get(ax1,'position');
            if strcmpi(powerSpectrumType,'none')
                aPos = [aPos(1)*0.8 aPos(2)+ aPos(4)*(1/3) aPos(3) aPos(4)*0.75];
            else
                aPos = [aPos(1)*0.8 aPos(2)+ aPos(4)*(1/3) aPos(3)*(0.95) aPos(4)*0.75];
            end
            set(ax1,'position', aPos,'drawmode','fast')
            sigmax = std(sigMat(:,chNum,file));
            sigmay = sigMat(chNum+1,file);
            plotwave(Wxy.sig(:,:,file,chNum),Wxy.time,period,coi,sig95, sigmax,sigmay)
            set(ax1,'color','k','xtick',[], 'xcolor','w','ycolor','k','drawmode','fast')
            xlim([Wxy.time(1) Wxy.time(end)]) %%%% This line is NECESSARY to ensure that x-axis is aligned with traces below
            xlabel('')
            ylims1 = get(ax1,'ylim');
            yticklabels_ax1 = str2num(get(ax1,'yticklabel'));
            dyticklabels_ax1 = diff(yticklabels_ax1);
            yScale = 'log';
            yl = ylabel('Frequency (Hz)','fontsize',14);
            ylpos = get(yl,'pos');
            aPos = get(ax1,'position');
            
            %# XW power spectrum axes
            if strcmpi(powerSpectrumType,'none')
            else
                ax2 = axes; hold on, box off
                aPos2 = get(ax2,'pos');
                aPos2 = [aPos(1) + aPos(3) aPos(2) aPos2(3)*(0.15) aPos(4)];
                set(ax2,'pos', aPos2, 'color','none','tickdir','out','fontsize',11,'drawmode','fast');
                xlabel([{'Normalized'}; {'Power'}])
                if imag(maxtab(:,2))
                    maxtab(:,2)= 1+i;
                    peakFreqs = 1+i;
                else
                    peakFreqs = Wxy.freq(maxtab(:,1));
                end
                if strcmpi(yScale,'linear')
                    logPeakFreqs = peakFreqs;
                    freq2 = Wxy.freq;
                else
                    logPeakFreqs = log2(peakFreqs);
                    freq2 = log2(Wxy.freq);
                end
                switch powerSpectrumType
                    case 'linear'
                        plot(norm_globalPowSpecow.(fStr).(chStr),freq2,'k','linewidth',2)
                        leg = {'Linear Spectrum'}
                    case 'log'
                        plot(normlog_globalPowSpecow.(fStr).(chStr),freq2,'k','linewidth',2)
                        set(ax2,'xtick',[],'xcolor','w','drawmode','fast')
                        leg = {'Log Spectrum'};
                        ax2b = axes('position',aPos2,'xaxislocation','bottom',...
                            'yaxislocation','right','color','none','xscale','log','ytick',[],'ycolor','k','fontsize',11);
                        xt = [0.25 0.5 1];
                        set(ax2b,'xtick',xt);
                        xlabel([{'Normalized'}; {'Power'}])
                    case 'both'
                        plot(globalPowSpec_maxnorm.(fStr).(chStr), freq2,'k','linewidth',2)
                        plot(normlog_globalPowSpec.(fStr).(chStr), freq2,'k:','linewidth',2,'parent',ax2)
                        leg ={'Linear'; 'Log'};
                        ax2b = axes('position',aPos2,'xaxislocation','top',...
                            'yaxislocation','right','tickdir','out','color','none','xscale','log','ytick',[],'ycolor','w','fontsize',11);
                        xt = [0.25 0.5 1];
                        set(ax2b,'xtick',xt,'xlim',[0 1]);
                        hold off
                end
                set(ax2,'ylim',ylims1,'xlim',[0 1],'xtick',[0.5 1],'ytick',[],'drawmode','fast')
                xvals = 0.35*ones(size(peakFreqs));
                clear txt
                for pf = 1:length(peakFreqs)
                    txt{pf} = [num2str(round(peakFreqs(pf)*100)/100) ' Hz'];
                end
                
                text(xvals,logPeakFreqs,txt,'fontsize',11,'color','r','parent',ax2);
                legend(ax2,leg,'fontsize',10) % This line needs to be here to legend can be moved by hand after fig is generated
                
            end
            
            %% TIME SERIES AXES
            ax3=axes; hold on, box off
            aPos3 = get(ax3,'position');
            aPos3 = [aPos(1) aPos3(2) aPos(3) aPos3(4)*(1/3)];
            set(ax3,'position',aPos3,'tickdir','out','color','w','ycolor','w','drawmode','fast');
            
            
            if strcmpi(traceType,'raw')
                fstr = num2str(file);
                tempSig = eval(['temp' num2str(file) '(:,chNum);']);
                tempSig = truncatedata(tempSig,time,[firstTime lastTime]);
                
                tempTime = time; % Adding firstTime to the time vector here
                % will set the time of the first stimulus in the stimulus
                % train to a value of zero
                
                plot(tempTime,tempSig + (yShifter/2)*max(tempSig),'k','linewidth',1.5)
                tempSig = eval(['temp' num2str(file) '(:,chNum+1);']);
                tempTime = linspace(firstTime,lastTime,length(tempSig));
                %                 tempSig = zscore(truncatedata(tempSig,time,[firstTime lastTime]));
                tempSig = truncatedata(tempSig,time,[firstTime lastTime]);
                plot(tempTime,tempSig-(yShifter/2)*max(tempSig),...
                    'k','linewidth',1.5)
                
            elseif strcmpi(traceType,'smooth')
                tempSig(:,ch) = zscore(sigMat(:,chNum,file));
                tempTime = linspace(firstTime,lastTime,length(tempSig)); % Adding firstTime to the time vector here
                % will set the time of the first stimulus in the stimulus
                % train to a value of zero
                %                 plot(tempTime,tempSig(:,chNum),colors(chNum),'linewidth',1.5)
                plot(tempTime,tempSig(:,ch)+2.5,'k','linewidth',1.5)
                tempSig(:,ch+1) = zscore(sigMat(:,chNum+1,file));
                plot(tempTime,tempSig(:,ch+1)-2.5,'k','linewidth',1.5)
            end
            yl2 = ylabel([{'Normalized'};{'Amplitude'}],'fontsize',14,'color','k');
            ylpos2 = get(yl2,'pos');
            ylpos2_mod = ylpos2;
            ylpos2_mod(1) = tempTime(1)-abs(ylpos(1))+ 0.65; ylpos2_mod(2)=0;
            set(yl2,'pos',ylpos2_mod);
            set(ax3,'ytick',[])
            axis([tempTime(1) tempTime(end) -inf inf])
            set(gca,'fontsize',14)
            switch lower(time_axis_xticks)
                case 'regular'
                    %                 xtick = fix(linspace(time(1),time(end),10));
                    %                 dxt = xtick(2)-xtick(1);
                    %                 xtick = xtick(1):dxt:xtick(end);
                    %                 set(ax3,'ytick',[],'xtick',xtick)
                    set(ax3,'ytick',[])
                case 'train'
                    stimtrain % Running this program generates 'tStimArts' a vector of stimulus artifact times
                    xlabel(['Stim Train (' num2str(stimDur) 'sec) @ ' num2str(stimFreq) ' Hz'],'fontsize',14)
                    set(ax3,'ytick',[],'xticklabel',[],'xtick',[tStimArts - tStimArts(1)])
                    hold off;
                    
            end
            
            
            %% TIME-VARYING MEAN FREQUENCIES AND XW POWERS
            
            [mfvec,pfvec] = instantaneouswavefreq(Wxy,Wxy.freq);
            
            eval(['time_varying_meanfreq_f' num2str(file) 'ch' num2str(ch)...
                num2str(ch+1) ' = mfvec;']);
            eval(['time_varying_pfreq_f' num2str(file) 'ch' num2str(ch)...
                num2str(ch+1) ' = pfvec;']);
            tvpower = instantaneouswavepow(Wxy);
            
            eval(['time_varying_power_f' num2str(file) 'ch' num2str(ch)...
                num2str(ch+1) ' = tvpower;']);
            %         eval(['tvpf_comb' num2str(fileNum) 'ch = tvpower;']);
            eval(['time_varying_power_mat = [time_varying_power_mat; time_varying_power_f' num2str(file) 'ch' num2str(ch)...
                num2str(ch+1) '];']);
            
            %% OPTION FOR DYNAMIC FREQ AND PHASE PLOTS
            
            switch plotfig
                case 'Yes'
                    dynamicfreqpowplot
                case 'No'
            end
            
            
            Wxy3d.(chStr)(:,:,file) = [Wxy];
            sigxy.(chStr)(file,:) = sigmaxy(file,ch)* sigmaxy(file,ch+1);
            
        end
    end
end


%% Averaged XW Plots
Wxy = W_coi_sig_alt;
[Wxy_avg,Wxy_iso,masterVar,W_gm,S,S_prelog,Wxy_avg_channels]  = avgxwt(Wxy,Wxy.freq,Wxy.time,coi,sigMat,isoThresh,ch);
avgCheck = 1;


%% Updating Statistics (statMat) and Storing in an Excel File
clear intraIsoVars
intraIsoVars = cell(size(statMat,1)+2,2*(size(ch,2)-1));
shift = 1;

for cp = 1:nChannelPairs
    
    scalars =[statMat{1,cp+1:nChannelPairs:endCol}]; %%% Peak Frequencies
    weights = [statMat{12,cp+1:nChannelPairs:endCol}]; %%% Weights = Total Power
    mu = circ_mean(scalars(:),weights(:)); %%% Weighted Mean
    sig = circ_std(scalars(:),weights(:)); %%% Weighted Std
    statMat{1,endCol+shift} = mu; %%% Mean of Peak Frequencies
    statMat{1,endCol+shift+1} = sig; % Std of Peak Frequencies
    intraIsoVars{2,shift} = masterVar.meanPeakFreq(cp);
    intraIsoVars{2,shift+1} = masterVar.stdPeakFreq(cp);
    
    scalars =[statMat{2,cp+1:nChannelPairs:endCol}]; %%% Mean Frequencies
    mu = circ_mean(scalars(:),weights(:));
    sig = circ_std(scalars(:),weights(:));
    statMat{2,endCol+shift} = mu; %%%% Mean of Mean Frequencies
    statMat{2,endCol+shift+1} = sig; %%% Std of Mean Frequencies
    intraIsoVars{3,shift} = masterVar.meanFreq(cp);
    intraIsoVars{3,shift+1} = masterVar.stdMeanFreq(cp);
    
    scalars =[statMat{4,cp+1:nChannelPairs:endCol}]; %%% Peak Phases
    scalars = scalars*(pi/180); %%% Converting to radians
    mu = circ_mean(scalars(:),weights(:)); mu =  mu*(180/pi);
    sig = circ_std(scalars(:),weights(:)); sig = sig*(180/pi);
    statMat{4,endCol+shift} = mu; %%%% Mean of Peak Phases
    statMat{4,endCol+shift+1} = sig; %%% Std of Peak Phases
    intraIsoVars{5,shift} = masterVar.meanPeakPhase(cp);
    intraIsoVars{5,shift+1} = masterVar.stdPeakPhase(cp);
    
    scalars =[statMat{5,cp+1:nChannelPairs:endCol}]; %%% Mean Phases
    scalars = scalars*(pi/180); %%% Converting to radians
    mu = circ_mean(scalars(:),weights(:)); mu = mu*(180/pi);
    sig = circ_std(scalars(:),weights(:)); sig = sig*(180/pi);
    statMat{5,endCol+shift} = mu; %%%% Mean of Mean Phases
    statMat{5,endCol+shift+1} = sig; %%% Std of Mean Phases
    intraIsoVars{6,shift} = masterVar.meanPhase(cp);
    intraIsoVars{6,shift+1} = masterVar.stdMeanPhase(cp);
    
    scalars =[statMat{7,cp+1:nChannelPairs:endCol}]; %%% Mean Powers
    mu = mean(scalars(:));
    sig = std(scalars(:));
    statMat{7,endCol+shift} = mu; %%%% Mean of Mean Powers
    statMat{7,endCol+shift+1} = sig; %%% Std of Mean Powers
    intraIsoVars{8,shift} = masterVar.meanPow(cp);
    
    scalars =[statMat{12,cp+1:nChannelPairs:endCol}]; %%% Total Powers
    mu = mean(scalars(:));
    sig = std(scalars(:));
    statMat{12,endCol+shift} = mu; %%%% Mean of Total Powers
    statMat{12,endCol+shift+1} = sig; %%% Std of Mean Powers
    intraIsoVars{13,shift} = masterVar.totPow(cp);
    
    chLabelMat{1,endCol+shift} = ['Cross-file Avg ' num2str(ch(cp)) ' & ' num2str(ch(cp+1))];
    chLabelMat{1,endCol+shift+1} = ['Cross-file STD ' num2str(ch(cp)) ' & ' num2str(ch(cp+1))];
    fNamesMat{1,endCol+shift} = ['Cross-file Avg ' num2str(ch(cp)) ' & ' num2str(ch(cp+1))];
    fNamesMat{1,endCol+shift+1} = ['Cross-file STD ' num2str(ch(cp)) ' & ' num2str(ch(cp+1))];
    intraIsoVars{1,shift} = ['Intra Iso Cross-file Avg ' num2str(ch(cp)) ' & ' num2str(ch(cp+1))];
    shift = shift + 2; %%% Since two cell cols are added in each loop
end
statMat =[chLabelMat; statMat; fNamesMat];
statMat = [statMat intraIsoVars];



% labCompSid = 'S-1-5-21-12604286-656692736-1848903544';
% % myLaptopSid = 'S-1-5-21-2395549063-1654931228-1756539298'; % Dell E1505
% myLaptopSid = 'S-1-5-21-3197369867-541179473-1092110829'; % Dell XPS
% sid =getsid;
% strangeCompCheck =0;
% switch sid
%     case labCompSid
%         [success, message] = ...
%             xlswrite('C:\Documents and Settings\pujalaav.nih\My Documents/temp.xls',statMat);
%         if success
%             'Data has been written to temp in My Documents'
%         else
%             errordlg('Data writing failed! Excel file must be closed for writing data')
%         end
%     case myLaptopSid
%         [success,message]= ...
%             xlswrite('C:\Users\Avi\Documents\temp.xls',statMat);
%
%         if success
%             'Data has been written to temp in My Documents'
%         else
%             errordlg('Data writing failed! Excel file must be closed for writing data')
%         end
%     otherwise
%         'This is not your lab or home computer, so data has not been written'
%         strangeCompCheck = 1;
% end

%% OPTION FOR PHASE PLOTS
ButtonName = questdlg('Plot phase distributions?', ...
    'Phase distribution', ...
    'Yes', 'No', 'No');
switch ButtonName,
    case 'Yes',
        plotphase
    case 'No',
end %


%% NEED FIXING --> The following lines of code need to be fixed to take multidimensional Wxy into acct
% [mfvec,pfvec] = instantaneouswavefreq(Wxy_iso(:,:,1,1),freq);
% tvpower = instantaneouswavepow(Wxy_iso(:,:,1,1));
% dynamicfreqpowplot
% plotphase



%% Creating a Master Structure Variable that Saves All the Important Variables

answer = questdlg('Save Data?','Saving the Master Variable','No','Yes','No');
answer = 'no';
if strcmpi(answer,'Yes')
    clear mName
    [master.Data, master.Time,master.Time_reduced, master.statMat] =...
        deal(dataStruct,timeAxisStruct,Wxy.time,statMat);
    
    [master.Wxy, master.Wxy_avg, master.meanPhaseVec, master.maxPhaseVec] =...
        deal(Wxy, Wxy_avg, meanPhaseVec, maxPhaseVec);
    
    [master.samplingInt,master.meanPowVec, master.meanFreqVec, master.maxFreqVec] = ...
        deal(samplingInt, tvpower, mfvec, pfvec);
    
    tm = repmat(Wxy.time(:),1,nFiles);
    tmat = sigMat;
    tmat(:,:,3) = tm;
    masterVar.Signals = tmat;
    [masterVar.Wxy, masterVar.statMat] = deal(Wxy,statMat);
    p = paths;
    p(union(strfind(p,':'),strfind(p,'\')))='_';
    mName = ['masterVariable_' p date];
    save(mName,'masterVar')
end




