% XWPLOTMD  Generates and plots crosswavelet spectra for pairs of signals
% This is a variant of XWPLOT in which data is stored as multidimensional
% arrays instead of as structures. File numbers change across the 3rd
% dimension and channel pair numbers vary across the 4th dimension.
% For e.g., Wxy(:,:,2,3) is the matrix of XW coefficients from the 2nd file
% and computed for the 3rd pair of channels being compared.

% Updated: 26-Aug-2014 (Modifications to allow for specification and overlaying of light channel on XWT)
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

wavelet_scale_resolution = 1/64; % (Must at least be 1/10)
number_of_phase_bins = 90; % Number of bins for phase histograms
motherWavelet = 'Morlet'; %%%('Morlet', 'Paul','DOG') - For now use only Morlet, other wavelets will give erroneous results
avgCheck =0;
peakDetectionThresh = 0.25;

%%  OPTIONS

Pad = 1; % Zero padding of signals (Pad = 1: Pad with zeroes; Pad = 0: Don't zero pad)
dj = wavelet_scale_resolution; % Wavelet scale resolution (Must be at least 1/10)
nPhaseBins = number_of_phase_bins; % Number of bins for phase histograms

time_axis_xticks = 'regular'; %('train' - displays the stimulus train, 'regular' - displays time regularly; default:'train')
figdisp = 'y'; %%% ('n' = does not display figures; [] = displays figs );


%% Wavelet Parameters - Derived and Fixed
fourier_factor = 1.0330; % This is the factor into which wavelet scales
% divide to yield frequency values. DO NOT CHANGE THIS VALUE (as long as
% you are using Morlet with wave# = 6)
scaleRange = 1./(freqRange*fourier_factor); % Scale range corresponding to frequency range.
S0 = min(scaleRange);
MaxScale = max(scaleRange);

%% XW Calculations
if freqRange(2) >= floor(1/samplingInt); errordlg('High frequency value above Nyquist limit'), end
% lpf = ceil(freqRange(2)/2.5);

%% Signal Reduction Parameters
newSamplingFrequency = max(freqRange(2)*2.5,20); % This ensures that the new sampling rate is well over twice the largest frequency.

%% Some processing
tempSig = eval(['signal' num2str(1) ';']);
sigMat = zeros(size(tempSig,1),nFiles,size(tempSig,2)); % Third dimension...
%... of sigMat contain channels, while second the files.
sigmas = zeros(size(sigMat,2),size(sigMat,3));
clear tempSig
for fileNum = 1:nFiles
    for chNum = 1:size(sigMat,3)
        eval(['sigMat(:,fileNum,chNum)= signal' num2str(fileNum) '(:,chNum);'])
        sigmas(fileNum,chNum) = std(sigMat(:,fileNum,chNum));
    end
end
%% Mean of Standard Deviation of all Signals
b = mean(sigmas,1);
c = repmat(b,size(sigmas,1),1);
sigmas_prenorm = sigmas;
sigmas = c;

%% Truncating the signal matrix and filtering
firstTime = time(1); lastTime = time(end);
time_reduced = time;


%% TO PLOT OR NOT TO PLOT
% plotfig  = questdlg('Would you like to plot time-varying frequencies for individual files?','To Plot or Not to Plot?','No','Yes','No');
plotfig = 'no';
%% Convert 3-D matrix of signals into a 2-D matrix
% Converts the 3-D matrix of signals (data points, files, channels) into a
% 2-D matrix of signals (data points, [(n file's channels of file1)
% (n + 1 file's channels)...(Nth file's channels)])

sigMat2 = zeros(size(sigMat,1),size(sigMat,3),size(sigMat,2));
for chNum = 1:size(sigMat,3)
    sigMat2(:,chNum,:) = sigMat(:,:,chNum);
end
sigMat2 =reshape(sigMat2,size(sigMat2,1),size(sigMat2,2)*size(sigMat2,3));
sigMat2 = [time_reduced(:) sigMat2];

%% SIGNAL REDUCTION PARAMETERS
sigMat = reducedata(sigMat,time_reduced,newSamplingFrequency); 
%% Normalizing Signals when Multiple Files are Loaded 
% This is so as to apply the same stringency (statistical threshold) to all signals. Important when computing \mathbb{S}^e.
% b = permute(sigmas,[3 1 2]);
% c = repmat(b,size(sigMat,1),1);
% sigMat = sigMat.*c;


%% Slow (Due to Light) Artifact Removal
% sigMat = slowartifactremove(sigMat) %%%% In case, I want to remove light
                                      %%%% artifacts from optogenetic
                                      %%%% trials.
time_reduced = linspace(time(1),time(end),size(sigMat,1));
dtr = time_reduced(2)-time_reduced(1);
endPt = time_reduced(1)+ (dtr*length(time_reduced));
time_reduced = time_reduced(1):dtr:endPt-dtr; %%% This step is necessary to ensure equal timesteps within time_reduced!!!
lenTime  = length(time_reduced);

%% STATISTICAL PARAMETERS
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

%% TRACE DISPLAY
if strcmpi(traceType,'smooth')|| strcmpi(traceType,'raw')
else
    errordlg('Input for the variable "traceType" not specified properly')
end
yShifter = 1.5; %(Factor by which the max value of a trace is multiplied...
% and shifted along the y-axis to prevent overlap; default = 1.5)

%% PHASE DISPLAY
if strcmpi(phaseType,'alt')|| strcmpi(phaseType,'synch') || strcmpi(phaseType,'all')
else
    errordlg('Input for the variable "phaseType" not specified properly')
end

%% POWER SPECTRUM DISPLAY
powerSpectrumType = 'both';  %('linear', 'log', 'both', 'none')
if strcmpi(powerSpectrumType,'linear')|| strcmpi(powerSpectrumType,'log') || strcmpi(powerSpectrumType,'both') || strcmpi(powerSpectrumType,'none')
else
    errordlg('Input for the variable "powerSpectrumType" not specified properly')
end
peakDetectionThresh = 0.1; % This determines the amplitude of a peak in the power spectrum for which a frequency value should be displayed. Setting to a low value detects low peaks and displace their corresponding freq values.

%% XWT
clear statMat chLabelMat fNamesMat
nChannelPairs = size(sigMat,3)-1;
nFiles = size(sigMat,2);
endCol = nFiles*nChannelPairs+1;
statMat = cell(15, nFiles*nChannelPairs + 1 + 2*nChannelPairs);
% statMat = cell(15, size(sigMat,2)*(size(sigMat,3)-1)+ 1 + 2*(size(ch,2)-1));
statMat(:,1)= deal({'Peak f','Mean f', 'Std f','Peak Ph', 'Mean Ph', 'Std Ph',...
    'Mean Pow','Std Pow','Synch Pow','Alt Pow','Alt Pow/Synch Pow' ,'Tot Pow','Tot Pow Ratio','Time Range','Freq Range'});
% chLabelMat =cell(1, size(sigMat,2)*(size(sigMat,3)-1)+1+2*(size(ch,2)-1));
chLabelMat =cell(1, nFiles*nChannelPairs + 1 + 2*nChannelPairs);
chLabelMat{1} = 'Channels';
% fNamesMat =cell(1, size(sigMat,2)*(size(sigMat,3)-1)+1+2*(size(ch,2)-1));
fNamesMat = cell(1,nFiles*nChannelPairs + 1 + 2*nChannelPairs);
fNamesMat{1}='File Names';
cellNum = 1;
time_varying_power_mat =[];
fileCounter = 0;
clear Wxy3d % 3D array with stacking Wxy matrices generated for each pair of channels from each file along the z-dimension
clear sigxy



lf = round(log2(MaxScale/S0)/dj)+10; % Aribitrarily chose 10
clear W W_coi W_coi_sig W_coi_sig_alt W_temp W_noncoi W_noncoi_sig
% W = zeros(lf,length(time_reduced),nFiles,length(ch));
% [W_coi,W_coi_sig,W_coi_sig_alt,W_temp] = deal(W);

for fileNum = 1:nFiles % File Number Loop # 1
    fStr =['f' num2str(fileNum)];
    fileCounter = fileCounter + 1;
    chNum =[];
    for chNum = 1:nChannelPairs % Channel Number Loop # 1
        chStr = ['ch' num2str(ch(chNum)) num2str(ch(chNum+1))];
        cellNum = cellNum+1;
        
        % Start: Normalizing time-series - Jun 06, 2012
        %                 sigMat(:,fileNum,chNum) = std(sigMat(:,fileNum,chNum))* normalizepdf(sigMat(:,fileNum,chNum));
        %                 sigMat(:,fileNum,chNum+1) = std(sigMat(:,fileNum,chNum+1))* normalizepdf(sigMat(:,fileNum,chNum+1));
        % End: Normalizing time-series
        
        [Wxy,period,scale,coi,sig95]= xwt([time_reduced(:) sigMat(:,fileNum,chNum)],[time_reduced(:) sigMat(:,fileNum,chNum+1)],...
            Pad, dj,'S0',S0, 'ms',MaxScale, 'Mother', motherWavelet);
        
        %         if fileNum == 1
        %             W(length(period)+1:lf,:,:,:)=[];
        %         end
        
        W(:,:,fileNum,chNum) = Wxy;
        
        %%%%% Normalization procedure so that the same stringency is
        %%%%% applied to all loaded signals
        xvar1 = sigmas(fileNum,chNum)* sigmas(fileNum,chNum+1);
        xvar2 = sigmas_prenorm(fileNum,chNum)* sigmas_prenorm(fileNum,chNum+1);
        multiplier = xvar2/xvar1;
        sig95(:,:,fileNum,chNum)= sig95*multiplier;
        
        % sig95(:,:,fileNum,chNum)= sig95; %%%% Uncomment this line and
        % comment the previous one if you want to establish significance
        % for each signal based on its own std, and not the group mean std.
        
        % eval(['sig95' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1) '= sig95;'])
        freq =1./period;
        ftmat = repmat(freq(:), 1, lenTime);
        coimat = repmat(1./coi(:)',length(freq), 1);
        
        W_temp = W(:,:,fileNum,chNum);       
        W_temp(ftmat<=coimat) = 0; %%% Removing regions outside of COI        
        W_coi(:,:,fileNum,chNum)= W_temp;
        W_noncoi(:,:,fileNum,chNum) = W(:,:,fileNum,chNum);
        
        %% SIGMA OR LEVEL-BASED THRESHOLDING
        if strcmpi(threshType,'sigma')
            W_temp(sig95(:,:,fileNum,chNum) < stringency) = 0; % (default: sig95 < 1)
            W_coi_sig(:,:,fileNum,chNum) = (W_temp);
            W_temp2 = W_noncoi(:,:,fileNum,chNum);
            W_temp2(sig95(:,:,fileNum,chNum) < stringency) = 0;
        elseif strcmpi(threshType,'level')
            W_temp(W_temp<level) = 0; % Keeps only those values
            % which are above a certain specified cross power 'level'
            W_coi_sig(:,:,fileNum,chNum) = W_temp;
        end
           Wxy = W_coi_sig(:,:,fileNum,chNum);
           W_noncoi_sig(:,:,fileNum,chNum) = W_temp2;
        
        %% PHASE-BASED FILTERING
        a  = angle(Wxy);
        if strcmpi(phaseType,'alt')
            Wxy(abs(a)<=(0.5*pi))=0; % Keeps only those matrix elements with angles >= (0.75*pi) = 135 deg
            W_coi_sig_alt(:,:,fileNum,chNum) = Wxy;  % abs(a)is necessary b/c angle(mat)outputs -ve values as well
            W_temp2(abs(a)<=(0.5*pi))=0; 
        elseif strcmpi(phaseType,'synch')
            Wxy(abs(a)>pi/2)=0; % Keeps only those matrix elements with angles <= pi/2
            W_coi_sig_alt(:,:,fileNum,chNum) = Wxy;  % abs(a)is necessary b/c angle(mat)outputs -ve values as well
            W_temp2(abs(a)>pi/2)=0;
        elseif strcmpi(phaseType,'all')
            W_coi_sig_alt(:,:,fileNum,chNum)= Wxy;
        else
            errordlg('PLEASE CHOOSE PROPER PHASETYPE!!!')
        end
        W_noncoi_sig(:,:,fileNum,chNum) = W_temp2;

        %% FREQUENCIES, PHASES, XW POWERS, ETC
        %       Wxy_bin =Wxy; Wxy_bin(Wxy~=0)=1; % Wxy_bin is a binary matrix created by substituting unity for all non-zero values in Wxy
        Wxy_prob = abs(Wxy./sum(Wxy(:))); % Each entry in Wxy is expressed as the probability of occurrence of that value within the distribution of all values of Wxy to create Wxy_prob
        if all(isnan(Wxy_prob(:)))
            Wxy_prob = zeros(size(Wxy_prob));
        end
        
        % Calculating mean & std of frequencies
        
        %         Wxy_s = sum(Wxy_prob,2);
        
        Fxy = ftmat; Fxy(Wxy==0)=0;
        Fxy_prob = Fxy.*Wxy_prob;
        meanFreqs.(fStr).(chStr) = round(sum(Fxy_prob(:))*100)/100;
        sigFreqs = Fxy;
        sigFreqs(Wxy==0)=[];
        stdFreqs.(fStr).(chStr) = round(std(sigFreqs)*100)/100;
        
        
        % Calculating XW power spectrum & Peak Frequency(frequency at which coherent power is highest)
        arbitraryDivisiveFactor = 1;
        globalPowSpec.(fStr).(chStr) = sum(abs(Wxy),2) * arbitraryDivisiveFactor;
        globalPowSpec_maxnorm.(fStr).(chStr) =...
            globalPowSpec.(fStr).(chStr)/...
            max(globalPowSpec.(fStr).(chStr)); % Power units indicate cumulative probability
        normlog_globalPowSpec.(fStr).(chStr) = ...
            log2(globalPowSpec.(fStr).(chStr));
        blah = normlog_globalPowSpec.(fStr).(chStr);
        blah(blah==-inf) = NaN;
        blah = (blah-min(blah))/(max(blah)-min(blah)); % Normalization by subtracting min and dividing by amplitude
        normlog_globalPowSpec.(fStr).(chStr) = blah;
        [maxtab,mintab] = peakdet(globalPowSpec_maxnorm.(fStr).(chStr),peakDetectionThresh);
        
        if numel(maxtab) == 0
            maxtab(:,1)= find(globalPowSpec_maxnorm.(fStr).(chStr) == max(globalPowSpec_maxnorm.(fStr).(chStr)));
            maxtab(:,2)= globalPowSpec_maxnorm.(fStr).(chStr)(maxtab(:,1));
        end
        
        pv = find(maxtab(:,2)== max(maxtab(:,2)));
        pf = round(freq(maxtab(pv,1))*100)/100;
        
        if isempty(pv) || isempty(pf)
            maxtab =[1 1+i];
            pv = 1;
            pf= 1+i;
        end
        
        % Calculating mean and std of phases within significant regions
        
        Wxy_nonzero_lin = Wxy; Wxy_nonzero_lin(Wxy==0)= []; % This removes all zero elements from the matrix and vectorizes it.
        if isempty(Wxy_nonzero_lin)
            Wxy_nonzero_lin = 0.001*(randn(1) + randn(1)*i);
        end
        Axy = angle(Wxy_nonzero_lin);
        
        nPhaseBins = min([numel(Axy), 90]); % Number of bins circular for phase histograms
        [ph_dist,th] = hist(Axy(:),nPhaseBins); % Unweighted phase histogram
        mag = abs(Wxy_nonzero_lin);
        [ph_dist3,vals] = hist3([Axy(:) mag(:)],[nPhaseBins,nPhaseBins]); % 3D bivariate (phases and power) histogram
        powmat = repmat(vals{2},size(ph_dist3,1),1);
        ph_dist_wt = ph_dist3.*powmat; % Scaling the number of elements in each bin by power (power-weighting)
        ph_dist_wt = sum(ph_dist_wt,2)'; % Power-weighted phase histogram
        
        mphase = angle(sum(Wxy_nonzero_lin)); % Mean phase is the angle of the resultant vector obtained by summing all the wavelet coefficients
        mphase(mphase<0) = mphase(mphase<0)+ 2*pi; % Addition of 2*pi to -ve values converts angle range from 0 to 360 rather than -180 to +180
        meanPhases.(fStr).(chStr) = mphase*180/pi; % Converts radians to degrees.
        meanPhases.(fStr).(chStr) = round(meanPhases.(fStr).(chStr)*100)/100;
        sphase = circ_std(Axy(:),abs(Wxy_nonzero_lin(:)));
        stdPhases.(fStr).(chStr) = sphase*180/pi;
        stdPhases.(fStr).(chStr) = round(stdPhases.(fStr).(chStr)*100)/100;
        
        if isempty(ph_dist)
            theta.(fStr).(chStr) = zeros(size(ph_dist_wt));
        else
            ph_dist = [ph_dist(:); ph_dist(1)]; % This will close the loop in polar plot by circularizing the vector.
            ph_dist_wt = [ ph_dist_wt(:); ph_dist_wt(1)];
            phase_dist.(fStr).(chStr) = ph_dist./max(ph_dist);
            phase_dist_weight.(fStr).(chStr) = ph_dist_wt./max(ph_dist_wt);
            theta.(fStr).(chStr) = [th(:);th(1)];
            phf = find(phase_dist_weight.(fStr).(chStr)==max(phase_dist_weight.(fStr).(chStr)));
            if numel(phf)~=0
                phf = phf(1);
            else phf = 1;
            end
        end
        peakPhase.(fStr).(chStr) = round(theta.(fStr).(chStr)(phf)*180/pi);
        if peakPhase.(fStr).(chStr)<0, peakPhase.(fStr).(chStr) = peakPhase.(fStr).(chStr)+360; end
        
        % Calculating mean and std of xw power in significant regions
        altPhases = find(angle(Wxy_nonzero_lin)>pi/2 | angle(Wxy_nonzero_lin)< -pi/2);
        synchPhases = find(angle(Wxy_nonzero_lin)<= pi/2 & angle(Wxy_nonzero_lin)>= -pi/2);
        sblah = abs(Wxy_nonzero_lin(synchPhases));
        synchPowers.(fStr).(chStr) = round(sum(sblah(:))*100)/100;
        ablah = abs(Wxy_nonzero_lin(altPhases));
        altPowers.(fStr).(chStr) = round(sum(ablah(:))*100)/100;
        meanPowers.(fStr).(chStr) = round(mean(abs(Wxy_nonzero_lin))*100)/100;
        stdPowers.(fStr).(chStr) = round(std(abs(Wxy_nonzero_lin))*100)/100;
        totalPowers.(fStr).(chStr) = round(sum(abs(Wxy_nonzero_lin)));
        altSynchPowRatio.(fStr).(chStr) = round((altPowers.(fStr).(chStr)/synchPowers.(fStr).(chStr))*100)/100;
        if round(totalPowers.(fStr).(chStr)-(synchPowers.(fStr).(chStr)+ altPowers.(fStr).(chStr)))>10
            errordlg('Synch Pow + Alt Pow ~= Tot Pow');
        end
        if cellNum==2
            firstPow = totalPowers.(fStr).(chStr);
        end
        totPowRatio.(fStr).(chStr) = round(100*totalPowers.(fStr).(chStr)/firstPow)/100;
        statMat(:,cellNum)= deal({pf, meanFreqs.(fStr).(chStr),...
            stdFreqs.(fStr).(chStr),peakPhase.(fStr).(chStr),meanPhases.(fStr).(chStr),...
            stdPhases.(fStr).(chStr),meanPowers.(fStr).(chStr),...
            stdPowers.(fStr).(chStr),synchPowers.(fStr).(chStr),...
            altPowers.(fStr).(chStr),altSynchPowRatio.(fStr).(chStr),totalPowers.(fStr).(chStr),...
            totPowRatio.(fStr).(chStr),num2str(timeRange),num2str(freqRange)});
        chLabelMat{cellNum} = ['f' num2str(fileNum) ' ch' num2str(ch(chNum))...
            num2str(ch(chNum+1))];
        fNamesMat{cellNum} = fNames(fileNum,:);
        
        
        
        %% PLOTTING FIGURES
        
        if lower(figdisp) =='y'
            
            figure('Name', ['XW Power, File ' num2str(fileNum) ', Channels '...
                num2str(ch(chNum)) ' vs ' num2str(ch(chNum+1))],'color','w','renderer','painter')
            figPos = get(gcf,'position');
            if strcmpi(powerSpectrumType,'none')
                set(gcf,'position',[figPos(1)*0.6 figPos(2)-figPos(4)/2 figPos(3)*1.3...
                    figPos(4)*2]);
            else
                set(gcf,'position',[figPos(1)*0.6 figPos(2)-figPos(4)/2 figPos(3)* 1.7...
                    figPos(4)*1.5]);
            end
            %% XWT AXES
            ax1 = axes; box off
            aPos = get(ax1,'position');
            if strcmpi(powerSpectrumType,'none')
                aPos = [aPos(1)*0.8 aPos(2)+ aPos(4)*(1/3) aPos(3) aPos(4)*0.75];
            else
                aPos = [aPos(1)*0.8 aPos(2)+ aPos(4)*(1/3) aPos(3)*(0.95) aPos(4)*0.75];
            end
            set(ax1,'position', aPos,'drawmode','fast')
        end
        
        
        
        sig95 = sig95(:,:,fileNum,chNum);
        Wxy = W_coi_sig_alt(:,:,fileNum,chNum);
        
        if lower(figdisp)=='y'
%             plotwave(Wxy,time_reduced,period,coi,sig95,...
%                 sigmas(fileNum,chNum), sigmas(fileNum,chNum+1))
sigmax = sigmas(fileNum,chNum);
sigmay = sigmas(fileNum,chNum+1);
          plotwave(W_temp2,time_reduced,period,coi,sig95,...
                sigmax, sigmay)
W_coi_sig_alt(:,:,fileNum,chNum) = W_coi_sig_alt(:,:,fileNum,chNum)/(sigmax*sigmay);
            set(ax1,'color','k','xtick',[], 'xcolor','w','ycolor','k','drawmode','fast')
            xlim([time_reduced(1) time_reduced(end)]) %%%% This line is NECESSARY to ensure that x-axis is aligned with traces below
            xlabel('')
            ylims1 = get(ax1,'ylim');
            yticklabels_ax1 = str2num(get(ax1,'yticklabel'));
            dyticklabels_ax1 = diff(yticklabels_ax1);
            %         if numel(dyticklabels_ax1)==1,yScale = 'log';
            %         elseif dyticklabels_ax1(1)== dyticklabels_ax1(2), yScale = 'linear';
            %         else yScale = 'log';
            %         end
            yScale = 'log';
            yl = ylabel('Frequency (Hz)','fontsize',14);
            ylpos = get(yl,'pos');
            aPos = get(ax1,'position');
            
            %% XW POWER SPECTRUM AXES
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
                    peakFreqs = freq(maxtab(:,1));
                end
                if strcmpi(yScale,'linear')
                    logPeakFreqs = peakFreqs;
                    freq2 = freq;
                else
                    logPeakFreqs = log2(peakFreqs);
                    freq2 = log2(freq);
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
                    % txt = {round(peakFreqs*100)/100};
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
                fstr = num2str(fileNum);
                tempSig = eval(['temp' num2str(fileNum) '(:,chNum);']);
                tempSig = truncatedata(tempSig,time,[firstTime lastTime]);
                
                tempTime = time; % Adding firstTime to the time vector here
                % will set the time of the first stimulus in the stimulus
                % train to a value of zero
                
                plot(tempTime,tempSig + (yShifter/2)*max(tempSig),'k','linewidth',1.5)
                tempSig = eval(['temp' num2str(fileNum) '(:,chNum+1);']);
                tempTime = linspace(firstTime,lastTime,length(tempSig));
                %                 tempSig = zscore(truncatedata(tempSig,time,[firstTime lastTime]));
                tempSig = truncatedata(tempSig,time,[firstTime lastTime]);
                plot(tempTime,tempSig-(yShifter/2)*max(tempSig),...
                    'k','linewidth',1.5)
                
            elseif strcmpi(traceType,'smooth')
                tempSig(:,chNum) = zscore(sigMat(:,fileNum,chNum));
                tempTime = linspace(firstTime,lastTime,length(tempSig)); % Adding firstTime to the time vector here
                % will set the time of the first stimulus in the stimulus
                % train to a value of zero
                %                 plot(tempTime,tempSig(:,chNum),colors(chNum),'linewidth',1.5)
                plot(tempTime,tempSig(:,chNum)+2.5,'k','linewidth',1.5)
                tempSig(:,chNum+1) = zscore(sigMat(:,fileNum,chNum+1));
                plot(tempTime,tempSig(:,chNum+1)-2.5,'k','linewidth',1.5)
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
        end
        
        
        %% TIME-VARYING MEAN FREQUENCIES AND XW POWERS
        
        [mfvec,pfvec] = instantaneouswavefreq(Wxy,freq);
        
        eval(['time_varying_meanfreq_f' num2str(fileNum) 'ch' num2str(chNum)...
            num2str(chNum+1) ' = mfvec;']);
        eval(['time_varying_pfreq_f' num2str(fileNum) 'ch' num2str(chNum)...
            num2str(chNum+1) ' = pfvec;']);
        tvpower = instantaneouswavepow(Wxy);
        
        eval(['time_varying_power_f' num2str(fileNum) 'ch' num2str(chNum)...
            num2str(chNum+1) ' = tvpower;']);
        %         eval(['tvpf_comb' num2str(fileNum) 'ch = tvpower;']);
        eval(['time_varying_power_mat = [time_varying_power_mat; time_varying_power_f' num2str(fileNum) 'ch' num2str(chNum)...
            num2str(chNum+1) '];']);
        
        %% OPTION FOR DYNAMIC FREQ AND PHASE PLOTS
        
        switch plotfig
            case 'Yes'
                dynamicfreqpowplot
            case 'No'
        end
        
        
        Wxy3d.(chStr)(:,:,fileNum) = [Wxy];
        sigxy.(chStr)(fileNum,:) = sigmas(fileNum,chNum)* sigmas(fileNum,chNum+1);
        
    end
end


%% Averaged XW Plots
 Wxy = W_coi_sig_alt;
[Wxy_avg,Wxy_iso,masterVar,W_gm,S,S_prelog,Wxy_avg_channels]  = avgxwt(Wxy,freq,time_reduced,coi,sigMat,isoThresh,ch);
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
% ButtonName = questdlg('Plot phase distributions?', ...
%     'Phase distribution', ...
%     'Yes', 'No', 'No');
% switch ButtonName,
%     case 'Yes',
%         plotphase
%     case 'No',
% end %


%% NEED FIXING --> The following lines of code need to be fixed to take multidimensional Wxy into acct
% [mfvec,pfvec] = instantaneouswavefreq(Wxy_iso(:,:,1,1),freq);
% tvpower = instantaneouswavepow(Wxy_iso(:,:,1,1));
% dynamicfreqpowplot
% plotphase



%% Creating a Master Structure Variable that Saves All the Important Variables

% answer = questdlg('Save Data?','Saving the Master Variable','No','Yes','No');
answer = 'no';
if strcmpi(answer,'Yes')
    clear mName
%     [master.Data, master.Time,master.Time_reduced, master.statMat] =...
%         deal(dataStruct,timeAxisStruct,time_reduced,statMat);
%     
%     [master.Wxy, master.Wxy_avg, master.meanPhaseVec, master.maxPhaseVec] =...
%         deal(Wxy, Wxy_avg, meanPhaseVec, maxPhaseVec);
%     
%     [master.samplingInt,master.meanPowVec, master.meanFreqVec, master.maxFreqVec] = ...
%         deal(samplingInt, tvpower, mfvec, pfvec);
    
    tm = repmat(time_reduced(:),1,nFiles);
    tmat = sigMat;
    tmat(:,:,3) = tm;
    masterVar.Signals = tmat;
    [masterVar.Wxy, masterVar.statMat] = deal(Wxy,statMat);
    p = paths;
    p(union(strfind(p,':'),strfind(p,'\')))='_';
    mName = ['masterVariable_' p date];
    save(mName,'masterVar')
end



%% Modifications log
%%%%% 06.19.13 - Created the Wxy3d structure array.

%% Pending fixes
% 1) Inscribing 'Power' on colorbar? - I can always do this using 'text'
% function
% 2) Ideally, the two signals need to be normalized using normalizepdf.m
%   provided by Grinsted if the histogram deviates severely from normal.
% 3) Verify logic for stdFreqs
% 4) Since sig regions are affected by length of signal it might be
% beneficial to select a region of the signals from which to calculate
% sigmax and sigmay


% Related Code
% STIMTRAIN - Generates a stim train
% PLOTPHASES - Plots of phase distributions