% XWFREQSPECPLOT  Plots cross wavelet frequency spectrum
%   xwfreqspecplot plots the cross wavelet power spectrum of two signals
%   by collapsing the time axis. The figure axes are rotated and placed
%   next to the xwt plot generated by 'session5c_mod' with the frequency
%   axis running along the y-axis for ease of visualization. The program
%   offers a choice of plotting the power on a linear or log scale.


session5 % This program has to be run first to load up necessary variables into the workspace

%% Wavelet parameters
fourier_factor = 1.0330; % This is the factor into which wavelet scales
% divide to yield frequency values
freq_limits_for_xwt = [0.125 3];
scale_limits_for_xwt = 1./(freq_limits_for_xwt*fourier_factor);
Pad = 1;
dj = 1/64;
S0 = min(scale_limits_for_xwt);
MaxScale = max(scale_limits_for_xwt);

%% Statistical parameters
stringency =1; % (default: stringency = 1)
level = 2^2; % (default: level = 2)
threshType = 'sigma'; %('sigma' - thresholds based on std; 'level' -
% thresholds based on level);
if strcmpi(threshType,'level')|| strcmpi(threshType,'sigma')
else
    errordlg('Input for the variable "threshType" not specified properly')
end
%% Signal Reduction Parameters
newSamplingFrequency = 40;
sigMat = reducedata(sigMat,time_mod,newSamplingFrequency);
time_reduced = linspace(0,time_reduced(end)-time_reduced(1),size(sigMat,1));
lenTime  = length(time_reduced);

%% Plot Display Parameters
%%% Trace Display
traceType = 'smooth'; % ('raw' - displays the traces under the xwt as raw...
% data highpassed at 50Hz; 'smooth' - dispays as rectified, lowpassed
% traces)
if strcmpi(traceType,'smooth')|| strcmpi(traceType,'raw')
else
    errordlg('Input for the variable "traceType" not specified properly')
end
yShifter = 1.5; %(Factor by which the max value of a trace is multiplied...
% and shifted along the y-axis to prevent overlap; default = 1.5)

%%% Phase Display
phaseType = 'all'; % (default = 'alt'; can also be 'synch' or 'all' wherein the displayed phases are in synch or all respectively)
if strcmpi(phaseType,'alt')|| strcmpi(phaseType,'synch') || strcmpi(phaseType,'all')
else
    errordlg('Input for the variable "phaseType" not specified properly')
end

%%% Power Spectrum Display
powerSpectrumType = 'both';  %('linear', 'log', or 'both')
if strcmpi(powerSpectrumType,'linear')|| strcmpi(powerSpectrumType,'log') || strcmpi(powerSpectrumType,'both')
else
    errordlg('Input for the variable "powerSpectrumType" not specified properly')
end

%% XWT
statMat = cell(size(sigMat,2)*(size(sigMat,3)-1)+1,7);
statMat(1,:)= deal({'Mean f', 'Std f', 'Mean Ph', 'Std Ph',...
    'Mean Pow','Std Pow', 'Tot Pow'});
chLabelMat =cell(size(sigMat,2)*(size(sigMat,3)-1)+1,1);
chLabelMat{1} = 'Channels';
fNamesMat =cell(size(sigMat,2)*(size(sigMat,3)-1)+1,1);
fNamesMat{1}='File Names';
cellNum = 1;


fileCounter = 0;

for fileNum = 1:size(sigMat,2) % File Number Loop # 1
    fStr =['f' num2str(fileNum)];
    fileCounter = fileCounter + 1;
    for chNum =1:size(sigMat,3)-1 % Channel Number Loop # 1
        chStr = ['ch' num2str(ch(chNum)) num2str(ch(chNum+1))];
        cellNum = cellNum+1;
        eval(['[Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
            ',period,scale,coi,sig95]'...
            '= xwt([time_reduced(:) sigMat(:,fileNum,chNum)]'...
            ',[time_reduced(:) sigMat(:,fileNum,chNum+1)]'...
            ',Pad, dj,''S0'',S0, ''ms'',MaxScale, ''Mother'', ''Morlet'');'])
        eval(['sig95' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1) '= sig95;'])
        freq =1./period;
        ftmat = repmat(freq(:), 1, lenTime);
        coimat = repmat(1./coi(:)',length(freq), 1);
        eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
            '_coi = Wf' num2str(fileNum) 'ch' num2str(chNum)...
            num2str(chNum+1) ';'])
        eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
            '_coi(ftmat<coimat) = 0;'])
        eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
            '_coi_sig = Wf' num2str(fileNum) 'ch'  num2str(chNum)...
            num2str(chNum+1) '_coi;'])
        
        %% Sigma or Level Based Thresholding
        if strcmpi(threshType,'sigma')
            
            %                 eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
            %                     '_coi_sig(sig95<(mean(sig95(:)))+ stringency*std(sig95(:))) = 0;']) % This is really not how it is supposed to be done. I just did this when I didn't know any better. The real way to do this is to change the value of Zv in xwt.m - AP note
            
            eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
                '_coi_sig(sig95 < 1) = 0;'])
        elseif strcmpi(threshType,'level')
            Wxy = eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
                '_coi_sig;']);
            eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
                '_coi_sig(abs(Wxy)<level) = 0;'])% Keeps only those values
            % which are above a certain specified cross power 'level'
            
        end
        Wxy = eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
            '_coi_sig;']);
        %% Angle Based Filtering
        a  = angle(Wxy);
        switch phaseType
            case 'alt'
                Wxy(abs(a)<=pi/2)=0; % Keeps only those matrix elements with angles <= pi/2
                eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
                    '_coi_sig_alt = Wxy;']) % abs(a)is necessary b/c angle(mat)outputs -ve values as well
            case 'synch'
                Wxy(abs(a)>pi/2)=0; % Keeps only those matrix elements with angles <= pi/2
                eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
                    '_coi_sig_alt = Wxy;']) % abs(a)is necessary b/c angle(mat)outputs -ve values as well
            case 'all'
                eval(['Wf' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1)...
                    '_coi_sig_alt = Wxy;'])
        end
        
        %% Calculating XW power spectrum
        power_distribution_along_freq_axis.(fStr).(chStr) = sum(abs(Wxy),2);
        norm_power_distribution_along_freq_axis.(fStr).(chStr) =...
            power_distribution_along_freq_axis.(fStr).(chStr)/...
            max(power_distribution_along_freq_axis.(fStr).(chStr)); % Power units indicate cumulative probability
        normlog_power_distribution_along_freq_axis.(fStr).(chStr) = ...
            log2(power_distribution_along_freq_axis.(fStr).(chStr));
        normlog_power_distribution_along_freq_axis.(fStr).(chStr) = ...
            normlog_power_distribution_along_freq_axis.(fStr).(chStr)/max(normlog_power_distribution_along_freq_axis.(fStr).(chStr));
        
        
        %% Plotting Figures
        figure('Name', ['XW Power, File ' num2str(fileNum) ', Channels '...
            num2str(ch(chNum)) ' vs ' num2str(ch(chNum+1))],'color','w')
        figPos = get(gcf,'position');
        set(gcf,'position',[figPos(1) figPos(2)-figPos(4)/2 figPos(3)* 1.5...
            figPos(4)+figPos(4)/2]);
        
        %% XWT Axes
        ax1 = axes; box off
        aPos = get(ax1,'position');
        aPos = [aPos(1)-0.02 aPos(2)+ aPos(4)*(1/3) aPos(3)*(0.8) aPos(4)*(2/3)];
        set(ax1,'position', aPos)
        sig95 = eval(['sig95' num2str(fileNum) 'ch' num2str(chNum) num2str(chNum+1) '']);
        Wxy = eval(['Wf' num2str(fileNum) 'ch' num2str(chNum)...
            num2str(chNum+1) '_coi_sig_alt']);
        plotwave(Wxy,time_reduced,period,coi,sig95,...
            sigmas(fileNum,chNum), sigmas(fileNum,chNum+1))
        set(ax1,'color','k','xtick',[], 'xcolor','w','ycolor','k')
        xlim([time_reduced(1) time_reduced(end)]) %%%% This line is NECESSARY to ensure that x-axis is aligned with traces below
        xlabel('')
        
        yl = ylabel('Frequency (Hz)','fontsize',14);
        ylpos = get(yl,'pos');
        aPos = get(ax1,'position');
        
        %% XW Power Spectrum Axes
        ax2 = axes; hold on, box off
        aPos2 = get(ax2,'pos');
        aPos2 = [aPos(1) + aPos(3) aPos(2) aPos2(3)*(0.3) (aPos2(4)*(2/3))];% The -0.0025 is necessitated by a quirk in the way matlab generates figures
        set(ax2,'pos', aPos2, 'color','w','tickdir','out','fontsize',14);
        xlabel([{'Normalized'}; {'Power'}])
        switch powerSpectrumType
            case 'linear'
                plot(norm_power_distribution_along_freq_axis.(fStr).(chStr),log2(freq),'k','linewidth',2)
                legend('Linear', 'fontsize',10)
            case 'log'
                plot(normlog_power_distribution_along_freq_axis.(fStr).(chStr),log2(freq),'k','linewidth',2)
                legend('Log','fontsize',10)
            case'both'
                plot(norm_power_distribution_along_freq_axis.(fStr).(chStr),log2(freq),'k','linewidth',2)
                hold on
                plot(normlog_power_distribution_along_freq_axis.(fStr).(chStr),log2(freq),'k:','linewidth',2)
                legend('Linear', 'Log')
        end
        ylim([-inf inf])
        yt = get(ax2,'ytick');
        ytl = num2str(2.^yt');
        set(ax2,'xtick',[0.5 1],'ytick',[]); % One tick at half-way pt and the other at peak value of 1 since power is normalized.
        maxPowFreq = freq(power_distribution_along_freq_axis.(fStr).(chStr)==max(power_distribution_along_freq_axis.(fStr).(chStr)));
        logMaxPowFreq = log2(maxPowFreq);
        roundmpf = round(maxPowFreq*100)/100;
        text(0.35,logMaxPowFreq,[num2str(roundmpf) 'Hz'],'fontsize',14,'color','r')
        
        
        %% Time Series Axes
        ax3=axes; hold on, box off
        aPos3 = get(ax3,'position');
        aPos3 = [aPos(1) aPos3(2) aPos(3) aPos3(4)*(1/3)];
        set(ax3,'position',aPos3,'tickdir','out','color','w','ycolor','w');
        switch traceType
            case 'raw'
                fstr = num2str(fileNum);
                tempSig = eval(['temp' num2str(fileNum) '(:,chNum);']);
                tempSig = zscore(truncatedata(tempSig,time,[firstTime lastTime]));
                
                tempTime = linspace(0,time_reduced(end),length(tempSig))...
                    + firstTime; % Adding firstTime to the time vector here
                % will set the time of the first stimulus in the stimulus
                % train to a value of zero
                
                plot(tempTime,tempSig + (yShifter/1.5)*max(tempSig),'k','linewidth',1.5)
                tempSig = eval(['temp' num2str(fileNum) '(:,chNum+1);']);
                tempSig = zscore(truncatedata(tempSig,time,...
                    [firstTime lastTime]));
                plot(tempTime,tempSig-(yShifter/1.5)*max(tempSig),...
                    'k','linewidth',1.5)
                
            case 'smooth'
                tempSig(:,chNum) = zscore(sigMat(:,fileNum,chNum));
                tempTime = linspace(0,time_reduced(end),length(tempSig))...
                    + firstTime; % Adding firstTime to the time vector here
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
        %                     xticks = get(ax3,'xtick');
        set(ax3,'ytick',[])
        axis([tempTime(1) tempTime(end) -inf inf])
        xlabel('Time (s)','fontsize',14)
        set(gca,'fontsize',14)
        set(ax3,'ytick',[])
        hold off;
    end
end

%% Bug Fixes
% 1) XWT plot time axis was misaligned by a slight bit with Time Series
% plot time axes. Fixed by setting xlims for each axis
% 2) Found a way to adjust position of ylabels so that they align for
% vertically juxtaposed axes (2.11.11)


%% Pending fixes
% 1) Sig95? How does this work?
% 2) Inscribing 'Power' on colorbar.
