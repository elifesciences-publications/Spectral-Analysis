

 function [Wxy_avg, Wxy_iso,masterVar] = avgxwt(Wxy,freq,time,coi,sigMat,isoThresh,ch)
% [Wxy_avg,Wxy_iso,statMat] = avgxwt(Wxy,freq,time,coi,sigMat,isoThresh,ch)
% Size of Wxy: m by n ny f by c, where m = number of frequency values, n =
% number of time points, and f = number of files, c = number of channel
% pairs

arrowDensity =[18 18]; % default [25 25]
ad = mean(arrowDensity);
arrowSize = 1*30.*0.03/ad;
% arrowSize = 0; % Uncommenting this line and commenting the prev one will
% cause arrows not to display
arrowHeadSize=0.9*arrowSize*180; % default 1*arrowsize*220
peakDetectionThresh = 0.03;

if isempty(isoThresh)
    isoThresh = 0.5;
end
period = 1./freq;
coi = 1./coi;
masterVar = [];
nFreqPts = size(Wxy,1);
nTimePts = size(Wxy,2);
nFiles = size(Wxy,3);
nChannelPairs = size(Wxy,4);

dt = time(2)-time(1);
freq = flipud(sort(freq(:)))';
signals = permute(sigMat,[3,1,2]); % Row # = Channel #, Cols # = Time Point
% 3rd dim # = File Num
% sigmaxy = std(signals,0,2);

%% Extracting Crossvariance for Each Pair of Signals
crossVar = [];
for fn = 1:nFiles
    for cp = 1:nChannelPairs
        crossVar(1,1,fn,cp) = std(signals(cp,:,fn))*std(signals(cp+1,:,fn));
    end
end
crossVarMat = repmat(crossVar,nFreqPts,nTimePts);

Wxy_norm = Wxy./crossVarMat; % Normalizing Wxy by the crossvariance of the underlying signal pair


%% Averaging Across Files, Channels, and Both
clear Wxy_avg_files Wxy_avg_channels Wxy_avg
Wxy_avg_files(:,:,1,:) = mean(Wxy_norm,3);
Wxy_avg_channels(:,:,:,1) = mean(Wxy_norm,4);
Wxy_avg = mean(mean(Wxy_norm,4),3);

clear tvmf tvpf tvpow
for cp = 1:nChannelPairs
    Wxy_temp = Wxy_avg_files(:,:,1,cp);
    absPow = abs(Wxy_temp);
    absPow(absPow==0) = []; % Ignoring non-significant regions
    [pCount,pVal] = hist(absPow,100);
    pProb = pCount/sum(pCount);
    pCdf = cumsum(pProb);
    cutoff = min(find(pCdf>=isoThresh));
    pCutoff = log2(pVal(cutoff));
    
    globPowSpec(:,cp) = mean(abs(Wxy_temp),2);
    amp1 = max(globPowSpec(:,cp))- min(globPowSpec(:,cp));
    normGlobPowSpec(:,cp) = globPowSpec(:,cp)/amp1;
    logGlobPowSpec(:,cp) = log2(globPowSpec(:,cp));
    blah = logGlobPowSpec(:,cp);
    blah(blah==-inf)= NaN; % log2(0) = -Inf
    amp2 = max(blah) - min(blah);
    isoThresh_lin = (2^pCutoff)/amp1;
    isoThresh_log = (pCutoff-min(blah))/amp2;
    blah  = (blah-min(blah))/amp2; % Normalization by min subtraction and amp division
    normLogGlobPowSpec(:,cp) = blah;
    
    Wxy_temp(abs(Wxy_temp)<(2^pCutoff))=0;
    Wxy_iso(:,:,1,cp) = Wxy_temp;
    
    Wxy_prob = abs(Wxy_temp)/sum(abs(Wxy_temp(:)));
    freqMat = repmat(freq(:),1,size(Wxy_temp,2));
    powWtedFreqMat = freqMat.*Wxy_prob;
%     meanFreq(cp,:) = sum(powWtedFreqMat(:));  %%% This is one way to
%     estimate mean frequency
    meanFreq(cp,:) = circ_mean(freqMat(:),abs(Wxy_temp(:)));
    stdMeanFreq(cp,:) = circ_std(freqMat(:),abs(Wxy_temp(:)));
    
    [tvmf(cp,:),tvpf(cp,:)] = instantaneouswavefreq(Wxy_temp, freq);
    tvpow(cp,:)= mean(abs(Wxy_temp),1);
    tvPeakPow(cp,:) = max(abs(Wxy_temp),[],1);
    tvpow_prob = tvpow(cp,:);
    tvpow_prob = tvpow_prob/sum(tvpow_prob(:));
%     meanPeakFreq(cp,:) = sum(tvpow_prob.*tvpf(cp,:));
    meanPeakFreq(cp,:) = circ_mean(tvpf(cp,:)',tvPeakPow(cp,:)');
    stdPeakFreq(cp,:) = circ_std(tvpf(cp,:)',tvPeakPow(cp,:)');
    
    meanPhase(cp,:) = angle(sum(Wxy_temp(:)))*(180/pi);
    meanPhase(cp,:)  = round(meanPhase(cp,:)*100)/100; % Rounding to 2 decimal places
    stdMeanPhase(cp,:) = circ_std(angle(Wxy_temp(:))*(180/pi),abs(Wxy_temp(:)));
    
    [tvmph(cp,:),tvpph(cp,:),stdph(cp,:)] = instantaneouswavephase(Wxy_temp);
%     meanPeakPhase(cp,:) = sum(tvpph(cp,:).*tvpow_prob);
    meanPeakPhase(cp,:) = circ_mean(tvpph(cp,:)',tvPeakPow(cp,:)');
    stdPeakPhase(cp,:) = circ_std(tvpph(cp,:)',tvPeakPow(cp,:)');
    
    totPow(cp,:) = sum(abs(Wxy_temp(:)));
    nonZeroPow = abs(Wxy_temp(:));
    nonZeroPow(nonZeroPow==0)=[];
    meanPow(cp,:) = mean(nonZeroPow);
    
    [maxtab,mintab] = peakdet(normGlobPowSpec(:,cp),peakDetectionThresh);
    if numel(maxtab) == 0
        maxtab(:,1)= find(normGlobPowSpec(:,cp) == max(normGlobPowSpec(:,cp)));
        maxtab(:,2)= normGlobPowSpec(maxtab(:,1),cp);
    end
    
    pv = find(maxtab(:,2)== max(maxtab(:,2)));
    pf = round(freq(maxtab(pv,1))*100)/100;
    
    if isempty(pv) || isempty(pf)
        maxtab =[1 1+i];
        pv = 1;
        pf= 1+i;
    end
    
    CData = Wxy_avg_files(:,:,1,cp); CData(CData==0)=nan;
    monSize = getMonitorSize;
    figName = [' XW Spectrum Averaged Across Files for Ch ' num2str(ch(cp)) ' vs ' num2str(ch(cp+1)) ];
    fh = figure('Name', figName,'color','w');
  
    set(fh,'position',[monSize(1)+50 monSize(2)+150 monSize(3)-150 monSize(4)-300]);
    %% XW Axes
    ax1 = gca; box off
    aPos = get(ax1,'position');
    aPos = [aPos(1)*0.8 aPos(2) aPos(3)*(0.95) aPos(4)];
    set(ax1,'position', aPos)
    
    Yticks = 2.^(fix(log2(min(freq))):fix(log2(max(freq))));
    H=imagesc(time,log2(freq),log2(abs(CData)));
    
    colormap(jet)
    clim=get(gca,'clim'); % center color limits around log2(1)=0
    
    clim = [-9 9]; % default= [-9 9]
    clim = fix(clim);
    set(gca,'clim',clim);
    hold on
    contour(time,log2(freq),log2(abs(CData)),[pCutoff pCutoff],'m','linewidth',2);
    % ch = contour(time(:),log2(freq(:)),log2(abs(CData)),[pCutoff pCutoff]);
    % set(ch,'color','k','linestyle',':','linewidth',2)
    
    figName = ['$$' figName];
   eqtxt = ' : {1\over n} \Sigma_{file = 1}^{n} |W_{ch';
%    eqtxt = ' : {1\over n} \sum\limits{file = 1}^{n} |W_{ch';
    eqtxt = [ eqtxt num2str(ch(cp)) 'ch' num2str(ch(cp+1))];  
    eqtxt = [eqtxt '}(file)| $$'];
  figName = [figName ' ' eqtxt];
   title(figName,'interpreter','latex','fontsize', 14)
    
    %% Y-Axis Adjustments
    set(gca,'YLim',log2([min(freq),max(freq)]), ...
        'YDir','normal', ...
        'YTick',log2(Yticks(:)), ...
        'YTickLabel',num2str(Yticks'), ...
        'layer','top','fontsize',14)
    ylabel('$$ Frequency (Hz) $$','interpreter','latex','fontsize',14)
    xlabel('$$ Time (sec) $$','interpreter','latex','fontsize', 14)
    hold on
    
    %% PHASE PLOTTING
    aWxy=angle(CData);
    phs_dt=round(length(time)/arrowDensity(1));
    tidx=max(floor(phs_dt/2),1):phs_dt:length(time);
    phs_dp=round(length(freq)/arrowDensity(2));
    pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
    phaseplot(time(tidx),log2(freq(pidx)),aWxy(pidx,tidx),arrowSize,arrowHeadSize);
    time = time(:);
    % freq = fliplr(freq(:)');
    %% COI
    tt=[time([1 1])-dt*.5;time;time([end end])+dt*.5];
    hcoi=plot(tt,log2([freq([end 1]) coi freq([1 end])]),'w','linewidth',1.5);
    hold off
    
    isoThreshLine_lin = isoThresh_lin*ones(size(freq));
    isoThreshLine_log = isoThresh_log*ones(size(freq));
    ax2 = axes;
    aPos2 = get(ax2,'position');
    aPos2 = [aPos(1) + aPos(3) aPos(2) aPos2(3)*(0.15) aPos(4)];
    set(ax2,'position', aPos2, 'color','none','tickdir','out','fontsize',11);
    xlabel([{'Normalized'}; {'Power'}])
    
    if imag(maxtab(:,2))
        maxtab(:,2)= 1+i;
        peakFreqs = 1+i;
    else
        peakFreqs = freq(maxtab(:,1));
    end
    
    plot(normGlobPowSpec(:,cp), log2(freq(:)),'k','linewidth',2), box off
    hold on
    plot(normLogGlobPowSpec(:,cp), log2(freq(:)),'k:','linewidth',2,'parent',ax2)
    plot(isoThreshLine_lin,log2(freq),'m:','linewidth',2,'parent',ax2)
    plot(isoThreshLine_log,log2(freq),'m:','linewidth',2,'parent',ax2)
    set(ax2,'xlim',[0 1])
    leg ={'Linear'; 'Log'; 'Isoline Threshold'};
    ax2b = axes('position',aPos2,'xaxislocation','top',...
        'yaxislocation','right','tickdir','out','color','none','xscale','log','ytick',[],'ycolor','w','fontsize',11);
    set(ax2b,'xtick',[0.25 0.5 1])
    hold off
    set(ax2,'YLim',log2([min(freq),max(freq)]), ...
        'YDir','normal', ...
        'YTick',[], ...
        'layer','top','fontsize',12,'xlim',[0 1],'xtick',[0.5 1])
    box off
    
    xvals = 0.35*ones(size(peakFreqs));
    clear txt
    for pf = 1:length(peakFreqs)
        txt{pf} = [num2str(round(peakFreqs(pf)*100)/100) ' Hz'];
        % txt = {round(peakFreqs*100)/100};
    end
    logPeakFreqs = log2(peakFreqs);
    text(xvals,logPeakFreqs,txt,'fontsize',11,'color','r','parent',ax2);
    
    legend(ax2,leg,'fontsize',10) % This line needs to be here to legend can be moved by hand after fig is generated
   
    %%%%% Plotting Dynamic Variables %%%%% 
    dynamicfreqpowplot
    plotphase
end

[masterVar.tvmf, masterVar.tvpf, masterVar.tvpow, masterVar.tvmph, masterVar.tvpph] = ...
    deal(tvmf,tvpf,tvpow,tvmph,tvpph);
[masterVar.meanFreq,masterVar.stdMeanFreq, masterVar.meanPeakFreq,masterVar.stdPeakFreq,...
    masterVar.meanPhase,masterVar.stdMeanPhase, masterVar.meanPeakPhase, masterVar.stdPeakPhase] = ...
    deal(meanFreq,stdMeanFreq, meanPeakFreq,stdPeakFreq, meanPhase,stdMeanPhase, meanPeakPhase,stdPeakPhase);
[masterVar.totPow, masterVar.meanPow] = deal(totPow, meanPow);


% absPow = abs(Wxy_avg(:));
% absPow(absPow==0) = []; % Ignoring non-significant regions
% [pCount,pVal] = hist(absPow,100);
% pProb = pCount/sum(pCount);
% pCdf = cumsum(pProb);
% cutoff = min(find(pCdf>=isoThresh));
% pCutoff = log2(pVal(cutoff));
%
% globPowSpec_avg = sum(abs(Wxy_avg),2)/length(time);
% normGlobPowSpec_avg = globPowSpec_avg/(max(globPowSpec_avg));
% logGlobPowSpec_avg = log2(globPowSpec_avg);
% blah = logGlobPowSpec_avg;
% blah(blah==-inf)= NaN; % log2(0) = -Inf
% isoThresh_lin = (2^pCutoff);
% isoThresh_lin = (isoThresh_lin)/max(globPowSpec_avg);
% isoThresh_log = (pCutoff-min(blah))/(max(blah)-min(blah));
% blah  = (blah-min(blah))/(max(blah)-min(blah)); % Normalization by min subtraction and amp division
% normlogGlobPowSpec_avg = blah;
% % sNormGlobPowSpec_avg = globPowSpec_avg(:)./log2(period(:));
%
% %% Creating a Wxy_iso: Wxy_iso(m,n) = {0, if Wxy_iso(m,n) < Isoline Threshold
% Wxy_iso = Wxy_avg;
% Wxy_iso(abs(Wxy_iso)<(2^pCutoff))=0;
%
% %% Creating a Statistics Matrix
% variableNames = {'Mean Freq'; 'Mean of Peak Freqs'; 'Mean Phase'; 'Mean of Peak Phases'; 'Total Pow'; 'Mean Pow'};
%
% %% Calculating Mean Freq & Phase from Wxy_iso
% Wxy_prob = abs(Wxy_avg)/sum(abs(Wxy_avg(:)));
% freqMat = repmat(freq(:),1,size(Wxy_iso,2));
% powWtedFreqMat = freqMat.*Wxy_prob;
% meanFreq_avg = sum(powWtedFreqMat(:));
% statMat{1,1} = meanFreq_avg;
%
% [mfvec,pfvec] = instantaneouswavefreq(Wxy_iso, freq);
% tvpow_prob = sum(abs(Wxy_iso));
% tvpow_prob = tvpow_prob/sum(tvpow_prob);
% meanFreq_avg = sum(tvpow_prob.*mfvec);
% peakFreq_avg = sum(tvpow_prob.*pfvec);
%
% statMat{1,2} = meanFreq_avg;
% statMat{2,1} = peakFreq_avg;
%
% meanPhase = angle(sum(Wxy_iso(:)))*(180/pi);
% meanPhase  = round(meanPhase*100)/100; % Rounding to 2 decimal places
%
% [meanPhaseVec,maxPhaseVec,stdPhaseVec] = instantaneouswavephase(Wxy_iso);
% meanPeakPhase = sum(maxPhaseVec.*tvpow_prob);
%
% statMat{3,1} = meanPhase;
% statMat{4,1} = meanPeakPhase;
%
% totPow = sum(abs(Wxy_iso(:)));
% nonZeroPow = abs(Wxy_iso(:));
% nonZeroPow(nonZeroPow==0)=[];
% meanPow = mean(nonZeroPow);
%
% statMat{5,1} = totPow;
% statMat{6,1} = meanPow;
%
% statMat = [variableNames, statMat];
%
%
%% Determining Peak Frequency from the Averaged XW Global Power Spectrum
% [maxtab,mintab] = peakdet(normGlobPowSpec_avg,peakDetectionThresh);
% if numel(maxtab) == 0
%     maxtab(:,1)= find(normGlobPowSpec_avg == max(normGlobPowSpec_avg));
%     maxtab(:,2)= normGlobPowSpec_avg(maxtab(:,1));
% end
%
% pv = find(maxtab(:,2)== max(maxtab(:,2)));
% pf = round(freq(maxtab(pv,1))*100)/100;
%
% if isempty(pv) || isempty(pf)
%     maxtab =[1 1+i];
%     pv = 1;
%     pf= 1+i;
% end

% function hcb=safecolorbar(varargin) % AP code
% colorbarlocation = 'NorthOutside'; % AP code. Default = EastOutside (type help colorbar for more location options)
% vv=sscanf(version,'%i.');
% %      hcb=colorbar(varargin{:}); % Grinsted code
% if strcmpi(colorbarlocation,'North')
%     hcb=colorbar(varargin{:},'Location',colorbarlocation,'Box','on','Xcolor','k','YColor','k');
% else
%     hcb=colorbar(varargin{:},'Location',colorbarlocation,'Xcolor','k');
% end
