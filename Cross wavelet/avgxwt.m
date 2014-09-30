% Updated: 26-Aug-2014 (Modifications to allow for specification and overlaying of light channel on XWT)
% Author: AP


function [Wxy_avg, Wxy_iso, masterVar, W_gm, S, S_prelog, Wxy_avg_channels] = avgxwt(Wxy,freq,time,coi,sigMat,isoThresh,ch);
% [Wxy_avg,Wxy_iso,masterVar,Wxy_gm, S] = avgxwt(Wxy,freq,time,coi,sigMat,isoThresh,ch)
% Size of Wxy: m by n ny f by c, where m = number of frequency values, n =
% number of time points, and f = number of files, c = number of channel
% pairs.
% Wxy_iso = Wxy in which only regions within isolines are kept
% W_gm = (Wxy(f,t)*Wyz(f,t)*Wza(f,t)).^(1/3);

arrowDensity =[18 18]; % default [25 25]
ad = mean(arrowDensity);
arrowSize = 1*30.*0.03/ad;
% arrowSize = 0; % Uncommenting this line and commenting the prev one will
% cause arrows not to display
arrowHeadSize=0.9*arrowSize*180; % default 1*arrowsize*220
peakDetectionThresh = 0.03;
global lightChannelIndex globData globTime;

% if isoThresh > 1
%     errordlg('Isothreshold value must be < 1')
%     return;
% end
%
if isempty(isoThresh)
    isoThresh = 1;
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

Wxy_norm = Wxy; % Now already normalized in xwplotmd;
% Wxy_norm = Wxy./crossVarMat; % Normalizing Wxy by the crossvariance of the underlying signal pair


%% Averaging Across Files, Channels, and Both
clear Wxy_avg_files Wxy_avg_channels Wxy_avg
Wxy_avg_files(:,:,1,:) = mean(Wxy_norm,3);
Wxy_avg_channels(:,:,:,1) = mean(Wxy_norm,4);
% for jj = 1:nChannelPairs-1
%     blah = Wxy_norm(:,:,:,jj).*conj(Wxy_norm(:,:,:,jj+1));
% end
% W_gm(:,:,:,1) = blah.^(1/nChannelPairs);
W_gm(:,:,:,1) = -(prod(-Wxy_norm,4)).^(1/nChannelPairs);
Wxy_avg = mean(mean(Wxy_norm,4),3);

S = size(nFiles,1);
S2 = size(nFiles,1);
S_prelog = size(nFiles,1);
nFiles = size(Wxy,3);
%% Matrix of Scaling Coefficients To Rectify for Frequency-Dependent Power due to Non-linear Sampling along the Frequency Dimension
[c,f] = hist(freq,fix(length(freq)/4));
cc = spline(f,c,freq);
cc = cc/max(cc);
fnmat = repmat(cc(:),1,size(Wxy,2));

%% Computing the Strength
for fnum = 1:nFiles
    blah = W_gm(:,:,fnum);
    blahb = Wxy_avg_channels(:,:,fnum);
    blahb = blahb./fnmat;
    blah = blah./fnmat; %%%% Normalizing the power values to rectify for frequency subsampling
    blah2 = log2(blah);
    blah = abs(sum(blah(:)));
    blahb = abs(sum(blahb(:)));
    %     blah2 = log2(abs(W_gm(:,:,fnum)));
   
    blah2(blah2 ==-inf) = 0;
%     S(fnum,1) = sum(blah(:));
    S(fnum,1) = sqrt(blah);
    S2(fnum,1) = sqrt(blahb);
%     S_prelog(fnum,1) = sum(blah2(:));
    S_prelog(fnum,1) = abs(sum(blah2(:)));
end
S = round(S);
S2  = round(S2);
S_prelog = round(S_prelog);
S = [S(:) S2(:) S_prelog(:)];

clear tvmf tvpf tvpow
for cp = 1:nChannelPairs
    Wxy_temp = Wxy_avg_files(:,:,1,cp);
    absPow = abs(Wxy_temp(:));
    absPow(absPow==0) = []; % Ignoring non-significant regions
    if isempty(absPow)
        absPow = 0.001;
    end
    [pCount,pVal] = hist(absPow,100);
    pProb = pCount/sum(pCount);
    pCdf = cumsum(pProb);
%     cutoff = find(pCdf>=isoThresh,1); % Commented out 25-Jul-2014
    cutoffPow = min(absPow) + isoThresh*std(absPow);
    %%%%% Preveting cutoffPow from falling outside the range of power values (absPow)
    if cutoffPow > max(absPow)
        cutoffPow = max(absPow);
    elseif cutoffPow < min(absPow);
        cutoffPow = min(absPow);
    end
    
    cutoff = find(absPow >= cutoffPow, 1); % Added 25-Jul-2014
%     if isempty(cutoff)
%         cutoff = pCdf(2);
%     end

%     pCutoff = log2(pVal(cutoff)); % Commented out 25-Jul-2014
    
     pCutoff = log2(cutoffPow);
    
    globPowSpec(:,cp) = mean(abs(Wxy_temp),2);
    amp1 = max(globPowSpec(:,cp))- min(globPowSpec(:,cp));
    normGlobPowSpec(:,cp) = globPowSpec(:,cp)/amp1;
    logGlobPowSpec(:,cp) = log2(globPowSpec(:,cp));
    blah = logGlobPowSpec(:,cp);
    blah(blah==-inf)= NaN; % log2(0) = -Inf
    amp2 = max(blah) - min(blah);
    doh = normGlobPowSpec(:,cp);
    isoThresh_lin = min(doh) + isoThresh*std(doh); % Added 25-Jul-2014
    blah  = (blah-min(blah))/amp2; % Normalization by min subtraction and amp division
    normLogGlobPowSpec(:,cp) = blah;
    doh = normLogGlobPowSpec(:,cp);
    isoThresh_log = min(doh) + isoThresh*std(doh); % Added 25-Jul-2014
    
    Wxy_temp(abs(Wxy_temp)< cutoffPow)=0; % Revised 25-Jul-2014
    Wxy_iso(:,:,1,cp) = Wxy_temp;
    
    doh = Wxy_temp./fnmat; % Added 25-Jul-2014
    Wxy_prob = abs(doh)/sum(abs(doh(:)));
    freqMat = repmat(freq(:),1,size(doh,2));
    powWtedFreqMat = freqMat.*Wxy_prob;
    %     meanFreq(cp,:) = sum(powWtedFreqMat(:));  %%% This is one way to
    %     estimate mean frequency
    meanFreq(cp,:) = circ_mean(freqMat(:),abs(doh(:)));
    stdMeanFreq(cp,:) = circ_std(freqMat(:),abs(doh(:)));
    
    [tvmf(cp,:),tvpf(cp,:)] = instantaneouswavefreq(doh, freq);
    
    
    tvpow(cp,:)= mean(abs(doh),1);
    tvPeakPow(cp,:) = max(abs(doh),[],1);
    tvpow_prob = tvpow(cp,:);
    tvpow_prob = tvpow_prob/sum(tvpow_prob(:));
    %     meanPeakFreq(cp,:) = sum(tvpow_prob.*tvpf(cp,:));
    meanPeakFreq(cp,:) = circ_mean(tvpf(cp,:)',tvPeakPow(cp,:)');
    stdPeakFreq(cp,:) = circ_std(tvpf(cp,:)',tvPeakPow(cp,:)');
    
    meanPhase(cp,:) = angle(sum(doh(:)))*(180/pi);
    meanPhase(cp,:)  = round(meanPhase(cp,:)*100)/100; % Rounding to 2 decimal places
    stdMeanPhase(cp,:) = circ_std(angle(doh(:))*(180/pi),abs(Wxy_temp(:)));
    
    [tvmph(cp,:),tvpph(cp,:),stdph(cp,:)] = instantaneouswavephase(doh);
    %     meanPeakPhase(cp,:) = sum(tvpph(cp,:).*tvpow_prob);
    meanPeakPhase(cp,:) = circ_mean(tvpph(cp,:)',tvPeakPow(cp,:)');
    stdPeakPhase(cp,:) = circ_std(tvpph(cp,:)',tvPeakPow(cp,:)');
    
%     totPow(cp,:) = sum(abs(Wxy_temp(:)));
    
    
    totPow(cp,:) = abs(sum(doh(:))); % Revised 25-Jul-2014
    nonZeroPow = abs(doh(:));
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
        maxtab =[1 1+1i];
        pv = 1;
        pf= 1+1i;
    end
    %     coif = 1./coi;
    %     plotxwt(Wxy_avg_files(:,:,1,cp),time,freq,coif);
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
    
    figName = ['XW\, Spectrum \, Averaged\, Across\, Files\, for\, Ch\,' num2str(ch(cp)) '\, vs \,' num2str(ch(cp+1))];
    figName = ['$$' figName];
    eqtxt = ' : n^{-1} \Sigma_{file = 1}^{n} |W_{ch';
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
    
    %% Overlaying Light Trace
    if lightChannelIndex ~= -1;
    aH = gca;
    OverlayLightChannel(aH,globData(:,end),globTime,'w','linewidth',2)
    end
    
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
%     dynamicfreqpowplot
%     plotphase
    
end

[masterVar.tvmf, masterVar.tvpf, masterVar.tvpow, masterVar.tvmph, masterVar.tvpph] = ...
    deal(tvmf,tvpf,tvpow,tvmph,tvpph);
[masterVar.meanFreq,masterVar.stdMeanFreq, masterVar.meanPeakFreq,masterVar.stdPeakFreq,...
    masterVar.meanPhase,masterVar.stdMeanPhase, masterVar.meanPeakPhase, masterVar.stdPeakPhase] = ...
    deal(meanFreq,stdMeanFreq, meanPeakFreq,stdPeakFreq, meanPhase,stdMeanPhase, meanPeakPhase,stdPeakPhase);
[masterVar.totPow, masterVar.meanPow] = deal(totPow, meanPow);



% function hcb=safecolorbar(varargin) % AP code
% colorbarlocation = 'NorthOutside'; % AP code. Default = EastOutside (type help colorbar for more location options)
% vv=sscanf(version,'%i.');
% %      hcb=colorbar(varargin{:}); % Grinsted code
% if strcmpi(colorbarlocation,'North')
%     hcb=colorbar(varargin{:},'Location',colorbarlocation,'Box','on','Xcolor','k','YColor','k');
% else
%     hcb=colorbar(varargin{:},'Location',colorbarlocation,'Xcolor','k');
% end
