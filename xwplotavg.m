
function [Wxy_avg, Wxy_iso, statMat] = xwplotavg(Wxy3d, sigsxy, time, period, coi)
% [Wxy_avg,Wxy_iso,statMat] = xwplotavg(Wxy3d, sigsxy,time,period,coi)


arrowDensity =[18 18]; % default [25 25]
ad = mean(arrowDensity);
arrowSize = 1*30.*0.03/ad;
% arrowSize = 0; % Uncommenting this line and commenting the prev one will
% cause arrows not to display
arrowHeadSize=0.9*arrowSize*180; % default 1*arrowsize*220
peakDetectionThresh = 0.1; % Peak detection in the XW global power spectrum

dt = time(2)-time(1);
isoThresh = 0.5;

sigsxy = sigsxy(:);
zsize = size(Wxy3d,3);

if zsize == 0
    errordlg('Input matrix must be 3-dimensional')
end
if length(sigsxy)~= zsize
    errordlg('Dimension mismatch in inputs. Check inputs!')
end
if size(sigsxy,2)>1
    errordlg('2nd input to the function must be a vector!')
end

freq = 1./period;
coi = 1./coi;

%% Cross variance normalization and averaging.
% for chloop = 2:length(ch)
%     chStr =['ch' num2str(ch(chloop-1)) num2str(ch(chloop))];
%     Wxy3d_norm.(chStr) = Wxy3d.(chStr)
% end
Wxy3d_norm = Wxy3d;
for nn = 1:zsize
    Wxy3d_norm(:,:,nn) = Wxy3d(:,:,nn)./sigsxy(nn); % Normalizing wavelet coefficients with signal pair crossvariance
end

Wxy_avg = mean(Wxy3d_norm,3); % Average of all inputted Wxy matrices
absPow = abs(Wxy_avg(:));
absPow(absPow==0) = []; % Ignoring non-significant regions
[pCount,pVal] = hist(absPow,200);
pProb = pCount/sum(pCount);
pCdf = cumsum(pProb);
cutoff = min(find(pCdf>=isoThresh));
pCutoff = log2(pVal(cutoff));

globPowSpec_avg = sum(abs(Wxy_avg),2)/length(time);
normGlobPowSpec_avg = globPowSpec_avg/(max(globPowSpec_avg));
logGlobPowSpec_avg = log2(globPowSpec_avg);
blah = logGlobPowSpec_avg;
blah(blah==-inf)= NaN; % log2(0) = -Inf
isoThresh_lin = (2^pCutoff);
isoThresh_lin = (isoThresh_lin)/max(globPowSpec_avg);
isoThresh_log = (pCutoff-min(blah))/(max(blah)-min(blah));
blah  = (blah-min(blah))/(max(blah)-min(blah)); % Normalization by min subtraction and amp division
normlogGlobPowSpec_avg = blah;
% sNormGlobPowSpec_avg = globPowSpec_avg(:)./log2(period(:));

%% Creating a Wxy_iso: Wxy_iso(m,n) = {0, if Wxy_iso(m,n) < Isoline Threshold
Wxy_iso = Wxy_avg;
Wxy_iso(abs(Wxy_iso)<(2^pCutoff))=0;

%% Creating a Statistics Matrix
variableNames = {'Mean Freq'; 'Mean of Peak Freqs'; 'Mean Phase'; 'Mean of Peak Phases'; 'Total Pow'; 'Mean Pow'};

%% Calculating Mean Freq & Phase from Wxy_iso
Wxy_prob = abs(Wxy_avg)/sum(abs(Wxy_avg(:)));
freqMat = repmat(freq(:),1,size(Wxy_iso,2));
powWtedFreqMat = freqMat.*Wxy_prob;
meanFreq_avg = sum(powWtedFreqMat(:));
statMat{1,1} = meanFreq_avg;

[mfvec,pfvec] = instantaneouswavefreq(Wxy_iso, freq);
tvpow_prob = sum(abs(Wxy_iso));
tvpow_prob = tvpow_prob/sum(tvpow_prob);
meanFreq_avg = sum(tvpow_prob.*mfvec);
peakFreq_avg = sum(tvpow_prob.*pfvec);

statMat{1,2} = meanFreq_avg;
statMat{2,1} = peakFreq_avg;

meanPhase = angle(sum(Wxy_iso(:)))*(180/pi);
meanPhase  = round(meanPhase*100)/100; % Rounding to 2 decimal places

[meanPhaseVec,maxPhaseVec,stdPhaseVec] = instantaneouswavephase(Wxy_iso);
meanPeakPhase = sum(maxPhaseVec.*tvpow_prob);

statMat{3,1} = meanPhase;
statMat{4,1} = meanPeakPhase;

totPow = sum(abs(Wxy_iso(:)));
nonZeroPow = abs(Wxy_iso(:));
nonZeroPow(nonZeroPow==0)=[];
meanPow = mean(nonZeroPow);

statMat{5,1} = totPow;
statMat{6,1} = meanPow;

statMat = [variableNames, statMat];


%% Determining Peak Frequency from the Averaged XW Global Power Spectrum
[maxtab,mintab] = peakdet(normGlobPowSpec_avg,peakDetectionThresh);
if numel(maxtab) == 0
    maxtab(:,1)= find(normGlobPowSpec_avg == max(normGlobPowSpec_avg));
    maxtab(:,2)= normGlobPowSpec_avg(maxtab(:,1));
end

pv = find(maxtab(:,2)== max(maxtab(:,2)));
pf = round(freq(maxtab(pv,1))*100)/100;

if isempty(pv) || isempty(pf)
    maxtab =[1 1+i];
    pv = 1;
    pf= 1+i;
end


%% Plotting XW Spectrum
CData = Wxy_avg; CData(CData==0)=nan;
monSize = getMonitorSize;
fh = figure('Name', 'XW Power Averaged Across Files','color','w');
set(fh,'position',[monSize(1)+50 monSize(2)+150 monSize(3)-150 monSize(4)-300]);
%% XW Axes
ax1 = axes; box off
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
contour(time(:),log2(freq(:)),log2(abs(Wxy_avg)),[pCutoff pCutoff],'m','linewidth',2);
% ch = contour(time(:),log2(freq(:)),log2(abs(CData)),[pCutoff pCutoff]);
% set(ch,'color','k','linestyle',':','linewidth',2)



%% Y-Axis Adjustments
set(gca,'YLim',log2([min(freq),max(freq)]), ...
    'YDir','normal', ...
    'YTick',log2(Yticks(:)), ...
    'YTickLabel',num2str(Yticks'), ...
    'layer','top','fontsize',14)
ylabel('Frequency (Hz)','fontsize',14)
xlabel('Time(s)')
hold on

%% PHASE PLOTTING
aWxy=angle(CData);
phs_dt=round(length(time)/arrowDensity(1));
tidx=max(floor(phs_dt/2),1):phs_dt:length(time);
phs_dp=round(length(freq)/arrowDensity(2));
pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
phaseplot(time(tidx),log2(freq(pidx)),aWxy(pidx,tidx),arrowSize,arrowHeadSize);
time = time(:);

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

plot(normGlobPowSpec_avg, log2(freq),'k','linewidth',2), box off
hold on
plot(normlogGlobPowSpec_avg, log2(freq),'k:','linewidth',2,'parent',ax2)
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
    'layer','top','fontsize',12,'xlim',[0 1])
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


% %% COLORBAR
% HCB=safecolorbar;
% loc = get(HCB,'Location');
% if strcmpi(loc,'North') || strcmpi(loc,'NorthOutside') || strcmpi(loc,'South')...
%         || strcmpi(loc,'SouthOutside') % Since horizontal arrangement of the colorbar changes yticks to xticks
%
%     set(HCB,'xtick',[clim(1):2:clim(2)]); % default(-9:3:9)
%     barylbls=rats(2.^(get(HCB,'xtick')'));
%     barylbls([1 end],:)=' ';
%     barylbls(:,all(barylbls==' ',1))=[];
%     set(HCB,'xticklabel',barylbls);
% else
%     set(HCB,'ytick',-9:3:9); % default(-9:9)
%     barylbls=rats(2.^(get(HCB,'ytick')'));
%     barylbls([1 end],:)=' ';
%     barylbls(:,all(barylbls==' ',1))=[];
%     set(HCB,'yticklabel',barylbls);
% end


function hcb=safecolorbar(varargin) % AP code
colorbarlocation = 'NorthOutside'; % AP code. Default = EastOutside (type help colorbar for more location options)
vv=sscanf(version,'%i.');
%      hcb=colorbar(varargin{:}); % Grinsted code
if strcmpi(colorbarlocation,'North')
    hcb=colorbar(varargin{:},'Location',colorbarlocation,'Box','on','Xcolor','k','YColor','k');
else
    hcb=colorbar(varargin{:},'Location',colorbarlocation,'Xcolor','k');
end



% plot(normlog_globalPowSpec.(fStr).(chStr), freq2,'k:','linewidth',2,'parent',ax2)
%                         leg ={'Linear'; 'Log'};
%                         ax2b = axes('position',aPos2,'xaxislocation','top',...
%                             'yaxislocation','right','tickdir','out','color','none','xscale','log','ytick',[],'ycolor','w','fontsize',11);
%                         xt = [0.25 0.5 1];
%                         set(ax2b,'xtick',xt,'xlim',[0 1]);
%                         hold off