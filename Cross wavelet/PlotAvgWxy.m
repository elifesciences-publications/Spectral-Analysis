function PlotAvgWxy(Wxy,time,freq)

arrowDensity =[18 18]; % default [25 25]
ad = mean(arrowDensity);
arrowSize = 1*30.*0.03/ad;
% arrowSize = 0; % Uncommenting this line and commenting the prev one will
% cause arrows not to display
arrowHeadSize=0.9*arrowSize*180; % default 1*arrowsize*220

pkDetThr = 0.1;


%%%%%% XWT Axes
aWxy = angle(Wxy);
Wxy(abs(aWxy)< pi/2) = 0;
cData = Wxy;
cData(cData==0) = NaN;

dt = time(2)-time(1);


powerSpectrum = sum(abs(Wxy),2);
normPowerSpectrum = powerSpectrum/max(powerSpectrum);
logPowerSpectrum = log2(powerSpectrum);
normLogPowerSpectrum = logPowerSpectrum/max(logPowerSpectrum);
[maxtab,~] = peakdet(normLogPowerSpectrum,pkDetThr);
if numel(maxtab) == 0
    maxtab(:,1) = 1;
    maxtab(:,2)= 0.5;
end



figPos = get(gcf,'position');
set(gcf,'position',[figPos(1) figPos(2) figPos(3)*1.2 figPos(4)])
ax1 = axes; box off
aPos = get(ax1,'position');
aPos = [aPos(1) aPos(2) aPos(3)*(0.8) aPos(4)];
set(ax1,'position', aPos)



Yticks = 2.^(fix(log2(min(freq))):fix(log2(max(freq))));
H = imagesc(time,log2(freq),abs(log2(cData)));
clim =get(gca,'clim'); % center color limits around log2(1)=0
clim = [0 0.9]*max(clim(2));
set(gca,'clim',clim);

HCB=safecolorbar;
loc = get(HCB,'Location');
if strcmpi(loc,'North') || strcmpi(loc,'NorthOutside') || strcmpi(loc,'South')...
        || strcmpi(loc,'SouthOutside') % Since horizontal arrangement of the colorbar changes yticks to xticks
    %
    %     set(HCB,'xtick',[clim(1):2:clim(2)]); % default(-9:3:9)
    %     set(HCB,'xtick',[-9:3:9]); % default(-9:3:9)
    %     barylbls=rats(2.^(get(HCB,'xtick')'));
    %     barylbls([1 end],:)=' ';
    %     barylbls(:,all(barylbls==' ',1))=[];
    %     set(HCB,'xticklabel',barylbls);
else
    set(HCB,'ytick',-9:3:9); % default(-9:9)
    barylbls=rats(2.^(get(HCB,'ytick')'));
    barylbls([1 end],:)=' ';
    barylbls(:,all(barylbls==' ',1))=[];
    set(HCB,'yticklabel',barylbls);
end

set(gca,'YLim',log2([min(freq),max(freq)]), ...
    'YDir','normal', ...
    'YTick',log2(Yticks(:)), ...
    'YTickLabel',num2str(Yticks'), ...
    'layer','top','fontsize',14)

ylabel('Frequency (Hz)','fontsize',14)
xlabel('Time(s)')
hold on

%%%%%%% Plotting phases
phs_dt=round(length(time)/arrowDensity(1));
tidx=max(floor(phs_dt/2),1):phs_dt:length(time);
phs_dp=round(length(freq)/arrowDensity(2));
pidx=max(floor(phs_dp/2),1):phs_dp:length(freq);
phaseplot(time(tidx),log2(freq(pidx)),aWxy(pidx,tidx),arrowSize,arrowHeadSize);

%%%%%% COI
% time = time';
% tt=[time([1 1])-dt*.5;time;time([end end])+dt*.5];
% hcoi=plot(tt,log2([freq([end 1]) coi freq([1 end])]),coiLineColor,...
%             'linewidth',1.5);

 hold off


set(ax1,'color','k','ycolor','k')
xlim([time(1) time(end)]) %%%% This line is NECESSARY to ensure that x-axis is aligned with traces below
ylims1 = get(ax1,'ylim');
yticklabels_ax1 = str2num(get(ax1,'yticklabel'));
dyticklabels_ax1 = diff(yticklabels_ax1);
if dyticklabels_ax1(1)== dyticklabels_ax1, yScale = 'linear'; else yScale = 'log'; end
yl = ylabel('Frequency (Hz)','fontsize',14);
ylpos = get(yl,'pos');
aPos = get(ax1,'position');
colormap('jet')
%% XW Power Spectrum Axes
ax2 = axes; hold on, box off
aPos2 = get(ax2,'pos');
aPos2 = [aPos(1) + aPos(3) aPos(2) aPos2(3)*(0.3) aPos(4)];
set(ax2,'pos', aPos2, 'color','w','tickdir','out','fontsize',14);
xlabel([{'Normalized'}; {'Power'}])
peakFreqs = freq(maxtab(:,1));
if strcmpi(yScale,'linear'), plot(normPowerSpectrum, freq,'k','linewidth',2);
    hold on
    plot(normLogPowerSpectrum, freq,'k:','linewidth',2)
    logPeakFreqs = peakFreqs;
else
    plot(normPowerSpectrum,log2(freq),'k','linewidth',2)
    hold on
    plot(normLogPowerSpectrum,log2(freq),'k:','linewidth',2)
    logPeakFreqs = log2(peakFreqs);
end
set(ax2,'ylim',ylims1,'xlim',[0 1],'xtick',[0.5 1],'ytick',[])
ylims2 = get(ax2,'ylim');

xvals = 0.2*ones(size(peakFreqs)); % X coordinates for text displaying peak frequencies

for pf = 1:length(peakFreqs)
    txt{pf} = [num2str(round(peakFreqs(pf)*100)/100) ' Hz'];
end

text(xvals,logPeakFreqs,txt,'fontsize',14,'color','r');


ax2b = axes('position',aPos2,'xaxislocation','top',...
    'yaxislocation','right','color','none','xscale','log','ytick',[],'ycolor','w','fontsize',14);
xt = [0.25 0.5 1];
set(ax2b,'xtick',xt,'xlim',[0 1]);
leg = {'Linear';'Log'};
legend(ax2,leg)



    function hcb=safecolorbar(varargin) % AP code
       colorbarlocation = 'NorthOutside'; % AP code. Default = EastOutside (type help colorbar for more location options)
     vv=sscanf(version,'%i.');
%      hcb=colorbar(varargin{:}); % Grinsted code
if strcmpi(colorbarlocation,'North')
       hcb=colorbar(varargin{:},'Location',colorbarlocation,'Box','on','Xcolor','k','YColor','k'); 
else
    hcb=colorbar(varargin{:},'Location',colorbarlocation,'Xcolor','k'); 
end
