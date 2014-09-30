function varargout = plotxwt(varargin)
% PLOTXWT  Displays the spectrum for specified XW coefficients using
% custom settings.
% plotxwt(Wxy_norm,time,freq,coi)
% fh = plotxwt(Wxy_norm,time,freq,coi)
% fh = figure handle
% Wxy_norm is the normalized version of Wxy
% Wxy_norm = (Wx.*conj(Wy))/(std(x)*std(y)); (see XWT.m by Grinsted for an explanation)
% Created by AP on 24-May-2014

%% Specifying Parameters
arrowDensity =[18 18]; % Default [25 25] - determines spacing of phase arrows on XWT
ad = mean(arrowDensity);
arrowSize = 1*30.*0.03/ad;
% arrowSize = 0; % Uncommenting this line and commenting the prev one will
% cause arrows not to display
arrowHeadSize=0.9*arrowSize*180; % default 1*arrowsize*220

peakDetectionThresh = 0.25; % Threshold for Peak Detection in the Global
% Power Spectrum (see help peakdet.m for details)
yAxis = 'log';

Wxy = varargin{1};
if nargin < 4
    errordlg('PLOTXWT must have at least 4 inputs!')
    return
elseif size(Wxy,3) > 1
    errordlg('The first input to PLOTXWT must have only 2 dimensions!')
    return
end

time = varargin{2};
freq = varargin{3};
coi = varargin{4};
coi = 1./coi; % This is because the coi is aligned with period, whereas
              %  XWT will be plotted with frequency on the y-axis.
period = 1./freq;
dt = time(2)-time(1);

%% Normalized Global Power Spectrum - Linear Scale
globPowSpec = sum(abs(Wxy),2)/size(Wxy,2);
amp1 = max(globPowSpec)- min(globPowSpec);
% globPowSpec_norm = globPowSpec/amp1;
globPowSpec_norm = globPowSpec;

%% Normalized Global Power Spectrum - Log Scale
globPowSpec_log = log2(globPowSpec);
blah = globPowSpec_log;
blah(blah==-inf)= NaN; % Since, log2(0) = -Inf
minLogPow = min(blah);
maxLogPow = max(blah);
blah = blah-minLogPow;
amp2= max(blah);
blah  = blah * (amp1/amp2); 
globPowSpec_log_norm = blah;
if amp1 > 0
    peakDetectionThresh = peakDetectionThresh*amp1;
end

%% Finding & Displaying Peaks in the Global Power Spectrum
[maxtab,mintab] = peakdet(globPowSpec_norm,peakDetectionThresh);
if (isempty(Wxy)| all(Wxy==0))
    maxtab =[1 1];
    pv = 1;
    pf= 1;
elseif isempty(maxtab) % In the event the peak detection threshold is too stringent
    maxtab(:,1)= find(globPowSpec_norm == max(globPowSpec_norm));
    maxtab(:,2)= globPowSpec_norm(maxtab(:,1));
end
pv = find(maxtab(:,2)== max(maxtab(:,2)));
pf = round(freq(maxtab(pv,1))*100)/100;


%% Plotting Wxy
CData = Wxy; CData(CData==0)=nan;
monSize = getMonitorSize;
figName = ['XW Spectrum'];
fh = figure('Name', figName,'color','w');
set(fh,'position',[monSize(1)+50 monSize(2)+150 monSize(3)-150 monSize(4)-300]);

%% XW Axes
ax1 = axes; box off
aPos = get(ax1,'position');
aPos = [aPos(1)*0.8 aPos(2) aPos(3)*(0.95) aPos(4)];
set(ax1,'position', aPos)

if strcmpi(yAxis,'lin')
Yticks = [fix(min(freq)):fix(max(freq))];
H=imagesc(time,freq,log2(abs(CData)));
else
Yticks = 2.^(fix(log2(min(freq))):fix(log2(max(freq))));
H=imagesc(time,log2(freq),log2(abs(CData)));    
end

colormap(jet)
clim=get(gca,'clim'); % center color limits around log2(1)=0

clim = [-9 9]; % default= [-9 9]
clim = fix(clim);
set(gca,'clim',clim);


%% Y-Axis Adjustments
if strcmpi(yAxis,'lin')
     set(gca,'YLim',[min(freq),max(freq)], ...
    'YDir','normal', ...
    'YTick',Yticks(:), ...
    'YTickLabel',num2str(Yticks'), ...
    'layer','top','fontsize',14)
else
    set(gca,'YLim',log2([min(freq),max(freq)]), ...
    'YDir','normal', ...
    'YTick',log2(Yticks(:)), ...
    'YTickLabel',num2str(Yticks'), ...
    'layer','top','fontsize',14)
end
ylabel('$$ Frequency\ (Hz) $$','interpreter','latex','fontsize',14)
xlabel('$$ Time\ (s) $$','interpreter','latex','fontsize', 14)

hold on

%% Colorbar - AP code
HCB=safecolorbar;
loc = get(HCB,'Location');
if strcmpi(loc,'North') || strcmpi(loc,'NorthOutside') || strcmpi(loc,'South')...
        || strcmpi(loc,'SouthOutside') % Since horizontal arrangement of the colorbar changes yticks to xticks
%     
    set(HCB,'xtick',[clim(1):2:clim(2)]); % default(-9:3:9)
%     set(HCB,'xtick',[-9:3:9]); % default(-9:3:9)
    barylbls=rats(2.^(get(HCB,'xtick')'));
    barylbls([1 end],:)=' ';
    barylbls(:,all(barylbls==' ',1))=[];
    set(HCB,'xticklabel',barylbls);
else
    set(HCB,'ytick',-9:3:9); % default(-9:9)
    barylbls=rats(2.^(get(HCB,'ytick')'));
    barylbls([1 end],:)=' ';
    barylbls(:,all(barylbls==' ',1))=[];
    set(HCB,'yticklabel',barylbls);
end

hold on
%% Phase Plotting
aWxy=angle(CData);
phs_dt=round(length(time)/arrowDensity(1));
tidx=max(floor(phs_dt/2),1):phs_dt:length(time);
phs_dp=round(length(freq)/arrowDensity(2));
pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
if all(isnan(aWxy(:)))
    blah =[];
else
    if strcmpi(yAxis,'lin')
        phaseplot(time(tidx),freq(pidx),aWxy(pidx,tidx),arrowSize,arrowHeadSize);
    else
      phaseplot(time(tidx),log2(freq(pidx)),aWxy(pidx,tidx),arrowSize,arrowHeadSize);
    end

end
time = time(:);

%% COI
tt=[time([1 1])-dt*.5;time;time([end end])+dt*.5];
if strcmpi(yAxis,'lin')
    hcoi=plot(tt,[freq([end 1]) coi freq([1 end])],'w','linewidth',1.5);
else
    hcoi=plot(tt,log2([freq([end 1]) coi freq([1 end])]),'w','linewidth',1.5);
end

figTitle = ['$\log_2|W_{xy}|$'];
title(figTitle,'interpreter','latex','fontsize', 14)
hold off

%% Plotting Global Power Spectrum

aPos = get(ax1,'position');

ax2 = axes;
aPos2 = get(ax2,'position');
aPos2 = [aPos(1) + aPos(3) aPos(2) aPos2(3)*(0.15) aPos(4)];
set(ax2,'position', aPos2, 'color','none','tickdir','out','fontsize',11);

% if imag(maxtab(:,2))
%     maxtab(:,2)= 1+i;
%     peakFreqs = 1+i;
% else
%     peakFreqs = freq(maxtab(:,1));
% end
peakFreqs = freq(maxtab(:,1));
plot(globPowSpec_norm, log2(freq(:)),'k','linewidth',2), box off
xl = ['$|W_{xy}|_{norm} $'];
xlabel(xl,'interpreter','latex')
hold on
plot(globPowSpec_log_norm, log2(freq(:)),'k:','linewidth',2,'parent',ax2)
maxLinPow = max(globPowSpec_norm);
% maxLogPow = max(globPowSpec_log_norm);
% minLogPow = min(globPowSpec_log_norm);
leg ={'Linear'; 'Log'};
ax2b = axes('position',aPos2,'xaxislocation','top',...
    'yaxislocation','right','tickdir','out','color','none','ytick',[],'ycolor','w','fontsize',11);
if maxLinPow ==0
    amp1 = 1;
    maxLinPow = 1;
    maxLogPow = 1;
    minLogPow = 0;
    amp2 = 1;
end
rf = ceil(1-log10(maxLinPow)); %%% No rounding for numbers with >= 1 digits.
rf = 10^rf;

xlim_lin = ([0 1.05])*maxLinPow;
xtick_lin = (round(([0.5 1])*maxLinPow*rf))/rf;
set(ax2,'YLim',log2([min(freq),max(freq)]), ...
    'YDir','normal', ...
    'YTick',[], ...
    'layer','top','fontsize',12,'xlim',xlim_lin,'xtick', xtick_lin)

% xlim_log = [0 log2(xlim_lin(2))];
xlim_log = [minLogPow maxLogPow];
amp2 = maxLogPow - minLogPow;
xtick_log = minLogPow + ([0.25 0.5 1])*amp2;

set(ax2b,'xlim', xlim_log,'xtick', xtick_log);
xticklabel_log = 2.^xtick_log(:); 
rf = ceil(1-log10(maxLogPow)); %%% No rounding for numbers with >= 1 digits.
rf = 10^rf;
xticklabel_log = (round(xticklabel_log*rf))/rf;
xticklabel_log = rats(xticklabel_log);
% xticklabel_log = num2str(2.^str2num(get(ax2b,'xticklabel')));
set(ax2b,'xticklabel',xticklabel_log,'ytick',[],'xdir','normal');


box off

xvals = 0.35*amp1*ones(size(peakFreqs));
clear txt
for pf = 1:length(peakFreqs)
    txt{pf} = [num2str(round(peakFreqs(pf)*100)/100) ' Hz'];
    % txt = {round(peakFreqs*100)/100};
end
logPeakFreqs = log2(peakFreqs);
pfLabelPos = 0.5*peakFreqs;
text(xvals,logPeakFreqs,txt,'fontsize',11,'color','r','parent',ax2);

legend(ax2,leg,'fontsize',10) % This line needs to be here to legend can be moved by hand after fig is generated
hold off

 function hcb=safecolorbar(varargin) % AP code
       colorbarlocation = 'NorthOutside'; % AP code. Default = EastOutside (type help colorbar for more location options)
     vv=sscanf(version,'%i.');
%      hcb=colorbar(varargin{:}); % Grinsted code
if strcmpi(colorbarlocation,'North')
       hcb=colorbar(varargin{:},'Location',colorbarlocation,'Box','on','Xcolor','k','YColor','k'); 
else
    hcb=colorbar(varargin{:},'Location',colorbarlocation,'Xcolor','k'); 
end
 end
end