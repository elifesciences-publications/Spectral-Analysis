function varargout = PlotWxy(varargin)
%PlotWxy  Function that accepts a matrix of crosswavelet coefficients and
%          plots its image
%
% PlotWxy(Wxy)
% PlotWxy(...,time)
% PlotWxy(...,freq)
% PlotWxy(...,[x,y])
% PlotWxy(...,clim)
% fh = PlotWxy(...)
%
% Inputs:
% Wxy - Crosswavelet coeff matrix
% time - Time vector corresponding to the col entries in Wxy
% freq - Freq vector corresponding to row entries in Wxy
% [x,y] - Signals that produced Wxy; in included, these will be plot in separate
%   ...axes under near Wxy.
% clim - Color limits; A 2-element array that specifies the lower and upper
%   ... limits of the color scale

Wxy = varargin{1};
if nargin < 2
    time = 1:size(Wxy,1);
    freq = 1:size(Wxy,1);
    sigs = [];
elseif nargin < 3
    time = varargin{2};
    freq = 1:size(Wxy,1);
    sigs = [];
elseif nargin < 4
    time = varargin{2};
    freq = varargin{3};
    sigs = [];
elseif nargin ==4
    time = varargin{2};
    freq = varargin{3};
    sigs = varargin{4};
elseif nargin == 5
    if strcmpi(varargin{4},'clim')
        time = varargin{2};
        freq= varargin{3};
        sigs =[];
        cLim = varargin{5};
    else
        error('Color limits not specified properly')
        return;
    end
elseif nargin ==6
    time = varargin{2};
    freq= varargin{3};
    sigs = varargin{4};
    cLim = varargin{6};
else
    error('Incorrect # of inputs')
    return
end


arrowDensity =[18 18]; % default [25 25]
ad = mean(arrowDensity);
arrowSize = 1*30.*0.03/ad;
% arrowSize = 0; % Uncommenting this line and commenting the prev one will
% cause arrows not to display
arrowHeadSize = 0.9*arrowSize*180; % default 1*arrowsize*220

pkDetThr = 0.1;


%%%%%% XWT Axes
cData = Wxy;
cData(cData==0) = NaN;
dt = time(2)-time(1);

powerSpectrum = sum(abs(Wxy),2)/size(Wxy,2);
normPowerSpectrum = powerSpectrum/max(powerSpectrum);
logPowerSpectrum = log2(powerSpectrum);
logPowerSpectrum(isinf(logPowerSpectrum))= nan;
normLogPowerSpectrum = logPowerSpectrum/max(logPowerSpectrum);
[maxtab,~] = peakdet(normLogPowerSpectrum,pkDetThr);
if isempty(maxtab)
    maxtab(1:2) = nan;
end

figPos = get(gcf,'position');
set(gcf,'position',[figPos(1) figPos(2) figPos(3)*1.2 figPos(4)])
hold on
% Colorbar
HCB = safecolorbar;
loc = get(HCB,'Location');
if strcmpi(loc,'North') || strcmpi(loc,'NorthOutside') || strcmpi(loc,'South')...
        || strcmpi(loc,'SouthOutside') % Since horizontal arrangement of the colorbar changes yticks to xticks
else
    set(HCB,'ytick',-9:3:9); % default(-9:9)
    barylbls=rats(2.^(get(HCB,'ytick')'));
    barylbls([1 end],:)=' ';
    barylbls(:,all(barylbls==' ',1))=[];
    set(HCB,'yticklabel',barylbls);
end

if ((nargin == 3) || isempty(sigs))
    aPos1 = [0.8 1 0 0]; % Wavelet axes
    aPos2 = [0.2 1 0.8 0]; % Global power spectrum axes - linear plot
    aPos2b  = aPos2; % Global power spectrum axes - log plot
    axHandles = CreateSubaxes(gcf, aPos1, aPos2, aPos2b);
elseif nargin > 3
    aPos1 = [0.8 0.8 0 0.2]; % Wavelet axes
    aPos2 = [0.2 0.8 0.8 0.2]; % Global power spectrum axes - linear plot
    aPos2b  = aPos2; % Global power spectrum axes - log plot
    aPos3 = [0.8 0.2 0 0]; % Timeseries axes
    axHandles = CreateSubaxes(gcf,aPos1, aPos2, aPos2b, aPos3);
    ax3 = axHandles(4);
end

ax1 = axHandles(1);
ax2 = axHandles(2);
ax2b = axHandles(3);

Yticks = 2.^(fix(log2(min(freq))):fix(log2(max(freq))));
axes(ax1), hold on
H = imagesc(time,log2(freq),abs(log2(cData)));
if (nargin)
end
if ~exist('cLim')
    cLim =get(gca,'clim'); % center color limits around log2(1)=0
    cLim = [0 0.9]*max(cLim(2));
end
safecolorbar;

if  isempty(sigs)
      set(gca,'clim',cLim,'YLim',log2([min(freq),max(freq)]), ...
        'YDir','normal', ...
        'YTick',log2(Yticks(:)), ...
        'YTickLabel',num2str(Yticks'), ...
        'layer','top')
else
      set(gca,'clim',cLim,'YLim',log2([min(freq),max(freq)]), ...
        'YDir','normal', ...
        'YTick',log2(Yticks(:)), ...
        'YTickLabel',num2str(Yticks'), ...
        'layer','top','xtick',[])
end

ylabel('Frequency (Hz)')
set(ax1,'color','k','ycolor','k')
if nargin == 3
    xlabel('Time')
end
xlim([time(1) time(end)]) %%%% This line is NECESSARY to ensure that x-axis is aligned with traces below
ylims1 = get(ax1,'ylim');
yticklabels_ax1 = str2num(get(ax1,'yticklabel'));
dyticklabels_ax1 = diff(yticklabels_ax1);
if dyticklabels_ax1(1)== dyticklabels_ax1, yScale = 'linear'; else yScale = 'log'; end
yl = ylabel('Frequency (Hz)');
ylpos = get(yl,'pos');
% aPos = get(ax1,'position');
colormap('jet')

% Plotting phases
% aWxy = angle(Wxy);
aWxy = angle(cData);
phs_dt=round(length(time)/arrowDensity(1));
tidx=max(floor(phs_dt/2),1):phs_dt:length(time);
phs_dp=round(length(freq)/arrowDensity(2));
pidx=max(floor(phs_dp/2),1):phs_dp:length(freq);
phaseplot(time(tidx),log2(freq(pidx)),aWxy(pidx,tidx),arrowSize,arrowHeadSize);

% COI
% time = time';
% tt=[time([1 1])-dt*.5;time;time([end end])+dt*.5];
% hcoi=plot(tt,log2([freq([end 1]) coi freq([1 end])]),coiLineColor,...
%             'linewidth',1.5);

hold off


%% XW Power Spectrum Axes
axes(ax2); hold on
set(ax2,'color','w');
xlabel('Power')
peakFreqs = freq(maxtab(:,1));
if strcmpi(yScale,'linear'), plot(powerSpectrum, freq,'k','linewidth',2);
    hold on
    plot(logPowerSpectrum, freq,'k:','linewidth',2)
    logPeakFreqs = peakFreqs;
else
    %     plot(normPowerSpectrum,log2(freq),'k','linewidth',2)
    plot(normPowerSpectrum,log2(freq),'k','linewidth',2)
    hold on
    plot(normLogPowerSpectrum,log2(freq),'k:','linewidth',2)
    logPeakFreqs = log2(peakFreqs);
end
% set(ax2,'ylim',ylims1,'xlim',[0 1],'xtick',[0.5 1],'ytick',[])
set(ax2,'ylim',ylims1,'xlim',[-inf inf],'ytick',[])
ylims2 = get(ax2,'ylim');

xvals = 0.2*ones(size(peakFreqs)); % X coordinates for text displaying peak frequencies

for pf = 1:length(peakFreqs)
    txt{pf} = [num2str(round(peakFreqs(pf)*100)/100) ' Hz'];
end

text(xvals,logPeakFreqs,txt,'color','r');
axes(ax2b)
set(ax2b,'xaxislocation','top','yaxislocation','right',...
    'color','none','xscale','log','ytick',[],'ycolor','w');
% xt = [0.25 0.5 1];
% set(ax2b,'xtick',xt,'xlim',[-inf inf]);
set(ax2b,'xlim',[min(logPowerSpectrum) max(logPowerSpectrum)]);
leg = {'Linear';'Log'};
legend(ax2,leg)

if ((nargin > 3) && (~isempty(sigs)))
    axes(ax3)
    plot(time,sigs(:,1),'b')
    if size(sigs,2) ==2
        hold on
        plot(time,sigs(:,2),'r')
    end
    xlim([time(1) time(end)])
    xlabel('Time(s)')
end


function hcb = safecolorbar(varargin)
colorbarlocation = 'North'; % Default = EastOutside (type help colorbar for more location options)
vv=sscanf(version,'%i.');
if strcmpi(colorbarlocation,'North')
    hcb=colorbar(varargin{:},'Location',colorbarlocation,'Box','on','Xcolor','k','YColor','k');
else
    hcb=colorbar(varargin{:},'Location',colorbarlocation,'Xcolor','k');
end
