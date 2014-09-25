function plotwave(Wxy, time, period, coi, sig95, sigmax, sigmay)
% plotwave(wave, time, period, coi, sig95, sigmax, sigmay)
% plotwave(wave, time, period, coi, sig95, sigmax, sigmay, colorbarlimits)
% Customized plot of the wavelet transform
% Adapted from xwt.m by Alsak Grinsted et al.
% Updated: 01-Mar-2011 02:40:47

%% Some variables
CData = Wxy; CData(CData==0)=nan;
arrowDensity =[18 18]; % default [25 25]
ad = mean(arrowDensity);
arrowSize = 1*30.*0.03/ad;
% arrowSize = 0; % Uncommenting this line and commenting the prev one will
% cause arrows not to display
arrowHeadSize=0.9*arrowSize*180; % default 1*arrowsize*220
ordinateType = 'log'; % 'linear' or 'log'
coiType = 'none'; % ('fill' - fills the coi with specified color;
%  'none' - does not fill the coi)
plotcoi = 'yes'; %('yes' - plots coi)
coiLineColor = 'w';
colorScheme = 'jet'; % AP code
plotMode = 'f'; %('f' for plotting Y-axis as frequency, 'p' for plotting
% ... Y-axis as period)
sig_level_for_phases = 0.5; % Added this line so that phases can only be displayed for sig regions - AP

if strcmpi (plotMode, 'f')
    freq = 1./period;
    coi = 1./coi;
elseif strcmpi(plotMode,'p')
    freq = period;
else
    error('Check plot mode: Frequency or period?')
end
time = time(:);
if nargin <7
    error('Not enough input arguments')
end
dt = time(2)-time(1);


%% Figure
if strcmpi(ordinateType,'linear')
     Yticks = fix(min(freq)):fix(max(freq));
%     H=imagesc(time,freq,log2(abs(CData/(sigmax*sigmay))));
     Wxy_variance_normalized = log2(abs(CData/(sigmax*sigmay)));
     Wxy_max_normalized = Wxy_variance_normalized/max(Wxy_variance_normalized(:));   
     H= contourf(time,freq,Wxy_max_normalized); shading flat
    colormap(colorScheme)
else
    Yticks = 2.^(fix(log2(min(freq))):fix(log2(max(freq))));
    H=imagesc(time,log2(freq),log2(abs(CData/(sigmax*sigmay))));

    colormap(colorScheme) % AP code
end
clim=get(gca,'clim'); % center color limits around log2(1)=0
% clim=[-1 1]*max(clim(2),6); %% ok
clim = [-9 9]; % default= [-9 9] - AP
clim = fix(clim);
set(gca,'clim',clim);

%% Colorbar
% HCB=safecolorbar;
% set(HCB,'ytick',-9:3:9); % default(-9:9)
% barylbls=rats(2.^(get(HCB,'ytick')'));
% barylbls([1 end],:)=' ';
% barylbls(:,all(barylbls==' ',1))=[];
% set(HCB,'yticklabel',barylbls);

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

%% Figure axes
if strcmpi(ordinateType,'linear')
    if strcmpi(plotMode,'f')
        set(gca,'YDir','normal','layer','top')
        ylabel('Frequency (Hz)')
    elseif strcmpi(plotMode,'p')
        set(gca,'YDir','reverse','layer','top')
        ylabel('Period (s)')
    end
else
    if strcmpi(plotMode,'f')
        set(gca,'YLim',log2([min(freq),max(freq)]), ...
            'YDir','normal', ...
            'YTick',log2(Yticks(:)), ...
            'YTickLabel',num2str(Yticks'), ...
            'layer','top','fontsize',14)
        ylabel('Frequency (Hz)','fontsize',14)
    elseif strcmpi(plotMode,'p')
        set(gca,'YLim',log2([min(freq),max(freq)]), ...
            'YDir','reverse', ...
            'YTick',log2(Yticks(:)), ...
            'YTickLabel',num2str(Yticks'), ...
            'layer','top','fontsize',14)
        ylabel('Period (s)')
    end
end
xlabel('Time(s)')
hold on
% aWxy=angle(CData_sig); 
% AP code. Added the next 3 lines so as to display arrows only in
% significant regions
CData_sig = CData;
CData_sig(sig95 < sig_level_for_phases)=nan;
aWxy=angle(CData_sig);

%% Adding phases to figure
if strcmpi(ordinateType,'linear')
    phs_dt=round(length(time)/arrowDensity(1));
    tidx=max(floor(phs_dt/2),1):phs_dt:length(time);
    phs_dp=round(length(freq)/arrowDensity(2));
    pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
    phaseplot(time(tidx),freq(pidx),aWxy(pidx,tidx),arrowSize,arrowHeadSize);
else
    phs_dt=round(length(time)/arrowDensity(1));
    tidx=max(floor(phs_dt/2),1):phs_dt:length(time);
    phs_dp=round(length(freq)/arrowDensity(2));
    pidx=max(floor(phs_dp/2),1):phs_dp:length(period);
    phaseplot(time(tidx),log2(freq(pidx)),aWxy(pidx,tidx),arrowSize,arrowHeadSize);
end

%% Contour plot for significant regions - Commented this part out on 1.26.2010 to simplify emf fig modification

% if strcmpi('ordinateType','linear')
%     if find(abs(Wxy)==0)
%         [c,h] = contour(time,freq,abs(Wxy),[0 0],'k');
%         set(h,'linewidth',0)
%     else
%         [c,h] = contour(time,freq,sig95,[1 1],'k');
%         set(h,'linewidth',2)
%     end
% else
%     if find(abs(Wxy)==0)
%         [c,h] = contour(time,log2(freq),abs(Wxy),[0 0],'k');
%         set(h,'linewidth',2)
%     else
%         [c,h] = contour(time,log2(freq),sig95,[1 1],'k');
%         set(h,'linewidth',2)
%     end
% end
% 
 %% Cone Of Influence
tt=[time([1 1])-dt*.5;time;time([end end])+dt*.5];

if strcmpi(plotcoi,'no')
coi = nan.*coi;
end

if strcmpi(ordinateType,'linear')
    if strcmpi(coiType,'fill')
        hcoi=fill(tt,[freq([end 1]) coi freq([1 end])],'k');
        set(hcoi,'alphadatamapping','direct','facealpha',.5)
    else
        hcoi=plot(tt,[freq([end 1]) coi freq([1 end])],coiLineColor,...
            'linewidth',1.5);
    end
else
    if strcmpi(coiType,'fill')
        hcoi=fill(tt,log2([freq([end 1]) coi freq([1 end])]),'w');
        set(hcoi,'alphadatamapping','direct','facealpha',.5)
    else
        hcoi=plot(tt,log2([freq([end 1]) coi freq([1 end])]),coiLineColor,...
            'linewidth',1.5);
    end
    hold off
end

    %% Helper functions

    % function [cout,H]=safecontourf(varargin)
    % vv=sscanf(version,'%i.');
    % if (version('-release')<14)|(vv(1)<7)
    %     [cout,H]=contourf(varargin{:});
    % else
    %     [cout,H]=contourf('v6',varargin{:});
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