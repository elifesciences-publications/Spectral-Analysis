function PlotPhaseHist(varargin)
% PlotPhaseHist Plots rose diagrams for phases obtained from XWT in
%   ComputeAndPlotXWTs. Displays both weighted(by XW power) and
%   unweighted rose diagrams
% PlotPhaseHist(phaseHist, phaseAngles)
% PlotPhaseHist(phaseHist, phaseHist_wt, phaseAngles)
% PlotPhaseHist(..., meanPhase)
% Inputs:
% phaseHist - Phase histogram
% phaseHist_wt - Power-weighted phase histogram
% phaseAngles - Phase angles over which histogram is computed
% meanPhase - Mean Phase in degrees
% 
% Avinash Pujala, 2015

%% Check input
if nargin == 2
    phaseHist = varargin{1};
    phaseAngles = varargin{2};
elseif nargin == 3
    phaseHist = varargin{1};
    phaseHist_wt = varargin{2};
    phaseAngles = varargin{3};
elseif nargin == 4
    phaseHist = varargin{1};
    phaseHist_wt = varargin{2};
    phaseAngles = varargin{3};
    meanPhase = varargin{4};
else
    errordlg('Incorrect # of inputs')
end

%% Display parameters
border_toggle = 'on';
ttdelta = 90; % Spacing in degrees between angular ticks
tl = [0:ttdelta:360];
ttl= [];
for jj = 1:length(tl)
    ttl{jj} = [num2str(tl(jj)) '\circ'];
    ttl = ttl(:);
end
rtl = [max(phaseHist_wt)/2 max(phaseHist_wt)];
fh = mmpolar(phaseAngles,phaseHist_wt,...
    'b-','border',border_toggle,'ttickdelta',ttdelta,...
    'rtickvalue',rtl,'border','on','tticklabel',ttl);
set(gcf,'color','w'), set(fh,'linewidth',2)
hold on
fh = mmpolar(phaseAngles,phaseHist,...
    'color',[0 0.5 0],'linestyle',':','border',border_toggle,'ttickdelta',ttdelta,...
    'rtickvalue',rtl,'border','on','tticklabel',ttl);
set(gcf,'color','w'), set(fh,'linewidth',2,'markersize',15)

if nargin == 4
meanPhase = meanPhase * pi/180
pol = polar(meanPhase*ones(1,2),[0 1]);
set(pol,'color','r','linewidth',2,'linestyle','--','linewidth',2)
legend('Weighted', 'Unweighted','Mean Phase')
else
   legend('Weighted', 'Unweighted')
end
hold off


