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
% 
% Wxy_nonzero_lin = Wxy_avg_files(:,:,1,cp);
% Wxy_nonzero_lin(Wxy_nonzero_lin==0)=[]; % This removes all zero elements from the matrix and vectorizes it.
% Axy = angle(Wxy_nonzero_lin);
% nPhaseBins = min([numel(Axy), 90]); % Number of bins circular for phase histograms
% [ph_dist,th] = hist(Axy(:),nPhaseBins); % Unweighted phase histogram
% mag = abs(Wxy_nonzero_lin);
% [ph_dist3,vals] = hist3([Axy(:) mag(:)],[nPhaseBins,nPhaseBins]); % 3D bivariate (phases and power) histogram
% powmat = repmat(vals{2},size(ph_dist3,1),1);
% ph_dist_wt = ph_dist3.*powmat; % Scaling the number of elements in each bin by power (power-weighting)
% ph_dist_wt = sum(ph_dist_wt,2)'; % Power-weighted phase histogram
% 
% mphase = angle(sum(Wxy_nonzero_lin)); % Mean phase is the angle of the resultant vector obtained by summing all the wavelet coefficients
% mphase(mphase<0) = mphase(mphase<0)+ 2*pi; % Addition of 2*pi to -ve values converts angle range from 0 to 360 rather than -180 to +180
% meanPhases.avg = mphase*180/pi; % Converts radians to degrees.
% meanPhases.avg = round(meanPhases.avg*100)/100;
% sphase = circ_std(Axy(:),abs(Wxy_nonzero_lin(:)));
% stdPhases.avg = sphase*180/pi;
% stdPhases.avg = round(stdPhases.avg*100)/100;

% if isempty(ph_dist)
%     theta.avg = zeros(size(ph_dist_wt));
%     phf = 1;
% else
%     ph_dist = [ph_dist(:); ph_dist(1)]; % This will close the loop in polar plot by circularizing the vector.
%     ph_dist_wt = [ ph_dist_wt(:); ph_dist_wt(1)];
%     phase_dist.avg = ph_dist./max(ph_dist);
%     phase_dist_weight.avg = ph_dist_wt./max(ph_dist_wt);
%     theta.avg = [th(:);th(1)];
%     phf = find(phase_dist_weight.avg == max(phase_dist_weight.avg));
%     if numel(phf)~=0
%         phf = phf(1);
%     else phf = 1;
%     end
% end
% peakPhase.avg = round(theta.avg(phf)*180/pi);
% if peakPhase.avg <0, peakPhase.avg = peakPhase.avg + 360; end


rtl = [max(phaseHist_wt)/2 max(phaseHist_wt)];

% fh = mmpolar(theta.avg,phase_dist_weight.avg,...
%     'b-','border',border_toggle,'ttickdelta',ttdelta,...
%     'rtickvalue',rtl,'border','on','tticklabel',ttl);
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


