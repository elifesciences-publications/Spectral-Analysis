% Created: 26-Aug-2014
% Author: AP
function OverlayLightChannel(figAxisHandle,lightChannel,timeVec,varargin)
% OverlayLightChannel(figAxisHandle,lightChannel,'r')
% OverlayLightChannel(figAxisHandle,lightChannel,'color',[0 0.5 0])
% OverlayLightChannel(figAxisHandle,...,'linewidth',2) 
% OverlayLightChannel(figAxisHandle,,'linestyle',':')
xLim = get(figAxisHandle,'xlim');
yLim = get(figAxisHandle,'yLim');

normalizedLightChannel = (lightChannel(:)/max(lightChannel(:))) - min(lightChannel(:));
baseY = yLim(1);
ampY = (yLim(2))*2.5;
adjustedLightChannel = (normalizedLightChannel*ampY) + baseY;
axes(figAxisHandle), hold on
if nargin == 4
    plot(timeVec,adjustedLightChannel, varargin{1})
elseif nargin ==5
    plot(timeVec,adjustedLightChannel, varargin{1},varargin{2});
elseif nargin ==6
    plot(timeVec,adjustedLightChannel, varargin{1},varargin{2},varargin{3})
elseif nargin == 7
    plot(timeVec,adjustedLightChannel,varargin{1},...
        varargin{2},varargin{3},varargin{4})
elseif nargin == 9
     plot(timeVec,adjustedLightChannel,varargin{1},...
        varargin{2},varargin{3},varargin{4},varargin{5},varargin{6})
else
    errordlg('Incorrect # of Input Arguments!')
end
