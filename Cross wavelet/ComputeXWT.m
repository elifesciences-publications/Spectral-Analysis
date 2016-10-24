function varargout = ComputeXWT(x,y,time,varargin)

% My custom written XWT based on xwt.m by Grinsted et al
% [Wxy,freq,coi,sig95] = ComputeXWT(x,y,time);
% [Wxy,...] = ComputeXWT(x,y,dt);
% [Wxy,...] = ComputeXWT(x,y,time,'freqRange',freqRange,'dj',dj,'stringency',stringency,...
%       'phaseType',phaseType,'sigmaXY',sigmaXY,'pad',pad,'freqScale',freqScale)

%% Fixed and variable parameters
fourier_factor      = 1.0330; % Conversion factor for changing wavelet scales to periods (true for wavenumber = 6)
pad                 = 1; % Zero padding of signals (Pad = 1: Pad with zeroes; Pad = 0: Don't zero pad)
mother              = 'Morlet'; %%%('Morlet', 'Paul','DOG') - For now use only Morlet, other wavelets will give erroneous results
dj                  = 1/12; %(default = 1/12 for log freq scale and 1 for linear freq scale);
stringency          = 1;
phaseType           = 'all'; %('all','alt','synch')
freqRange           = []; % Automatically computes if []
sigmaXY             = []; % Automatically computes if []
freqScale           = 'log'; %('log' or 'lin')
noiseType           = 'white'; % ('white' or 'red')


if nargin < 3
    error('At least 3 inputs required. If using for single timeseries, input the same timeseries twice!')
end

for jj = 1:numel(varargin)
    if ischar(varargin{jj})
        switch lower(varargin{jj})
            case 'pad'
                pad = varargin{1};
            case 'mother'
                mother = varargin{jj+1};
            case 'dj'
                dj = varargin{jj+1};
            case 'phaseType'
                phaseType = varargin{jj+1};
            case 'stringency'
                stringency = varargin{jj+1};
            case 'sigmaxy'
                sigmaXY = varargin{jj+1};
            case 'freqrange'
                freqRange = varargin{jj+1};
            case 'freqscale'
                freqScale = varargin{jj+1};
            case 'noisetype'
                noiseType = varargin{jj+1};
        end
    end
end

x = x(:);
y = y(:);
time = time(:);
if numel(time)==1 && numel(x)>1
    dt = time;
    time = (0:length(x)-1)*(1/dt);
else
    dt = time(2)-time(1);
end

if isempty(freqRange)
    freqRange(1) = 1/length(x);
    freqRange(2) = 1/dt;
else
    if freqRange(2) > floor(1/dt)
        errordlg('High frequency value above Nyquist limit, please respecify')
    end
end

if isempty(sigmaXY)
    sigmaXY = std(x)*std(y);
end

scaleRange = 1./(freqRange*fourier_factor); % Scale range corresponding to frequency range.
strcmpi(freqScale,'lin')
S0 = min(scaleRange);
maxScale = max(scaleRange);

x = [time(:) x(:)];
y = [time(:) y(:)];
[Wxy,period,~,coi,sig95]= xwt(x,y,'dj',dj,'S0',S0, 'ms', maxScale, 'Mother', mother,'pad',pad,'noiseType',noiseType,'freqScale',freqScale);
freq =1./period;
Wxy = Wxy/sigmaXY;


%% Removing regions outside COI
ftmat = repmat(freq(:), 1, length(time));
coimat = repmat(1./coi(:)',length(freq), 1);
% Wxy(ftmat < coimat) = 0; % Removing regions outside of COI

%% Removing insignificant points
% Wxy(sig95 < stringency) = 0;
Wxy(abs(Wxy) < stringency)= 0;

%% Phase-based filtering
a  = angle(Wxy);
if strcmpi(phaseType,'alt')
    Wxy(abs(a)<=(0.5*pi))=0; % Keeps only those matrix elements with angles >= (0.75*pi) = 135 deg
elseif strcmpi(phaseType,'synch')
    Wxy(abs(a)>pi/2)=0; % Keeps only those matrix elements with angles <= pi/2
elseif strcmpi(phaseType,'all')
    % Do nothing
else
    errordlg('Please input proper phase type')
end

varargout{1} = Wxy;
varargout{2} = freq;
varargout{3} = coi;
varargout{4} = sig95;


end
