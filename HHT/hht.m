function varargout =  hht(x,dt,varargin)
% Get and plot imds and HHT.
% H = hht(x,Ts)
% [H, emd, energies] = hht(x,Ts,'freqRange',freqRange,'dF',dF,'plotBool',plotBool);
% :: Syntax
%    The array x is the input signal and dt is the sampling period.
%    Example on use: [x,Fs] = wavread('Hum.wav');
%                    plot_hht(x(1:6000),1/Fs);
% Func : emd

% Get HHT.
%
% Modified by Avinash Pujala, JRC/HHMI, 2017

freqRange = [];
dF = 0.5;

win = gausswin(10);
gker  = win(:)*win(:)';
gKer = gker/100;


for jj = 1:numel(varargin)
    if ischar(varargin{jj})
        switch lower(varargin{jj})
            case 'freqrange'
                freqRange = varargin{jj+1};
            case 'df'
                dF = varargin{jj+1};
            case 'plotbool'
                plotBool = varargin{jj+1};
        end
    end
end

if isempty(freqRange)
    minF = 1/length(x);
    maxF = 1/(0.5*dt);
else
    minF = freqRange(1);
    maxF = freqRange(end);
end
fVec = round((minF:dF:maxF-dF)/dF)*dF;
 
%% IMF and frequencies
imf = emd(x);
N = length(x);
d = zeros(length(imf),N-1);
for k = 1:length(imf)
    b(k) = sum(imf{k}.*imf{k});
    th   = angle(hilbert(imf{k}));
    d(k,:) = diff(th)/dt/(2*pi);
end
[~,v] = sort(-b);
energies = b/max(b);
b = 1-energies;

%% HHT

nVec = 1:N;
H = zeros(length(fVec),length(x));
% d = abs(d);
d((d<minF) | (d>maxF))=0;
d = round(d/dF)*dF;
kVec = find(energies>0.1);
for k = kVec
    blah = d(k,:)/dF;
    blah(blah>size(H,1)) = size(H,1);
    nzInds = find(blah~=0);
    foo = H*0;
    inds = sub2ind(size(H),blah(nzInds),nVec(nzInds));
    foo(inds) = abs(imf{k}(nVec(nzInds)));
    H = H + foo;
end
H = conv2(H,gker,'same');

% Set time-frequency plots.
% N = length(x);
% c = linspace(0,(N-2)*dt,N-1);
% for k = v(1:2)
%    figure, plot(c,d{k},'k.','Color',b([k k k]),'MarkerSize',3);
%    set(gca,'FontSize',8,'XLim',[0 c(end)],'YLim',[0 1/2/dt]);
%    xlabel('Time'),
%    ylabel('Frequency');
% end

%% IMF plots
if plotBool
    M = length(imf);
    c = linspace(0,(N-1)*dt,N);
    nSub = 5;
    for k1 = 0:nSub:M-1
        figure('Name','IMFs');
        for k2 = 1:min(nSub,M-k1)
            if k2 ==1
                ax =  subplot(4,1,k2);
            else
                subplot(nSub,1,k2);
            end
            plot(c,imf{k1+k2});
            box off
            set(gca,'FontSize',8,'XLim',[0 c(end)], 'tickdir','out');
            ylabel(['imf ' num2str(k1+k2)])
            ylim([min(x),max(x)])
            if k2 ~= min(nSub,M-k1)
                set(gca,'xtick',[])
            end
        end
        xlabel('Time');
        axes(ax)
        title(['IMFs: ' num2str(k1+1) '-' num2str(k1+k2)])
    end
end
%% Outputs

varargout{1} = H;
varargout{2} = imf;
varargout{3} = fVec;
varargout{4} = energies;

end
