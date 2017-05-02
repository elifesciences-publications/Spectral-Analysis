function varargout =  hht_freq(x,dt,varargin)
% Get and plot imds and HHT.
% H = hht(x,Ts)
% [M, H, imf, freq, energies] = hht(x,dt,'freqRange',freqRange,'dF',dF,
%   'plotBool',plotBool,'tKerWid',tKerWid,'fKerWid',fKerWid,
%   'energyThr',energyThr,'nComps',nComps);
% :: Syntax
%    The array x is the input signal and dt is the sampling period.
%    Example on use: [x,Fs] = wavread('Hum.wav');
%                    plot_hht(x(1:6000),1/Fs);
% Requisites : emd
% Inputs:
% x - Signal to get HHT for.
% dt - Sampling interval.
% freqRange - A 2 element vector that defines the frequency range of
%   interest.
% dF - Frequency resolution. Default = 0.5 Hz.
% plotBool - Boolean that determines whether or not the imfs from emd
%   should be plotted
% tKerWid - Width of temporal kernel for smoothing HHT [Default: 1].
% fKerWid - Width of frequency kernel for smoothing HHT [Default: 5].
% energyThr - Energy threshold as a fraction of the IMF with maximal
%   energy. IMFs with energy less than this threshols are ommitted from
%   consideration when generating the HHT
% nComps - Number of IMFs to extract
% emdType = 'full' or 'partial' [default]. If 'partial', emd stops after
%   first iteration for a component without fulfilling the criterion 
%   wherein the number of extrema must match the number of zerocrossings
% extraPks - Extra peaks that can be specified to improve interpolation; 2
%   cell array. First cell is for maxPks and 2nd is for minPks
%
% Modified by Avinash Pujala, JRC/HHMI, 2017

freqRange = [];
dF = 0.5;
tKerWid = 1;
fKerWid = 5;
energyThr = 0;
plotBool = false;
nComps = 20;
emdType = 'full'; 
extraPks = []; % Not yet implemented
freqVec = [];

for jj = 1:numel(varargin)
    if ischar(varargin{jj})
        switch lower(varargin{jj})
            case 'freqrange'
                freqRange = varargin{jj+1};
            case 'df'
                dF = varargin{jj+1};
            case 'plotbool'
                plotBool = varargin{jj+1};
            case 'tkerwid'
                tKerWid = varargin{jj+1};
            case 'fkerwid'
                fKerWid = varargin{jj+1};
            case 'energythr'
                energyThr = varargin{jj+1};
            case 'ncomps'
                nComps = varargin{jj+1};
            case 'emdtype'
                emdType = varargin{jj+1};
            case 'extrapks'
                extraPks = varargin{jj+1};
            case 'freqvec'
                freqVec = varargin{jj+1};
        end
    end
end

if ~isempty(freqVec)
    fVec = freqVec;
    [minF, maxF] = deal(min(freqVec),max(freqVec));
else
    if isempty(freqRange)
        minF = 1/length(x);
        maxF = 1/(0.5*dt);
    else
        minF = min(freqRange);
        maxF = max(freqRange);
    end
    fVec = sort(round((minF:dF:maxF-dF)/dF)*dF,'descend');
end
fKerWid = max(1,ceil(fKerWid/dF));
tKerWid = max(1,ceil(tKerWid/dt));
tKer = gausswin(tKerWid);
fKer = gausswin(fKerWid);
gKer  = fKer(:)*tKer(:)';
gKer = gKer/sum(gKer(:));

%% IMF and frequencies
if strcmpi(emdType,'full')
    imf = emd(x,'nComps',nComps);
else
    imf_all = MyEMD(x,nComps);
    imf = cell(length(imf_all),1);
    for jj = 1:length(imf_all)
        imf{jj} = imf_all(jj).comp;
    end
end

N = length(x);
d = zeros(length(imf),N);
[a,h,m,a_env] = deal(d);
b = zeros(length(imf),1);
for k = 1:length(imf)
    b(k) = sum(imf{k}.*imf{k});
    blah = MyEMD(imf{k},1);   
    pkInds = union(blah.mx{1},blah.mn{1});
    dPks = 2*gradient(pkInds)*dt;
    m(k,pkInds) = dPks;
    if ~isempty(pkInds)
        temp = interp1([0; pkInds(:); max(pkInds(end),N)],[max(dPks); dPks(:); mean(dPks)],1:N, 'spline');
    else
        temp = m(k,:);
    end
    m(k,:) = temp;
    a_env(k,:) = blah.maxEnv - blah.minEnv;
    a(k,:) = abs(hilbert(imf{k}));
    th   = angle(hilbert(imf{k}));
    h(k,:) = gradient(th)/dt/(2*pi);    
end
m = 1./m;
energies = 1 -b/max(b);

%% HHT
nVec = 1:N;
H = zeros(length(fVec),length(x));
M = H;
h((h<minF) | (h>maxF))=0;
m((m<minF) | (m>maxF))=0;
if isempty(freqVec)
    h = ceil(h/dF)*dF;
    m = ceil(m/dF)*dF;
end

% kVec = find(energies>energyThr);
kVec = 1:length(imf);
for k = kVec(:)'
    A = a(k,:);
    [f_hist,binInds] = histc(h(k,:),fVec(:));
    zerInds = binInds == 0;
    binInds(zerInds) = [];
    nVec_temp = nVec;
    nVec_temp(zerInds)=[];
    inds = sub2ind(size(H),binInds,nVec_temp);    
    foo = H*0;
    foo(inds) = A(nVec_temp);
    H = H + foo;
    
    A_env = a_env(k,:);    
    [f_hist,binInds] = histc(m(k,:),fVec(:));
    zerInds = binInds == 0;
    binInds(zerInds) = [];
    nVec_temp = nVec;
    nVec_temp(zerInds)=[];
    inds = sub2ind(size(M),binInds,nVec_temp);    
    foo = M*0;
    foo(inds) = A_env(nVec_temp);
    M = M + foo;    
end
% H = flipud(conv2(H,gKer,'same'));
% M = flipud(conv2(M,gKer,'same'));

H = conv2(H,gKer,'same');
M = conv2(M,gKer,'same');
M(M<0) = 0;

%% IMF plots
if plotBool
    K = length(imf);
    c = linspace(0,(N-1)*dt,N);
    nSub = 5;
    for k1 = 0:nSub:K-1
        figure('Name','IMFs');
        for k2 = 1:min(nSub,K-k1)
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
            if k2 ~= min(nSub,K-k1)
                set(gca,'xtick',[])
            end
        end
        xlabel('Time');
        axes(ax)
        title(['IMFs: ' num2str(k1+1) '-' num2str(k1+k2)])
    end
end
%% Outputs
varargout{1} = M;
varargout{2}= H;
varargout{3} = imf;
varargout{4} = fVec;
varargout{5} = energies;

end
