function varargout =  hht_dy(x,dt,varargin)
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
energyThr = 0.1;
plotBool = false;
nComps = 10;
emdType = 'partial'; 
extraPks = []; % Not yet implemented

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
        end
    end
end


if isempty(freqRange)
    minF = 1/length(x);
    maxF = 1/(0.5*dt);
else
    minF = min(freqRange);
    maxF = max(freqRange);
end
fVec = sort(round((minF:dF:maxF-dF)/dF)*dF,'descend');

fKerWid = ceil(fKerWid/dF);
tKerWid = ceil(tKerWid/(dt));
tKer = gausswin(tKerWid);
fKer = gausswin(fKerWid);
gKer  = fKer(:)*tKer(:)';
gKer = gKer/sum(gKer(:));
% gKer = Standardize(gKer);


%% IMF and frequencies

if strcmpi(emdType,'full')
%     imf = emd(x,'interp','spline');
%     for jj = 1:length(imf)
%         delInds = zeros(length(imf),1);
%         c = corrcoef(imf{jj},x);
%         if c(2) <0.5
%             delInds(jj) = 1;
%         end
%     end
%     imf(find(delInds)) = [];
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
[a,h,m] = deal(d);
for k = 1:length(imf)
    b(k) = sum(imf{k}.*imf{k});
    blah = MyEMD(imf{k},1);
    %     blah = imf_all(k);
    pkInds = union(blah.mx{1},blah.mn{1});
    dPks = 2*gradient(pkInds)*dt;
    m(k,pkInds) = dPks;
    if ~isempty(pkInds)
%         temp = interp1([0; pkInds(:); max(pkInds(end),N)],[max(dPks); dPks(:); mean(dPks)],1:N, 'cubic');
        temp = interp1([0; pkInds(:); max(pkInds(end),N)],[max(dPks); dPks(:); mean(dPks)],1:N, 'spline');
    else
        temp = m(k,:);
    end
    m(k,:) = temp;
    a(k,:) = blah.maxEnv - blah.minEnv;
    th   = angle(hilbert(imf{k}));
    foo = MyEMD(gradient(th)/dt/(2*pi),1);
    h(k,:) = foo.maxEnv;
    blah = MyEMD(gradient(imf{k})/dt/(2*pi),1);
    d(k,:) = blah.maxEnv;
    %     d(k,:) = (blah.maxEnv-blah.minEnv)/2;
    d(k,:) = d(k,:)./a(k,:);
end
m = 1./m;
[~,v] = sort(-b);
energies = b/max(b);
b = 1-energies;

%% HHT
nVec = 1:N;
H = zeros(length(fVec),length(x));
D = H;
M = H;
d((d<minF) | (d>maxF))=0;
d = round(d/dF)*dF;

h((h<minF) | (h>maxF))=0;
h = round(h/dF)*dF;

m((m<minF) | (m>maxF))=0;
m = round(m/dF)*dF;

kVec = find(energies>energyThr);
for k = kVec
    A = a(k,:);
    
    blah = (d(k,:)/dF)-(min(fVec)/dF);
    blah(blah>size(D,1)) = size(D,1);
    nzInds = find(blah>0);
    foo = D*0;
    inds = sub2ind(size(D),blah(nzInds),nVec(nzInds));
    foo(inds) = A(nVec(nzInds));
    D = D + foo;
    
    blah = (h(k,:)/dF) - (min(fVec)/dF);
    blah(blah>size(H,1)) = size(H,1);
    nzInds = find(blah>0);
    foo = H*0;
    inds = sub2ind(size(H),blah(nzInds),nVec(nzInds));
    foo(inds) = A(nVec(nzInds));
    H = H + foo;
    
    blah = (m(k,:)/dF) - (min(fVec)/dF);
    blah(blah>size(M,1)) = size(M,1);
    nzInds = find(blah>0);
    foo = M*0;
    inds = sub2ind(size(M),blah(nzInds),nVec(nzInds));
    foo(inds) = A(nVec(nzInds));
    M = M + foo;
    
end
H = flipud(conv2(H,gKer,'same'));
D = flipud(conv2(D,gKer,'same'));
M = flipud(conv2(M,gKer,'same'));

% H = flipud(H);
% D = flipud(D);
% M = flipud(M);

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
