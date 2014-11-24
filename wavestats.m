
function varargout = wavestats(Wxy,freq,time)

nFreqPts = size(Wxy,1);
nTimePts = size(Wxy,2);
nFiles = size(Wxy,3);
nChannelPairs = size(Wxy,4);

tvmf = zeros(nFiles+2,nTimePts);
tvpf = zeros(size(tvmf));
for fn = 1:size(Wxy,3)
    [tvmf(fn,:),tvpf(fn,:)] = instantaneouswavefreq(Wxy(:,:,fn),freq);
end
tvmf(fn+1,:) = mean(tvmf(1:fn,:));
tvpf(fn+1,:) = mean(tvpf(1:fn,:));

varargout{1} = tvmf;
varargout{2} = tvpf;

end
