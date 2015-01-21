function crestInds = findpeaks_hht2(x, thresh)
% Find peaks.
% n = findpeaks(x)

% n    = find(diff(diff(x) > 0) < 0);
% u    = find(x(n+1) > x(n));
% n(u) = n(u)+1;

[maxT,~] = peakdet(x,thresh);
crestInds = maxT(:,1);
% troughInds = minT(:,1);