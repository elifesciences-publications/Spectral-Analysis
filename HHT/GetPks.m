function pks = GetPks(x)
% GetPks - Given a signal x, returns all peaks in the signal
% pks = findpeaks(x)

pks    = find(diff(diff(x) > 0) < 0);
u    = find(x(pks+1) > x(pks));
pks(u) = pks(u)+1;
end
