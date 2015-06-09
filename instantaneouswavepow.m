
function varargout = instantaneouswavepow(Wxy)
% INSTANTANEOUSWAVEPOW - Computes the total XW power across all frequencies
%                        at each instant in time.
% meanPow  = instantaneouswavepow(Wxy);
% [meanPow, maxPow] = instantaneouswavepow(Wxy);
% Inputs:
% Wxy - matrix of crosswavelet coefficients across of size w x n, where
%       w is number of frequency points, and n is the number of time pts.
% Outputs:
% meanPow - Mean power across all frequencies at each point time point
% maxPow - Max power ... at each time point

<<<<<<< HEAD
totPow = mean(abs(Wxy),1);
totPow_maxNorm = totPow/max(totPow);


varargout{2} = totPow;
varargout{1} = totPow_maxNorm;
varargout{3} = max(abs(Wxy),[],1);
=======
varargout{1} = mean(abs(Wxy),1);
varargout{2} = max(abs(Wxy),[],1);
>>>>>>> 1396fd4d24fe41e99ae1e7eb2fdfdc51759d45f1
