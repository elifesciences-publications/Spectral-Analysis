
function varargout = instantaneouswavepow(Wxy)
% INSTANTANEOUSWAVEPOW - Computes the total XW power across all frequencies
%                        at each instant in time.
% totPow_maxNorm  = instantaneouswavepow(Wxy);
% [totPow_maxNorm, totPow, maxPow] = instantaneouswavepow(Wxy);
% Inputs:
% Wxy - matrix of crosswavelet coefficients across of size w x n, where
%       w is number of frequency points, and n is the number of time pts.
% Outputs:
% totPow - Pan-spectral total power at each point n, normalized by # of
%          wavelet scales
% totPow_maxNorm - Max-normalized version of totPow
% maxPow - 

totPow = mean(abs(Wxy),1);
totPow_maxNorm = totPow/max(totPow);


varargout{2} = totPow;
varargout{1} = totPow_maxNorm;
varargout{3} = max(abs(Wxy),[],1);