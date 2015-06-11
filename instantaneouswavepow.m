
function varargout = instantaneouswavepow(Wxy)
% INSTANTANEOUSWAVEPOW - Computes the total XW power across all frequencies
%                        at each instant in time.
% meanPow  = instantaneouswavepow(Wxy);
% [meanPow, maxPow,stdPow] = instantaneouswavepow(Wxy);
% Inputs:
% Wxy - matrix of crosswavelet coefficients across of size w x n, where
%       w is number of frequency points, and n is the number of time pts.
% Outputs:
% meanPow - Mean power across all frequencies at each point time point
% maxPow - Max power ... at each time point

varargout{1} = mean(abs(Wxy),1);
varargout{2} = max(abs(Wxy),[],1);
varargout{3} = std(abs(Wxy),[],1);

