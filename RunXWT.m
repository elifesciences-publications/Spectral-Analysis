% RunWavelet - Allows for selection of multiple data files to load for susequent wavelet analysis
% Runs mfload, newpreprocessdata, and xwplotmd (in that order). Many other
% custom-written scripts are also executed in the process.

% addpath(genpath('C:\Users\pujalaa\Documents\MATLAB\General'));
% addpath(genpath('C:\Users\pujalaa\Documents\MATLAB\Spectral-Analysis'));

%% Load Multiple Files at Once and Create a Single Large Data Structure
% Use loadfiles instead of mfload to singly load data from files into
% MATLAB workspace. The only requirement of mfload is that all the files
% must be on the same path (i.e. in the same directory or folder)

LoadFiles

%% Preprocesses/Conditions the Loaded Data Prior to Wavelet Transformation
% (1) Aligns signals with respect to an event (such as onset of electrical 
% stimulation) for averaging or comparison across mutliple trials (in 
% this case, a single trial = a single file.
% (2) Filters the signals and subsamples as specified for subsequent WT.
% (3) Removes electrial or light-related artifacts if specified

Preprocess


%% Runs Wavelet Routines on Conditioned Time Series Signals

ComputeAndPlotXWTs


