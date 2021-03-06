%% The McsHDF5 Matlab tools
% This and the following documents gives a short explanation of the usage
% of the McsHDF5 Matlab tools.
%

%% Purpose
% The McsHDF5 Matlab toolbox can be used to import HDF5 files generated by
% software from Multi Channel Systems MCS GmbH. Importing the HDF5 files
% via this toolbox allows easy access to all data and metadata in the HDF5
% files, as well as basic plotting functions for the streams and datasets
% in the file.

%% Prerequisites
% This Toolbox has been tested for Matlab 2009a and newer. Matlab versions
% older than R2011a might have reduced functionality, because they need to
% use the old hdf5info functions instead of the newer h5info functions. In
% particular, it might not be possible to open result files generated by
% CMOS-MEA-Tools in Matlab R2010b and older.

%% Compatible files
% This toolbox opens HDF5 files generated by the Multi Channel DataManager
% software. It alsow opens files generated by the CMOS-MEA-Control software
% version 2 and newer, as well as the CMOS-MEA-Tools software version 2 and
% newer. Please see the respective sections in this Help file browser for
% either files generated by the DataManager or the CMOS-MEA software.

%% Documentation
% In addition to the usual help functions of the Matlab command line,
% documentation for all functions and classes in this package is provided
% by the Matlab help browser:
%
%   doc McsHDF5
%
% For more information the also see:
%%
% 
% * <helpInstalling.html Installing the Toolbox>
% * <helpDataManager.html MATLAB Tools for data generated by Multi Channel DataManager>
% * <helpCmosMea.html MATLAB Tools for data generated by CMOS-MEA software>
% 
