

function [stimsamplerates, fs, stimsamples, Nfiles ]  = getMCDsamplingrates(mcdpath)



if nargin < 1
   mcdpath = uigetdir('D:/','Select experiment folder with mcd files'); 
end
%get mcd filenames
mcdfilenames = dir([mcdpath,filesep,'*.mcd']);
[~, reindex]=sort(str2double(regexp(({mcdfilenames(:).name}),'\d+','match','once')));
mcdfilenames={mcdfilenames(reindex).name}'; 
Nfiles=numel(mcdfilenames);
%--------------------------------------------------------------------------
%load the dll file
[dllpath,libtoload] = getMCSdllPath();
nsresult=mexprog(18, [dllpath, filesep, libtoload]);  %set dll library
%--------------------------------------------------------------------------
%get information about the recording time
stimsamples=zeros(numel(mcdfilenames),1);
for imcd=1:numel(mcdfilenames)
    mcdpathname = [mcdpath,filesep,mcdfilenames{imcd}]; %get mcd path
    [nsresult, hfile] = mexprog(1, mcdpathname); %open file
    [nsresult, mcdfileInfo] = mexprog(3, hfile); %get file info 
    if nsresult ~= 0
        error('Aint nobody find any mex dlls here, good luck with getting sampling rates!');
    end
    stimsamples(imcd)=mcdfileInfo.TimeSpan/mcdfileInfo.TimeStampResolution;
    nsresult = mexprog(14, hfile);%close file
end


%NchanTOT=mcdfileInfo.EntityCount;
stimsamples=floor(stimsamples);

stimsamplessum = cumsum(stimsamples);
stimsamplerates = [[1;stimsamplessum(1:end-1)+1],stimsamplessum];


fs = 1/mcdfileInfo.TimeStampResolution; % sampling frequency
clear mexprog; %unload DLL
%ops.stimsamples=stimsamples;

end


function [dllpath,libtoload] = getMCSdllPath()
%GETMCSDLLPATH Summary of this function goes here

dlllocation = which('load_multichannel_systems_mcd');
dllpath = fileparts(dlllocation);

switch computer()
    case 'PCWIN'; libtoload = 'nsMCDLibraryWin32.dll';
    case 'GLNX86'; libtoload = 'nsMCDLibraryLinux32.so';
    case 'PCWIN64'; libtoload = 'nsMCDLibraryWin64.dll';
    case 'GLNXA64'; libtoload = 'nsMCDLibraryLinux64.so';
    case 'MACI64'; libtoload = 'nsMCDLibraryMacIntel.dylib';
    otherwise
        disp('Your architecture is not supported'); return;
end
end


