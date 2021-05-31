

function out = InesAanalyzeExperiments(datapath)


totaltime =tic;
if (nargin < 1),    datapath = uigetdir();  end

addpath(genpath('D:\ztest\Ines\functions'));
rawdatapath = [datapath,filesep,'Data Analysis',filesep,'Raw Data'];

if ~exist(rawdatapath,'dir')
    loadExperimentData(datapath,1,4);
end

out.spontaneous = Analyze_Spontaneous_Activity(datapath);
[out.burst, out.pulse, out.elecstim] = Analyze_Multi_Electrical_Pulses(datapath);

sound(struct2array(load('gong.mat','y')));
disp(seconds2human (toc(totaltime)))

end