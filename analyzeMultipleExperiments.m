
function analyzeMultipleExperiments()

t=tic;
expdp = getExperimentPath('D:\Ines',2,false);
% finished analysis
% pathlist = {'D:\Ines\20210408_pulsestimulation_no_treatment_ctr_EVs\exp41_chamb41_no_treatment_ctr';
%     'D:\Ines\20210408_pulsestimulation_no_treatment_ctr_EVs\exp42_chamb42_no_treatment_ctr';
%     'D:\Ines\20201124_pulsestimulation_EctoExoCtr\exp15_chamb15_ectoctr';
%     'D:\Ines\20201124_pulsestimulation_EctoExoCtr\exp16_chamb16_ectoctr';
%     'D:\Ines\20201124_pulsestimulation_EctoExoCtr\exp17_chamb17_ectoctr';
%     'D:\Ines\20201124_pulsestimulation_EctoExoCtr\exp18_chamb18_exoCtr';
%     'D:\Ines\20201124_pulsestimulation_EctoExoCtr\exp19_chamb19_exoCtr';
%     'D:\Ines\20201124_pulsestimulation_EctoExoCtr\exp20_chamb20_exoCtr';
%     'D:\Ines\20201201_pulsestimulation_EctoExocells\exp25_chamb25_ectoCtr'};

% new analysis
% pathlist = {'D:\Ines\20201201_pulsestimulation_EctoExocells\exp25_chamb25_ectoCtr'};
pathlist = expdp.exppaths;

addpath(genpath ('C:\MEA\Matlab_functions\MEAforCellCulture'));

for ii = 1:size(pathlist,1)
    
    disp(repmat('=',1,100));
    disp(pathlist{ii});
    disp(repmat('=',1,100));
    InesAanalyzeExperiments(pathlist{ii});
    disp('Analysis is finito');
end

seconds2human(toc(t))


end