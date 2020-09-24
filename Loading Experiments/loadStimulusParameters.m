

function [stimpara, stimlist] = loadStimulusParameters(datapath)

stimlistpath = [datapath,filesep,'stimuli'];

if not(exist(stimlistpath,'dir')) 
    error('there aint no stimuli list here, you should have a stimuli folder with list of stimuli!'); 
end

fprintf('loading stimuli parameters....');
stimlist = dir([stimlistpath,filesep,'*.txt']);
stimlist = {stimlist.name}';
stimpara = cell(size(stimlist,1),1);

for ii = 1: size(stimlist,1)
    
    fid = fopen([stimlistpath,filesep,stimlist{ii}]);
    stimclass = fgetl(fid);
    stimparaperstim.stimulus = stimclass;
    while ischar(stimclass)
        tline = fgetl(fid);
        if tline == -1, break; end
        tparam = textscan(tline,'%s%f','Delimiter','=','TreatAsEmpty',' ');
        
        stimparaperstim.(strtrim( tparam{1}{1})) = tparam{2};
    end
    fclose(fid);
    stimpara{ii} = stimparaperstim;
end

%stimpara = cell2mat(stimpara);
stimlist = extractBefore(stimlist,'.txt');
fprintf('done\n');

end