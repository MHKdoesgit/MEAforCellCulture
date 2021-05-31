

function electStimpara =  getElectStimParameters(h5path, h5stimname, varargin)


stimpath = [fileparts(h5path),filesep,'stimuli',filesep,extractBefore(h5stimname,'.h5'),'.txt'];

if not(exist(stimpath,'file'))
    warning(['Da FAQ, there is no stimulus with the name: ',extractBefore(h5stimname,'.h5'),...
        ' in the stimuli folder. check your stimuli names!']);
    electStimpara = {[]}; 
    return; 
end

fid = fopen(stimpath,'r');
stimclass = fgetl(fid);
stimparaperstim.stimulus = stimclass;
while ischar(stimclass)
    tline = fgetl(fid);
    if tline == -1, break; end
    tparam = textscan(tline,'%s%s','Delimiter','=','TreatAsEmpty',' ');
    if isnan(str2double(tparam{2}))
    
    stimparaperstim.(strtrim( tparam{1}{1})) = cell2mat(tparam{2});
    else
        stimparaperstim.(strtrim( tparam{1}{1})) = str2double (tparam{2});
    end
end
fclose(fid);
if ~ ismember(stimparaperstim.unit,{'s','sec','ms','microsec','ns'})
    stimparaperstim.unit = 'microsec';
end
if strcmpi(stimparaperstim.stimulus, 'cathodic_anodic_pulse')
    stimparaperstim.interpulseinterval = 4800;
end

electStimpara = stimparaperstim;


end



