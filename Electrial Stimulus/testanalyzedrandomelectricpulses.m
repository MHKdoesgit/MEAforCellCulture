


function testanalyzedrandomelectricpulses()
exp1 = load('I:\MEA\20200903_pulsematlabsript_wtcells\exp1_chamb1\h5\electrical_stimuli\02_electricPulseStimulus_30.0repeats_at0.50Hz_electrical_stimuli.mat');
dum1 = load('I:\MEA\20200903_pulsematlabsript_dummy\exp1_chamb1_dummy\h5\electrical_stimuli\02_electricPulseStimulus_30.0repeats_at0.50Hz_electrical_stimuli.mat');

ogstim  = load('I:\MEA\Electrical Stimulus\electricPulseStimulus_30.0repeats_at0.50Hz.mat'); 


%ogstim
[ca, cb, cc] = ndgrid(ogstim.pulsedurations, ogstim.pulseamplitudes, ogstim.pulsetypes);
possiblecombs = [ca(:), cb(:), cc(:)];
numcombs = size(possiblecombs,1);

stimrepeat = cell(numcombs,1);
for ii = 1:ogstim.numrepeats
    
    randstimorder =  psudoRandPerm(ogstim.seeds(ii), 1:numcombs);
    stimrepeat{ii} = possiblecombs(randstimorder,:);
    
end
stimall = cell2mat(stimrepeat);

nChan  = 60;
[el, eltvec, eldum , eldumtvec] = deal(cell(nChan,numcombs));


for ii = 1:numcombs
    
    idx = stimall(:,1)==possiblecombs(ii,1) & stimall(:,2)==possiblecombs(ii,2) & stimall(:,3)==possiblecombs(ii,3);
    for jj = 1:nChan        
        el{jj,ii} = squeeze(exp1.pulsedata(jj,idx,:));
        eltvec{jj,ii} = squeeze(exp1.pulsetimes(jj,idx,:));
        
        eldum{jj,ii} = squeeze(dum1.pulsedata(jj,idx,:));
        eldumtvec{jj,ii} = squeeze(dum1.pulsetimes(jj,idx,:));
    end
end


for n = 1:28; cla; plot(el{32,n}','r'); hold on; plot(eldum{32,n}','k'); title(n); pause; end

col = hsv(numcombs);
for n = 1:numcombs
   subplot_tight(4,7,n)
   plot(el{33,n}','color',col(n,:)); 
end

end
