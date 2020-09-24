

function pulselecstim = pulsestimulusforCellCulture(pulsedurations, stimfreq, pulseamplitudes, pulsetypes, numrepeats, savingpath)


if nargin < 6, savingpath = []; end

% pulsedurations = [50, 100, 150, 250];
% pulseamplitudes = [160,200,400,600,800,1000];
% pulsetypes = [1,2,3,4];
% numrepeats = 5;

[ca, cb, cc] = ndgrid(pulsedurations, pulseamplitudes, pulsetypes);
possiblecombs = [ca(:), cb(:), cc(:)];

pulsefreq = (1/stimfreq)*1e6;       % amplifier units are in ?S
iter = pulsefreq;                   % start one stimfrq late
numcombs = size(possiblecombs,1);

%seedset = randi([-100000,-100],nrepeats,1);
seedset = linspace(-10000,-1000,numrepeats);
%seed = -1000;

stimvec = cell(1,numrepeats);
stimrandoder = zeros(numrepeats, numcombs);

for ii = 1:numrepeats
    
    stimrandoder(ii,:) = psudoRandPerm(seedset(ii), 1:numcombs);
    randcombs = possiblecombs(stimrandoder(ii,:),:);
    
    
    stimsteps = cell(1,numcombs);
    for jj = 1:numcombs
        
        stimpoint = generateElectricPulses(randcombs(jj,1),randcombs(jj,2),0,randcombs(jj,3));
        %tdelays = iter - max(stimpoint(1,:));
        tdelays = iter;
        
        stimsteps{jj} = [[tdelays;0],[tdelays+1;0]+stimpoint];
        %stimsteps{jj} = [[tdelays;0]+stimpoint];
        
        iter = iter+pulsefreq;
        
    end
    if ii== 1
        stimvec{ii} = [[0;0],cell2mat(stimsteps)];
    else
        %stimvec{ii} = [[iter;0],cell2mat(stimsteps)];
         stimvec{ii} = cell2mat(stimsteps);
    end
end

stimvecall = cell2mat(stimvec); % int64 for the amplifier
stimvecall(2,:) = stimvecall(2,:) * 1e3; % convert to milivolt
stimvecall = int64 (stimvecall);

recdur = double(max(stimvecall(1,:))/60e6);
fprintf('total required recording time:\t\t%0.2f minutes (%.2f hours).\n',recdur, recdur/60);

pulselecstim.stimvec = stimvecall;
pulselecstim.stimorder = stimrandoder;
pulselecstim.stimcombinations = possiblecombs;
pulselecstim.seeds = seedset;
pulselecstim.pulsedurations = pulsedurations;
pulselecstim.pulseamplitudes = pulseamplitudes;
pulselecstim.pulsetypes = pulsetypes;
pulselecstim.numrepeats = numrepeats;


if ~isempty(savingpath)
    wrfilename = sprintf('electricPulseStimulus_%.1frepeats_at%.2fHz',numrepeats,stimfreq);
    writeMEAstimulusTxt(stimvecall', savingpath, wrfilename);
    % save everthing to mat file to be sure!
    save([savingpath,filesep,wrfilename,'.mat'],'-v7.3','-struct','pulselecstim');
end



end


