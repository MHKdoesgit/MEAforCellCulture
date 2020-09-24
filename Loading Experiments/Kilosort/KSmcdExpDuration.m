

function varargout = KSmcdExpDuration(mcdpath, savingpath)


if nargin < 1
    mcdpath = uigetdir('D:/','Select experiment folder with mcd files');
end
if nargin < 2,    savingpath = [];  end
if isnumeric(savingpath), savingpath = mcdpath; end
if ~isempty(savingpath), savetxt = true; else, savetxt = false; end

[stimsamplerates, fs, stimsamples ]  = getMCDsamplingrates(mcdpath);


stimtimes = stimsamplerates ./ fs;
stimdur = char(duration(0, 0, stimsamples ./  fs));
stpoint = num2str(stimtimes(:,1),6);
endpoint = num2str(stimtimes(:,2),6);

% get the names for meta data
mcdlists = dir([mcdpath,filesep,'*.mcd']);

if isempty(mcdlists)
    error('Hey yo!, there aint no recoreded MC data in this folder! good luck with analysis');
end

mtdat.mcdfilenames = {mcdlists.name}';
mtdat.mcdfilesize = [mcdlists.bytes]'/(2^10^3);
mtdat.totalexpsize = sum(mtdat.mcdfilesize);
mtdat.exptime = {mcdlists.date}';
[str,dt, strtym] = deal(cell(size(mcdlists,1),1));


for jj = 1:size(mcdlists,1)
    
    expsize = num2str(mtdat.mcdfilesize(jj),'%.3g');
    if length(expsize)>4, expsize = expsize(1:4); end
    if length(expsize)<4, expsize = [expsize,repmat(' ',1,4-length(expsize))]; end %#ok
    
    str{jj} = [num2str(jj,'%02d'),':',repmat(' ',1,5),mtdat.mcdfilenames{jj}(1:15),' ... ',...
        mtdat.mcdfilenames{jj}(end-3:end), repmat(' ',1,5),'size: ',expsize,...
        ' GB', repmat(' ',1,5),'recorded at: ', mtdat.exptime{jj}];
    gp = strfind(mtdat.exptime{jj},':');
    dt{jj} = mtdat.exptime{jj}(1:gp(1)-4);
    
    strtym{jj} = [num2str(jj,'%02d'),':',repmat(' ',1,5),mtdat.mcdfilenames{jj}(1:15),' ... ',...
        mtdat.mcdfilenames{jj}(end-3:end), repmat(' ',1,5),'started: ',stpoint(jj,:),repmat(' ',1,5),...
        'end: ',endpoint(jj,:),' sec', repmat(' ',1,5),'duration: ',stimdur(jj,:)];
end
mtdat.expdate = cell2mat(unique(dt));
mtdat.label = str;
mtdat.timelabel = strtym;
mtdat.samplerate = fs;
mtdat.stimsamplerates = stimtimes;
mtdat.stimduration = stimsamples ./  fs;

% printing shit on screen
l = 100;
fprintf('\n');
fprintf([repmat('=',1,l),'\n']);
fprintf([repmat(' ',1,l/2-9),' Experiment data ',repmat(' ',1,l/2-9),'\n']);
fprintf([repmat('=',1,l),'\n']);
fprintf('\n');

disp(str);
disp([repmat('-',1,40),'> Total size: ', repmat(' ',1,5), num2str(mtdat.totalexpsize,'%.3g'),' GB']);
disp([repmat('-',1,40),'> Experiment data: ',repmat(' ',1,5), mtdat.expdate]);

fprintf('\n');
fprintf([repmat('=',1,l),'\n']);
fprintf([repmat(' ',1,l/2-10),' Experiment timing ',repmat(' ',1,l/2-10),'\n']);
fprintf([repmat('=',1,l),'\n']);
fprintf('\n');

disp(strtym);
disp([repmat('-',1,40),'> Experiment length: ',repmat(' ',1,5), char(duration(0, 0, sum(stimsamples ./fs)))]);

% prepare the text string for the saved text file
namelen = cellfun(@length,mtdat.mcdfilenames);
l = 150;
if savetxt
    fid=fopen([savingpath,filesep,'Experiment timing info.txt'],'w');
    
    fprintf(fid,'\n');
    fprintf(fid,[repmat('=',1,l),'\n']);
    fprintf(fid, [repmat(' ',1,l/2-10),' Experiment timing ',repmat(' ',1,l/2-10),'\n']);
    fprintf(fid,[repmat('=',1,l),'\n']);
    fprintf(fid,'\n');
    for jj = 1:size(mcdlists,1)
        expname = [mtdat.mcdfilenames{jj},repmat(' ',1,max(namelen)-namelen(jj)+5)];
        fprintf(fid,[num2str(jj,'%02d'),':',repmat(' ',1,5),expname,repmat(' ',1,5),...
            'started: ',stpoint(jj,:),repmat(' ',1,5),...
            'end: ',endpoint(jj,:),' sec', repmat(' ',1,5),'duration: ',stimdur(jj,:),'\n']);
    end
    
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,[repmat('=',1,l),'\n']);
    fprintf(fid, [repmat(' ',1,l/2-9),' Experiment data ',repmat(' ',1,l/2-9),'\n']);
    fprintf(fid,[repmat('=',1,l),'\n']);
    fprintf(fid,'\n');
    for jj = 1:size(mcdlists,1)
        expname = [mtdat.mcdfilenames{jj},repmat(' ',1,max(namelen)-namelen(jj)+5)];
        expsize = num2str(mtdat.mcdfilesize(jj),'%.3g');
        if length(expsize)>4, expsize = expsize(1:4); end
        if length(expsize)<4, expsize = [expsize,repmat(' ',1,4-length(expsize))]; end %#ok
        fprintf(fid,[num2str(jj,'%02d'),':',repmat(' ',1,5),expname,repmat(' ',1,5),'size: ',...
            expsize,' GB', repmat(' ',1,5),'recorded at: ', mtdat.exptime{jj},'\n']);
    end
    
    fclose(fid);
end

if nargout > 1
    varargout{1} = mtdat;
end

end