
function [out, out_all]= getExperimentPath(mainfolder, guiflag, warnflag)

if nargin < 1, mainfolder = 'D:\Ines'; guiflag = 2; warnflag = true; end
if nargin > 1, guiflag = 2; warnflag = true; end
if nargin > 2, warnflag = false; end
if isempty(mainfolder), mainfolder = 'D:\Ines'; end

expfolders10s = dir([mainfolder,filesep,'201*']);
expfolders10s = {expfolders10s.name}';
expfolders20s = dir([mainfolder,filesep,'202*']);
expfolders20s = {expfolders20s.name}';
expfolders = [expfolders10s;expfolders20s];

exptypes = 'pulsestimulation';

expdates = extractBefore(expfolders,['_',exptypes]);
expconditions = extractAfter(expfolders,[exptypes,'_']);

% Analyzed folder names
afn = {'Cell_Responses_to_Electric_Pulse_Analysis','Electric_Pulse_Burst_Analysis',...
    'Electrical_Stimuluation_Effects_Analysis','Spontaneous_Activity_Analysis'};

estimNexps = 300;
[allexppaths, analyzedpaths, rawdatapaths, expnums] = deal(cell(estimNexps,1));
afnmatfile = cell(estimNexps, numel(afn));

iter = 1;
for ii = 1: size(expfolders)
    
    exps = dir([mainfolder,filesep,expfolders{ii},filesep]);
    exps = exps(3:end);     exps = exps([exps.isdir]); % to select only the folders
    exps = {exps.name}';
     
    for jj = 1:numel(exps)
        
        thispath = [mainfolder,filesep,expfolders{ii},filesep,exps{jj},filesep];
        allexppaths{iter} = thispath;
        expnums{iter} = exps{jj};
        
        if exist([thispath,'Data Analysis'],'dir')
            % Data Analysis folder
            datpath = [thispath,'Data Analysis',filesep];
            analyzedpaths{iter} =  [thispath,'Data Analysis',filesep];
            % Raw Data folder
            if exist([datpath,'Raw Data'],'dir')
                rawdatapaths{iter} =  [datpath,'Raw Data',filesep];
            else
                rawdatapaths{iter} = NaN;
                if warnflag,  warning(['No analysis is done for ',thispath]) ; end
            end
            
            for kk = 1:numel(afn)  % get the paths to the analyzed folders
                afnpath = [datpath,afn{kk}];
                if exist(afnpath,'dir')
                    afnmfile =  dir([afnpath,filesep,'*.mat']);
                    afnmatfile{iter, kk} = [afnmfile.folder,filesep,afnmfile.name];
                else
                    afnmatfile{iter, kk} = NaN;
                    if warnflag,   warning(['There aint no analysis here for ',afn{kk}, ' at ', thispath(1:end-1)]); end
                end
            end
            
        else
            analyzedpaths{iter} = NaN;
            if warnflag,  warning(['No analysis is done for ',thispath]) ; end
        end
        iter = iter +1;
    end
end

realexppathes = ~cellfun('isempty',analyzedpaths);

allexppaths = allexppaths(realexppathes,:);
expnums = expnums(realexppathes,:);
analyzedpaths = analyzedpaths(realexppathes,:);
rawdatapaths = rawdatapaths(realexppathes,:);
afnmatfile = afnmatfile(realexppathes,:);

switch guiflag
    case 1
        promtstr = 'Select an experiment: ';
        selmode = 'single';
    case 2
        promtstr = 'Select multiple experiments: ';
        selmode = 'multiple';
end

if guiflag > 0
    guistr = extractAfter(allexppaths,[mainfolder,filesep]);
    guistr = strrep(cellfun(@(x)x(1:end-1),guistr,'un',0),filesep,' ----------> ');
    
    selectidx = listdlg('PromptString',promtstr,'Name','Time for some data analysis!',...
        'SelectionMode',selmode,'ListString',guistr,'ListSize',[700,600],'OKString','Select');
else
    selectidx = 1: numel(analyzedpaths);
end

% first output only for the selected cells
out.mainfolder = mainfolder;
out.exppaths = allexppaths(selectidx,:);
out.expnumbers = expnums(selectidx,:);
out.expdates = extractBetween(out.exppaths,[mainfolder,filesep],['_',exptypes]);
out.conditions = extractBetween(out.exppaths,[exptypes,'_'],filesep);
out.experimenttype = exptypes;
out.DataAnalysisPaths = analyzedpaths(selectidx,:);
out.RawDataPaths = rawdatapaths(selectidx,:);
out.matFilePaths = afnmatfile(selectidx,:);
out.AnalysisFolderNames = afn;

% second output for all the experiments
out_all.mainfolder = mainfolder;
out_all.exppaths = allexppaths;
out_all.expnumbers = expnums;
out_all.expdates = expdates;
out_all.conditions = expconditions;
out_all.experimenttype = exptypes;
out_all.DataAnalysisPaths = analyzedpaths;
out_all.RawDataPaths = rawdatapaths;
out_all.matFilePaths = afnmatfile;
out_all.AnalysisFolderNames = afn;

end