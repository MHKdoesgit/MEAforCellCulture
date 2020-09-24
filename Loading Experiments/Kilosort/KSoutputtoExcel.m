

function KSoutputtoExcel(ksdir, expdir)

% load ks data
ksdat = struct2array(load([ksdir,'/Kilosortedrasters.mat']));
% get array size for selecting correct xlsx file
meaelctnum = size(ksdat.sortinginfos{1}.templates,3);

% path to xlsx folder
excelpath = 'D:/Functions/Excel files/';
slpos = strfind(expdir,filesep);

switch lower(expdir(slpos(1)+1:slpos(2)-1))
    
    case '1-color'
        if eq(meaelctnum,252)
            exlfile = 'Color_252_electrode_sortmodelrating.xlsx';
        elseif eq(meaelctnum,60)
            exlfile = 'Color_60electrode_sortmodelrating.xlsx';
        end
        
    case '2-marmoset'
        if eq(meaelctnum,252)
            exlfile = 'Marmoset_252electrode_sortmodelrating.xlsx';
        elseif eq(meaelctnum,60)
            exlfile = 'Marmoset_60electrode_sortmodelrating.xlsx';
        end
        
    case '3-rd1-opto-mglur6'
        if eq(meaelctnum,252)
            exlfile = 'OptomGluR6_252electrode_sortmodelrating.xlsx';
        elseif eq(meaelctnum,60)
            exlfile = 'OptomGluR6_60electrode_sortmodelrating.xlsx';
        end
        
    case '4-contrast adaptation'
        if eq(meaelctnum,252)
            exlfile = '252_electrode_sortmodelrating.xlsx';
        elseif eq(meaelctnum,60)
            exlfile = '60electrode_sortmodelrating.xlsx';
        end
        
    case '5-center surround'
        if eq(meaelctnum,252)
            exlfile = '252_electrode_sortmodelrating.xlsx';
        elseif eq(meaelctnum,60)
            exlfile = '60electrode_sortmodelrating.xlsx';
        end
        
    otherwise
        Error('Aint nobody found a valid experiment path, correct expdir and come back again!');
        
end
% copy and xlsx file to the experiment folder and rename it accordingly.
copyfile([excelpath,filesep,exlfile],expdir);
newxlname = strrep(expdir(slpos(2):end),filesep,'');
newxlname = strrep(newxlname,'_','');
newxlname = ['_kilosorted_',strrep(newxlname,newxlname(9:strfind(newxlname,'MEA')+2),'')];
savexlname = [expdir,filesep,[exlfile(1:end-5),newxlname,exlfile(end-4:end)]];
movefile([expdir,filesep,exlfile],savexlname);

sortinfos = cell2mat(ksdat.sortinginfos)';

% these are the selected parameters for showing in the excel file.
chnums = [sortinfos.channelnumber]' +1;     % channel numbers (+1 added to match matlab indexing)
clusnums = [sortinfos.clusternumber]';      % cluster number for each channel
clids = [sortinfos.clusterid]';             % kilosort sorting ID
clusquality = cell2mat({sortinfos.cluster_quality}');   % user-defined values for the quality of spike sorting
spknum = cell2mat({sortinfos.channelid}');
spknum = spknum(:,2);                       % number of detected spikes per cluster
tempamp = {sortinfos.templateAmplitude}';
numtemps = cellfun(@length,tempamp);        % number templates per cluster
%maxtempamp = cellfun(@max, tempamp);

sfun = @(x)(cellfun(@num2str,num2cell(x),'un',0));      % function to prepare matrix for xlsx writing

% xlsx header
headxl = [{'.mcd Order'}, {''},{''},{''},{'Cluster'}, {'Rating'},{'Sorting ID'}, {'Spike number'},{'Template number'}];

% body of xlsx file
valstoxl = [sfun(chnums), repmat({''},length(chnums),3), sfun(clusnums), sfun(clusquality(:,2)),...
    sfun(clids),sfun(spknum),sfun(numtemps)];

xlswrite (savexlname, headxl, 1,'B39' );    % writing header, check the starting index
xlswrite (savexlname, valstoxl, 1,'B40' );  % writing body text

end