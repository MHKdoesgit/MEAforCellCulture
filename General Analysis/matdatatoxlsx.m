

function matdatatoxlsx()


expdp = getExperimentPath('D:\Ines',2,false);
savingpath = [expdp.mainfolder,filesep,'Population Analysis',filesep,'Excel Files',filesep];
if ~exist(savingpath, 'dir'), mkdir(savingpath); end

stimtoxlx = {'Electric_Pulse_Burst_Analysis','Spontaneous_Activity_Analysis'};

stimnames = {'01_spontaneous_activity_5min';
    '02_cathodic-anodicpulse_200dur_100ppamp_50repeats';
    '03_cathodic-anodicpulse_200dur_160ppamp_50repeats';
    '04_cathodic-anodicpulse_200dur_200ppamp_50repeats';
    '05_cathodic-anodicpulse_200dur_400ppamp_50repeats';
    '06_cathodic-anodicpulse_200dur_600ppamp_50repeats';
    '07_cathodic-anodicpulse_200dur_800ppamp_50repeats';
    '08_cathodic-anodicpulse_200dur_1000ppamp_50repeats';
    '09_spontaneous_activity_5min'};

for stimid = 1:numel(stimtoxlx)
    tic
    anatype = (ismember(expdp.AnalysisFolderNames, stimtoxlx{stimid}));
    
    if strcmpi(stimtoxlx{stimid}, 'Spontaneous_Activity_Analysis')
        sn = stimnames([1,9]);
    else
        sn = stimnames(2:8);
    end
    
    mf = expdp.matFilePaths(:,anatype);
    expwithdata = (~cellfun('isempty',mf));
    expwithoutdata = find(~expwithdata);
    if any(~expwithdata)
        warning('The following selected expriments had no analyzed data folder in them. Run the analysis first!')
        for ii = 1:sum(~expwithdata)
            warning(expdp.exppaths{expwithoutdata(ii)});
            warning(repmat('-',1,100));
        end
    end
    mf = mf(expwithdata);
    numexps = size(mf,1);
    
    dat = cell(numexps,1);
    for ii = 1: numexps
        try
            dat{ii} = struct2array( load(mf{ii}));
        catch
            dat{ii} = load(mf{ii});
        end
    end
    dat = reshape(cell2mat(dat),[],numexps);
    
    numstims = size(dat,1);
    
    %fn = fieldnames(dat(1,1).indices);
    %idxtowrite = cell(numstims,numel(fn));
    for ii = 1 : numstims
        
        idx = [dat(ii,1:end).indices];
        
        fn = fieldnames(idx);
        %fn = fn(~ismember(fn,{'pearsoncoeff','explainedvar','Rsq_allcells','pearsoncoeff_allcells','explainedvar_allcells'}));
        for jj = 1:numel(fn)
            popidx.(fn{jj}) =  CelltoMatUE({idx.(fn{jj})});
        end
        popidx.numbursts = [dat(ii,1:end).numbursts];
        popidx.peakthreshold = [dat(ii,1:end).peakthreshold];
        %popidx.spkamps = CelltoMatUE({dat(ii,1:end).meanamplitudes});
        fn = fieldnames(popidx);
        
        idxtowrite = cell(1,numel(fn));
        
        for jj = 1: numel(fn)
            idxvalues = popidx.(fn{jj});
            switch lower(fn{jj})
                case 'burstpermin', lgtxt = 'Bursts/minute';
                case 'firingrate',  lgtxt = 'Firing rate (Hz)';
                case 'ibi',  lgtxt = 'Inter-burst-intervals (sec)';
                case 'intraburstfrq',  lgtxt = 'Intera-burst spiking frequency (Hz)';
                case 'spkinburstpercent',  lgtxt = 'Spikes in burst (%)';
                case 'burstduration',  lgtxt = 'Burst duration (ms)';
                case 'numbursts',  lgtxt = 'Total detected bursts';
                case 'peakthreshold', lgtxt = 'PSTH dectection threshold';
                case 'rsq',  lgtxt = 'Reliability R^2';
                case 'totalspikes',   lgtxt = 'Total number of spikes (all cells)';
                case 'totalburstyspikes',   lgtxt = 'Total number of spikes (bursty cells)';
                case 'allchanfiringrate',   lgtxt = 'Firing rate (all cells)';
                case 'allmeanisi',   lgtxt = 'Average inter-spike-intervals';
                case 'meanisibursty',   lgtxt = 'Average inter-spike-intervals (bursty cells)';
                case 'allchansmeanamplitudes',   lgtxt = 'Average spike amplitudes (all cells)';
                case 'meanamplitudes',   lgtxt = 'Average spike amplitudes (bursty  cells)';
                otherwise, lgtxt = fn{jj};
            end
            
            gaprow = repmat({''},1,numexps);
            f2w = [expdp.expdates(expwithdata)'; expdp.conditions(expwithdata)'; expdp.expnumbers(expwithdata)';...
                gaprow; repmat({lgtxt},1,numexps);gaprow; num2cell(idxvalues)];
            idxtowrite{1,jj} =  f2w;
        end
        
        xlsname  = [savingpath,filesep, sn{ii},'.xlsx'];
        
        for qq = 1:numel(fn)
            warning( 'off', 'MATLAB:xlswrite:AddSheet' ) ;
            xlsdat = idxtowrite{1,qq};
            xlrange = [xlRC2A1(6,4),':',xlRC2A1(size(xlsdat,1)+6,size(xlsdat,2)+4)];
            writecell(xlsdat, xlsname,'Sheet',qq,'Range',xlrange);
            
        end
        warning( 'on', 'MATLAB:xlswrite:AddSheet' ) ;
    end
    toc
end





%
% fn = fieldnames(matdat);
% nsheets = numel(fn);
% filename = [dp,filesep, excelname,'.xlsx'];
%
% for ii = 1: nsheets-5
%
%     n = fn{ii};
%     d = matdat.(fn{ii});
%     if isrow(d), d=d'; end
%
%     xlrange = [xlRC2A1(6,4),':',xlRC2A1(size(d,1)+6,size(d,2)+4)];
%
%     %dd = {'CellNumber',[];1:size(d,1),d};
%     if iscell(d)
%
%         %writecell(dd,filename,'Sheet',ii,'Range',xlrange);
%     elseif isstruct(d)
%
%     else
%         writematrix('cell number',filename,'Sheet',ii,'Range',xlRC2A1(4,2));
%         writematrix(n,filename,'Sheet',ii,'Range',xlRC2A1(4,10));
%         writematrix(transpose(1:size(d,1)),filename,'Sheet',ii,'Range',[xlRC2A1(6,2),':',xlRC2A1(size(d,1)+6,2)]);
%         writematrix(d,filename,'Sheet',ii,'Range',xlrange);
%     end
% end

end



function xlCELL = xlRC2A1(xlROW,xlCOL)
% Returns the column characters of Excel given a certain column number
% Input ROW : row number
%       COL : number of column
% Output CHAR : Character combination in Excel
if xlCOL <= 26                        % [A..Z]
    xlCHAR = char(mod(xlCOL-1,26)+1+64);
    xlCELL = [xlCHAR num2str(xlROW)];
elseif xlCOL <= 702                   % [AA..ZZ]
    xlCOL = xlCOL-26;
    xlCHAR1 = char(floor((xlCOL-1)/26)+1+64);
    xlCHAR0 = char(mod(xlCOL-1,26)+1+64);
    xlCHAR = [xlCHAR1 xlCHAR0];
    xlCELL = [xlCHAR num2str(xlROW)];
elseif xlCOL <= 16384                 % [AAA..XFD]
    xlCOL = xlCOL-702;
    CHAR2 = char(floor((xlCOL-1)/676)+1+64);
    xlCOL=xlCOL-(floor((xlCOL-1)/676))*676;
    xlCHAR1 = char(floor((xlCOL-1)/26)+1+64);
    xlCHAR0 = char(mod(xlCOL-1,26)+1+64);
    xlCHAR = [CHAR2 xlCHAR1 xlCHAR0];
    xlCELL = [xlCHAR num2str(xlROW)];
else
    disp('Column does not exist in Excel!');
end
end