

function InesDatatofolder(datapath)

%datapath= 'P:\Ines\MEA_2019_05_21_INES';
tic;
msrdnames = dir([datapath, '/*.msrd']);         msrdnames = {msrdnames.name}';

for ii = 1:numel(msrdnames)
    
    expname = msrdnames{ii}(1:end-length('.msrd'));
    expfolder = [datapath,'/',expname];
    
    if ~exist(expfolder,'dir'), mkdir(expfolder); end
    
    expfiles = dir([datapath,'/*',expname,'*']);
    expfiles = {expfiles.name}';
    for jj= 1:numel(expfiles)
        thisexpfile = [datapath,'/',expfiles{jj}];
        if isfolder(thisexpfile), continue; end
        ismoved = movefile(thisexpfile, [expfolder,'/',expfiles{jj}]);
        
        if ~ismoved
            error('Some shit files or shitty format  could not be moved!');
        end
    end
    
end
toc;

end