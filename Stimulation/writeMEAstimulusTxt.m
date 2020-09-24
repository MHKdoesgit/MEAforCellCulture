

function writeMEAstimulusTxt(stimtowrite, stimfilelocation, stimfilename, varargin)

if ~isinteger(stimtowrite)
    error('Hey Yo, the input should be integer, convert it first then push it for writing!');
end

if size(stimtowrite,2) ~= 2
    error('Weird shity input format, the input should have two columns [time voltage] according to manual!');
end

fileID = fopen([stimfilelocation,filesep,stimfilename,'.txt'],'W');
for ii = 1:size(stimtowrite,1)  % for loop to ensure proper behavior of fprintf
    fprintf(fileID,'%i\t%i\r\n',stimtowrite(ii,:));
end
% fwrite(fileID,stimtowrite(1:10,1),'int64');
fclose(fileID);

end