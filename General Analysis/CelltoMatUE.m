

function OutputMat = CelltoMatUE (inputcell, varargin)
%
% This function is an extension to MATLAB original cell2mat function but it
% works with cells with different dimentions. If the cell size is un-even
% it generate the zero in its place.
%
% ===============================Inputs====================================
%
%   inputcell : any cell.
%
%================================Output====================================
%
%   OutputMat : matrix of inputcell.
%
%
% written by Mohammad, 23.09.2013.
% added catch part on 5.06.2015
% added transpose part on 24.05.2016.
% update to faster version by removing try, catch and slower cellfun
% functions on 28.11.2017.
% added parfor option and speed option for better performance on
% on 02.03.2018.


if nargin > 1
    analyzeoption = varargin{1};
else
    analyzeoption = 2;
end

inptsize = size(inputcell);
maxLength=max(cellfun('length',inputcell));


switch analyzeoption
    case 1  % cell to fun option
        if inptsize(2) > inptsize(1)
            OutputMat = cell2mat(cellfun(@(x)cat(1,x,nan(maxLength-length(x),1)),inputcell,'un',0));
        else
            OutputMat = cell2mat(cellfun(@(x)cat(2,x',nan(1,maxLength-length(x))),inputcell,'un',0));
        end
        
    case 2  % fastest option
        
        if inptsize(2) > inptsize(1)
            OutputMat = nan(maxLength,length(inputcell));
            for ii = 1: length(inputcell)
                OutputMat(1:length(inputcell{ii}),ii) = inputcell{ii};
            end
        else
            OutputMat = nan(length(inputcell),maxLength);
            for ii = 1: length(inputcell)
                OutputMat(ii,1:length(inputcell{ii})) = inputcell{ii};
            end
        end
        
    case 3  % parfor option, it is generally slower than normal for loop because of variable generation
        % inside the for loop, only use it for bigger data
        if inptsize(2) > inptsize(1)
            OutputMat = nan(maxLength,length(inputcell));
            parfor ii = 1: length(inputcell)
                thisinput = nan(1,maxLength);
                thisinput(1:length(inputcell{ii})) = inputcell{ii};
                OutputMat(:,ii) = thisinput;
            end
            
        else
            OutputMat = nan(length(inputcell),maxLength);
            parfor ii = 1: length(inputcell)
                thisinput = nan(1,maxLength);
                thisinput(1:length(inputcell{ii})) = inputcell{ii};
                OutputMat(ii,:) = thisinput;
            end
            
        end
end

end
