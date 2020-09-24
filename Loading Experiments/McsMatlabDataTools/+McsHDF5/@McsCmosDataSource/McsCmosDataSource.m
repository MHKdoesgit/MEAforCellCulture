classdef McsCmosDataSource < handle
% Stores the contents of a generic data source for CMOS-MEA files. 
%
% This should only be used if a data source of unknown type is used. It
% reads all datasets in the data source, so it should not be used if the
% datasets are very large!
%
% (c) 2017 by Multi Channel Systems MCS GmbH

    properties (SetAccess = protected)
        Label           % (string) The name of the data source
        
        % Datasets - (cell array) The data sets in the data source. For
        % each dataset, its Name, Attributes and Content are read from the
        % file
        Datasets = {};  
        
        % Info - (struct) Information about the stream
        % The fields depend on the stream type and hold information about
        % the datasource
        Info
    end
    
    methods
        function str = McsCmosDataSource(filename, strStruct, varargin)
        % Contains the contents data source of unknown type
        %
        % function str = McsCmosDataSource(filename, strStruct, varargin)
        %
        % This function reads a data source group from a CMOS-MEA file, if
        % the type of the data source does not fit to any of the known
        % types. It stores the attributes of the group in its Info field
        % and reads all the datasets in to the Datasets array. If the data
        % source contains any nested groups, these will not be read.
        %
        %   filename    -   (string) Name of the HDF5 file
        %   strStruct   -   The data source subtree of the structure 
        %                   generated by the h5info command
            if exist('h5info')
                mode = 'h5';
            else
                mode = 'hdf5';
            end
            
            str.Label = McsHDF5.McsH5Helper.GetFromAttributes(strStruct, 'ID.Instance', mode);
            str.Info = McsHDF5.McsH5Helper.ReadInfoFromAttributes(strStruct, mode);
            
            for di = 1:length(strStruct.Datasets)
                dataset = [];
                dataset.Name = McsHDF5.McsH5Helper.GetFromAttributes(strStruct.Datasets(di), 'ID.Instance', mode);
                info = [];
                for ai = 1:length(strStruct.Datasets(di).Attributes)
                    [name, value] = McsHDF5.McsH5Helper.AttributeNameValueForStruct(strStruct.Datasets(di).Attributes(ai), mode);
                    info.(name) = value;
                end
                dataset.Info = info;
                if strcmp(mode, 'h5')
                    dataset.Content = h5read(filename, [strStruct.Name '/' strStruct.Datasets(di).Name]);
                elseif strcmp(mode, 'hdf5')
                    dataset.Content = hdf5read(filename, [strStruct.Name '/' strStruct.Datasets(di).Name]);
                end
                
                str.Datasets = [Datasets {dataset}];
            end
        end
    end
end