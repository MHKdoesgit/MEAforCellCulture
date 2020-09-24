classdef McsCmosRecording < handle
% Stores a single recording for HDF5 files produced by the CMOS-MEA
% software
%
% The different data sources present in the file are represented by the
% different properties of this object. Typically, a file generated by
% CMOS-MEA-Control will only contain an Acquisition group, while files
% generated by CMOS-MEA-Tools can have other data sources as well.
%
% (c) 2017 by Multi Channel Systems MCS GmbH

    properties (SetAccess = private)
        Acquisition             % McsCmosAcquisitionSource object or empty if the file has no Acquisition group
        FilterTool              % McsCmosFilterSource object or empty if the file has no FilterTool group 
        SpikeExplorer           % McsCmosSpikeStream object or empty if the file has no SpikeExplorer group
        NetworkExplorer         % McsCmosNetworkExplorerSource object or empty if the file has no STAExplorer group
        SpikeSorter             % McsCmosSpikeSorterSource object or empty if the file has no SpikeSorter group
        UnknownDataSources = {} % List of McsCmosDataSource objects, one for each unknown data source in the file
    end
    
    methods
        function rec = McsCmosRecording(filename, recStruct, cfg)
        % Reads a single recording inside a HDF5 file.
        %
        % function rec = McsCmosRecording(filename, recStruct, cfg)
        %
        % Input:
        %   filename    -   (string) Name of the HDF5 file
        %   recStruct   -   The recording subtree of the structure 
        %                   generated by the h5info command
        %   cfg     -   (optional) configuration structure, contains one or
        %               more of the following fields:
        %               'dataType': The type of the data, can be one of
        %               'double' (default), 'single' or 'raw'. For 'double'
        %               and 'single' the data is converted to meaningful
        %               units, while for 'raw' no conversion is done and
        %               the data is kept in ADC units. This uses less
        %               memory than the conversion to double, but you might
        %               have to convert the data prior to analysis, for
        %               example by using the getConvertedData function.
        %               'timeStampDataType': The type of the time stamps,
        %               can be either 'int64' (default) or 'double'. Using
        %               'double' is useful for older Matlab version without
        %               int64 arithmetic.  
        %               'readUnknown': For CMOS-MEA data files, this 
        %               decides, whether data sources with unknown type are 
        %               read in. Default is false
            if exist('h5info')
                mode = 'h5';
            else
                mode = 'hdf5';
            end
            
            knownDataSources = containers.Map(...
            {...
                '650d88ce-9f24-4b20-ac2b-254defd12761', ...
                '2f8c246f-9bab-4193-b09e-03aefe17ede0', ...
                '941c8edb-78b3-4275-a5b2-6876cbcdeffc', ...
                'c6a37148-fa9e-42f2-9d38-eea0434851e2', ...
                '7263d1b7-f57a-42de-8f51-5d6326d22f2a', ...
            }, ...
            { ...
                'Acquisition', ...
                'FilterTool', ...
                'NetworkExplorer', ...
                'SpikeExplorer', ...
                'SpikeSorter', ...
            });
        
            cfg = McsHDF5.checkParameter(cfg, 'readUnknown', false);

            for gi = 1:length(recStruct.Groups)
                attributes = recStruct.Groups(gi).Attributes;
                for att = 1:length(attributes)
                    if McsHDF5.McsH5Helper.AttributeNameEquals(attributes(att), 'ID.TypeID', mode)
                        id = McsHDF5.McsH5Helper.AttributeValue(attributes(att), mode);
                        if knownDataSources.isKey(id)
                            name = knownDataSources(id);
                            if strcmp(name, 'Acquisition')
                                rec.(name) = McsHDF5.McsCmosAcquisitionSource(filename, recStruct.Groups(gi), cfg);
                            elseif strcmp(name, 'FilterTool')
                                rec.(name) = McsHDF5.McsCmosFilterSource(filename, recStruct.Groups(gi), cfg);
                            elseif strcmp(name, 'SpikeExplorer')
                                rec.(name) = McsHDF5.McsCmosSpikeStream(filename, recStruct.Groups(gi), cfg);
                            elseif strcmp(name, 'NetworkExplorer')
                                rec.(name) = McsHDF5.McsCmosNetworkExplorerSource(filename, recStruct.Groups(gi), cfg);
                            elseif strcmp(name, 'SpikeSorter')
                                rec.(name) = McsHDF5.McsCmosSpikeSorterSource(filename, recStruct.Groups(gi), cfg);
                            end
                        elseif cfg.readUnknown
                            source = McsHDF5.McsCmosDataSource(filename, recStruct.Groups(gi), cfg);
                            rec.UnknownDataSources = [rec.UnknownDataSources {source}];
                        end
                    end
                end
            end

            for li = 1:length(recStruct.Links)
                link = recStruct.Links(li);
                if strcmp(link.Name, 'Acquisition')
                    source = McsHDF5.McsCmosLinkedDataSource(link.Value{1}, link.Value{2});
                    rec.Acquisition = source;
                end
            end
        end
    end
end