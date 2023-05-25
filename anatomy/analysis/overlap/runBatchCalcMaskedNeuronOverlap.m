% -------------------------------------------------------------------------
% script to calculate masked neuron overlap over multiple neuron sets
% -------------------------------------------------------------------------
%% path info
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ;
figDirectory = fileparts(mfilePath) ;
parentDirectory = fileparts(figDirectory) ;
dataRoot = fullfile(parentDirectory, 'data', 'processed_anatomy_data') ;
savePath = fullfile(dataRoot, 'overlap_calculations','masked') ;

% options
overWriteFlag = false ;

% neuron and data type
% neuronTypes = {'IN', 'DN', 'MN', 'VUM'} ;
neuronTypes = {'IN', 'MN'} ;
dataMode = 'coords' ;

% -------------------------------------------------------------------
%% get all combinations of neurons
[x1, x2] = meshgrid(1:length(neuronTypes), 1:length(neuronTypes));
dataTypeList1 = neuronTypes(x1(:)) ;
dataTypeList2 = neuronTypes(x2(:)) ;

% -------------------------------------------------------------------
%% loop over neuron combinations and get overlap data
for ii = 1:length(dataTypeList1)
    % current dataType1
    dataType1 = dataTypeList1{ii} ;
    
    for jj = 1:length(dataTypeList2)
        % current dataType2
        dataType2 = dataTypeList2{jj} ;
        
        % check if we've already calculated this
        save_fn = sprintf('%s_%s_overlap_table.mat', dataType1, dataType2) ;
        if exist(fullfile(savePath, save_fn),'file') && ~overWriteFlag
            fprintf('Already calculated overlap for %s -- skipping \n',...
                save_fn)
            continue
        end
        
        % if we get here, need to calculate overlap for current neuron
        % combination. do so
        overlap_data = batchCalcMaskedNeuronOverlap(dataRoot, ...
            dataType1, dataType2, savePath, dataMode) ;
    end
end