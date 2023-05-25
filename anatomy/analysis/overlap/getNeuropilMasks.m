% -------------------------------------------------------------------------
% function to get input/output masks for a given neuron, where masks here
% are the neuropil regions where the given neuron has inputs or outputs
%
% INPUTS:
%   - neuron_name: name of neuron that we're trying to get masks for
%       (string)
%   - annotation_vals: input/output annotation matrix. obtained from
%       "load_annotation_xlsx.m"
%   - neuron_labels: names of all neurons for which we have annotations.
%       obtained from "load_annotation_xlsx.m"
%   - target_labels: names of neuropil regions for which we have
%       annotations. obtained from "load_annotation_xlsx.m"
%   - neuropil_dir: directory struct containing registered/segmented
%       neuropil images/data
%   - output_type: string that determines form of mask output. either full
%       images ('bw') or list of coordinate indices ('coords')
%
% OUTPUTS:
%   - neuropil_masks: 1x2 cell array containing input and output masks for
%       current neuron (1 = input, 2 = output)
%
% -------------------------------------------------------------------------
function neuropil_masks = getNeuropilMasks(neuron_name, annotation_vals,...
    neuron_labels, target_labels, neuropil_dir, output_type, verbose_flag)
% ----------------------
%% inputs and params
% where to find segmented/registered neuropil images/data
if ~exist('neuropil_dir', 'var') || isempty(neuropil_dir)
    neuropil_dir = dir(['D:\Fly Imaging\Erica Interneuron Stacks\' ...
        'VNC neuropil\segmented_images\registered\*.mat']) ;
end
% do we want masks as full images or list of coordinate indices?
if ~exist('output_type', 'var') || isempty(output_type)
    output_type = 'coords' ;
end
% print when we're done?
if ~exist('verbose_flag', 'var') || isempty(verbose_flag)
    verbose_flag = false ;
end

% params
INPUT = 1 ;
OUTPUT = 2 ;
BOTH = 3 ;

fileSuffixStr = '' ; 
if strcmp(output_type,'coords') 
   fileSuffixStr = [fileSuffixStr '_coords'] ; 
end
% --------------------------------------
%% get annotations for current neuron
% first find row index for current neuron ("neuron_name")
match_idx = cellfun(@(y) strcmpi(neuron_name, y), neuron_labels) ;

% check that we just get one match
if sum(match_idx) ~= 1
    fprintf('Error finding match for %s \n', neuron_name)
    keyboard
end

% get annotation vals for this neuron
annotation_vals_curr = annotation_vals(match_idx,:) ;

% separate inputs/outputs
input_annotations = (annotation_vals_curr == INPUT) | ...
    (annotation_vals_curr == BOTH) ;
output_annotations = (annotation_vals_curr == OUTPUT) | ...
    (annotation_vals_curr == BOTH) ;
annotations_cell = {input_annotations, output_annotations} ;

% -------------------------------------------------
%% get masks
% initialize storage
neuropil_masks = cell(length(annotations_cell),1) ;

% loop over input/output
for k = 1:length(annotations_cell)
    % find targets for current neuron and condition (input vs output)
    targets_curr = target_labels(annotations_cell{k}) ;
    
    % initialize a variable for current mask
    mask_curr = [] ; 
    
    % loop over targets and build mask for neuron/condition
    for m = 1:length(targets_curr)
        % filename that we should search for (based on current target)
        search_fn = [targets_curr{m} fileSuffixStr '.mat'] ;
        
        % find neuropil mat file that matches current target
        file_idx = arrayfun(@(x) strcmpi(search_fn, x.name), neuropil_dir) ;
        
        % make sure we get a match
        if sum(file_idx) ~= 1
            fprintf('Error: cannot find a match for %s \n', search_fn)
            keyboard
        end
        
        % load data for current target neuropil
        np_data = importdata(fullfile(neuropil_dir(file_idx).folder, ...
            neuropil_dir(file_idx).name)) ; 
        
        % -----------------------------------------------------------
        % read out relevant data, depending on selected output type
        switch output_type
            case 'bw'
                % get total number of objects (binary images) for current 
                % neuropil
                N_obj = length(np_data.bw_cell) ; 
                
                % loop over neuropil objects
                for n = 1:N_obj
                    if isempty(mask_curr)
                        mask_curr = np_data.bw_cell{n} ; 
                    else
                        mask_curr = mask_curr || np_data.bw_cell{n} ; 
                    end
                end
            case 'coords'
                % get total number of objects (coordinate index lists) for 
                % current neuropil
                N_obj = length(np_data) ; 
                
                % loop over neuropil objects
                for n = 1:N_obj
                    if isempty(mask_curr)
                        mask_curr = np_data{n} ; 
                    else
                        mask_curr = union(mask_curr, np_data{n}) ; 
                    end
                end
            otherwise 
                fprintf('Invalid output type: %s \n', output_type)
                keyboard
        end
    end
    
    % ------------------------------------------------------------------
    % add current mask to output (having looped over all target labels)
    neuropil_masks{k} = mask_curr ; 
    
end

% ----------------
%% print update?
if verbose_flag
   fprintf('Obtained masks for neuron: %s \n', neuron_name) 
end

end