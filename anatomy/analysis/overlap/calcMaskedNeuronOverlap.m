% -------------------------------------------------------------------------
% function to calculate volume overlap between two neurons using the manual
% annotations of input/output sites along with the segmented/registered
% images of the VNC neuropil regions
%
% INPUTS:
%   - bw_cell: Nx1 cell array containing binary images of neurons 1 thru N
%   - mask_cell: Nx2 cell array containing masks for input and output
%       neuropils of neurons 1 thru N 
%           *mask_cell{i,1} = input mask for neuron i
%           *mask_cell{i,2} = output mask for neuron i
%   -bw_type: which representation of binary images we're using. either
%       'coords' or 'bw' (i.e. index list or full image, respectively)
%
% OUTPUTS:
%   - vols: NxN array containing volume overlaps.
%       *vols(i,j) = overlap in regions where neuron i has an output
%           and neuron j has an input (i.e. connection i -> j)
%   - overlap_coords: NxN cell array containing coordinates where images 
%       overlap (i.e. overlap_coords{i,j} is the list of voxel coordinate
%       indices corresponding to the overlap volume in (i,j)
% -------------------------------------------------------------------------
function [vols, overlap_coords] = calcMaskedNeuronOverlap(bw_cell, ...
    mask_cell, bw_type)
% inputs
if ~exist('bw_type', 'var') || isempty(bw_type)
   Ndim = numel(size(bw_cell{1})) ; 
   if Ndim > 2
       bw_type = 'bw' ; 
   else
       bw_type = 'coords' ; 
   end
end

% initialize storage for volume outputs
N_bw = length(bw_cell) ;
vols = zeros(N_bw,N_bw) ; 
overlap_coords = cell(N_bw,N_bw) ; 

% loop over the two neurons
for ii = 1:N_bw
    % we don't actually need to explicitly loop over jj here so long as
    % we're only dealing with 2 neurons at a time, but this should make the
    % function generalizable to more neurons
    for jj = 1:N_bw
        % skip when the indices match (i.e. same neuron)
        if ii == jj
            continue
        end
        
        % if we're here, it means we're comparing two different neurons.
        % use different method for overlap calculation depending on input
        % type
        switch bw_type
            case 'bw'
                % read out images
                bw_out = bw_cell{ii} ; % bw image where we want to look at outputs
                mask_out = mask_cell{ii, 2} ; % this is the mask for bw_out outputs
                
                bw_in = bw_cell{jj} ; % bw image where we want to look at inputs
                mask_in = mask_cell{jj, 1} ; % this is the mask for bw_in inputs
                
                % get overlap by taking intersection of all
                overlap = bw_out && mask_out && bw_in && mask_in ; 
                
                % get coordinates of overlap and total overlap volume
                overlap_coords{ii,jj} = find(overlap(:)) ; 
                vols(ii, jj) = numel(overlap_coords{ii,jj}) ; 
                
            case 'coords'
                % read out coordinates
                coords_out = bw_cell{ii} ; 
                mask_coords_out = mask_cell{ii,2} ; 
                
                coords_in = bw_cell{jj} ; 
                mask_coords_in = mask_cell{jj,1} ; 
                
                % mask coordinates
                coords_out_masked = intersect(coords_out, mask_coords_out) ; 
                coords_in_masked = intersect(coords_in, mask_coords_in) ; 
                
                % get overlap of masked coordinates
                overlap_coords{ii,jj} = intersect(coords_out_masked, ...
                    coords_in_masked) ; 
                vols(ii, jj) = numel(overlap_coords{ii,jj}) ; 
                
            otherwise 
                fprintf('Invalid input type -- skipping \n')
                continue
        end
    end
end
end