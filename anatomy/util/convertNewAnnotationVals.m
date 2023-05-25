% -------------------------------------------------------------------------
% function to  convert new annotation characters to integers. New
% annotations have a pair of characters for each neuron/neuropil
% combination with the following meanings:
%   - d or D: weak or strong input (dendrites)
%   - a or A: weak or strong output (axon)
%   - m or M: weak or strong mixed input and output (mixed)
%   - p or P: weak or strong partitioned input and output (partitioned)
%   - "-": nothing
%
% NB: not differentiating between partitioned and mixed at this point
%
% NB2: we now have ipsi and contralateral annotations. so output matrix can
% potentially be Nx(4*M) where M is the number of neuropil regions (we get
% a factor of two for input vs output, and another for contra vs ipsi)
% -------------------------------------------------------------------------
function [vals_out, vals_out_3D] = convertNewAnnotationVals(vals_in, ...
    capitalFactor, useIpsiContraFlag, rmWeakAnnotFlag, baseVal)
% -----------------------------
%% inputs/params
if ~exist('capitalFactor','var') || isempty(capitalFactor)
    % difference between upper and lowercase letters in new annotation
    % format -- Erica uses e.g. D and d to signify major and minor presence
    % of dendrites in a region. This factor determines the numerical
    % difference between "D" and "d"
    capitalFactor = 2 ;
end
if ~exist('useIpsiContraFlag','var') || isempty(useIpsiContraFlag)
    % differentiate between ipsi- and contralateral annotations (if true),
    % or just  lump them together (if false)?
    useIpsiContraFlag = true ;
end
if ~exist('rmWeakAnnotFlag','var') || isempty(rmWeakAnnotFlag)
    % remove weak (i.e. lower case) annotations?
    rmWeakAnnotFlag = false ;
end
if ~exist('baseVal','var') || isempty(baseVal)
    % base value to indicate presence of an input or output
    baseVal = 1 ; 
end

% list of characters signifying presence of input or output
input_list = {'d','D','m','M','p','P'} ; 
output_list =  {'a','A','m','M','p','P'} ;

% get number of neurons and number of neuropils
[N_neurons, N_neuropils] = size(vals_in) ;

% if we're removing all weak (lower case annotations) make sure
% capitalFactor > 1 (this will allow us to differentiate between weak and
% strong)
if rmWeakAnnotFlag && (capitalFactor <= 1)
   capitalFactor = 2 ;  
end

% -----------------------------------------------
%% loop over entries in vals_in and convert
% first, initialize output. start by putting ipsi/contra and in/out as 3rd
% dimension, so the order is [ipsi_in, contra_in, ipsi_out, contra_out]
vals_out = zeros(N_neurons, N_neuropils, 4) ;

% loop over rows (neurons) first, then columns (neuropils)
for rr = 1:N_neurons
    for cc = 1:N_neuropils
        % read out annotation for current string
        vals_curr = vals_in{rr,cc} ;
        
        % if the entry is empty, then there are no inputs or outputs.
        % otherwise, need to fill in vals_out
        if ~isempty(vals_curr)
            % loop over ipsi/contralateral
            for k = 1:length(vals_curr)
                % current ipsi or contra value
                val = vals_curr(k) ; 
                
                % ---------------------------------------------------
                % fill in values for INPUTS in vals_out
                if ismember(val,input_list)
                    % if current annotation value indicates input, get
                    % proper index for entry in vals_out
                    idx_in = sub2ind(size(vals_out), rr, cc, 1 + (k-1)) ; 
                    
                    % enter value, but check if it is upper or lower case
                    if val == upper(val)
                        vals_out(idx_in) = capitalFactor*baseVal ;
                    else
                        vals_out(idx_in) = baseVal ;
                    end
                end
                
                % ---------------------------------------------------
                % fill in values for Outputs in vals_out
                if ismember(val,output_list)
                    % if current annotation value indicates input, get
                    % proper index for entry in vals_out
                    idx_out = sub2ind(size(vals_out), rr, cc, 3 + (k-1)) ; 
                    
                    % enter value, but check if it is upper or lower case
                    if val == upper(val)
                        vals_out(idx_out) = capitalFactor*baseVal ;
                    else
                        vals_out(idx_out) = baseVal ;
                    end
                end
            end
        end
    end
end

% ---------------------------------------------------------------
%% use the option to take only strong (upper case) annotations?
if rmWeakAnnotFlag
   % weak annotations should equal baseVal and strong ones should equal
   % capitalFactor*baseVal. So just threshold by capitalFactor*baseVal to
   % retain only strong connections
   vals_out_logical = (vals_out > baseVal) ; % now we have strong annotations = 1, weak = 0
   vals_out = double(vals_out_logical) ; % convert from logical back to double
   
end

% ---------------------------------------------------------------
%% reshape matrix
% if differentiating between ipsi and contralateral, just unwrap third
% dimension. otherwise, need to take max of sides and combine
if useIpsiContraFlag
    % store copy of 3D array
    vals_out_3D = vals_out ; 
    % but convert main output to 2D array
    vals_out = reshape(vals_out, [N_neurons, 4*N_neuropils]) ; 
else
    % store copy of 3D array
    vals_out_3D = cat(3, max(vals_out(:,:,1:2),[],3), ...
        max(vals_out(:,:,3:4),[],3)); 
    % but convert main output to 2D array
    vals_out = [max(vals_out(:,:,1:2),[],3), max(vals_out(:,:,3:4),[],3)] ; 
end

end
