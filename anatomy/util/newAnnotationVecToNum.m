% -------------------------------------------------------------------------
% function to convert new annotation representation (1x2 or 1x4 vector per
% neuron) into the 0-3 system previously used for the old annotations. this
% should help make the code i've already written easy to use with the new
% annotation values
%
% NB: annot_vec should always be 1x2. If we're using ipsi- and
% contralateral annotations separately, we'll ultimately have a vector of
% the form:
%   [ipsi_input, contra_input, ipsi_output, contra_output]
% but we can just call this function on vec([1,3]) and vec([2,4]) to get
% separate inputs and output values for ipsi- and contralateral
% annotations, respectively
%
% -------------------------------------------------------------------------
function annot_val = newAnnotationVecToNum(annot_vec, capitalScale, ...
    baseVal, rmWeakAnnotFlag)
% ------------------------
% inputs and params
if ~exist('capitalScale','var') || isempty(capitalScale)
    capitalScale = 2 ;
end
if ~exist('baseVal','var') || isempty(baseVal)
    baseVal = 1 ;
end
if ~exist('rmWeakAnnotFlag','var') || isempty(rmWeakAnnotFlag)
    % remove weak (i.e. lower case) annotations?
    rmWeakAnnotFlag = false ;
end

% possible values for annot_vec entries -- note that this depends on
% whether or not we're including weighted (i.e. upper vs lower case)
% annotations
if rmWeakAnnotFlag
    % in this case we only have "1" or "0"
    low_val = 1 ;
    high_val = 1 ;
else
    % in this case, we have different values corresponding to different
    % annotation "strengths"
    high_val = capitalScale*baseVal ;  % corresponds to upper case letter in new annotations
    low_val = baseVal ; % corresponds to lower case letter in new annotations
end

% ------------------------------------------------------------------
% check cases for values of annot_vec
if ((annot_vec(1) == high_val) && (annot_vec(2) >= low_val)) || ...
        ((annot_vec(1) >= low_val) && (annot_vec(2) == high_val))
    % strong input and output
    annot_val = 3.0 ;
elseif (annot_vec(1) == low_val) && (annot_vec(2) == low_val)
    % weak input and output
    annot_val = 2.5 ;
elseif (annot_vec(1) == high_val) && (annot_vec(2) < low_val)
    % strong input
    annot_val = 1 ;
elseif (annot_vec(1) == low_val) && (annot_vec(2) < low_val)
    % weak input input
    annot_val = 0.5 ;
elseif (annot_vec(1) < low_val) && (annot_vec(2) == high_val )
    % strong output
    annot_val = 2 ;
elseif (annot_vec(1) < low_val) && (annot_vec(2) == low_val)
    % weak output
    annot_val = 1.5 ;
elseif (annot_vec(1) < low_val) && (annot_vec(2) < low_val)
    % nothing
    annot_val = 0 ;
else
    % don't know how this would come up
    keyboard
end

end