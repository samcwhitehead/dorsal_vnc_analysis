% -------------------------------------------------------------------------
% function to convert Erica's annotation format (1 = input, 2 = output, 3 =
% both) into a different format for clustering
%
% currently, for an annotation vector of length N, converting it to a
% vector of length 2*N, where input and output each get their own dimension
% -------------------------------------------------------------------------
function mat_out = reshapeAnnotationMat(mat_in, target_dim) 
% -----------------
%% inputs
if ~exist('target_dim','var') || isempty(target_dim)
   target_dim = 2 ; % we assume that the matrix is MxN, where M is the number of neurons and N is the number of targets
end

% ------------------------------
%% perform reshaping
[nr, nc] = size(mat_in) ; 
if (target_dim == 1)
    % double matrix size along target dimension
    mat_out = repmat(mat_in, 2, 1) ; 
    
    % give "1" for either "1" or "3" in first half of matrix, then "1" for
    % either "2" or "3" in second half (first and second halves are input
    % and output, respectively)
    mat_out(1:nr,:) = (mat_out(1:nr,:) == 1) | (mat_out(1:nr,:) == 3) ; 
    mat_out((nr+1):end,:) = (mat_out((nr+1):end,:) == 2) | ...
        (mat_out((nr+1):end,:) == 3) ; 
elseif (target_dim == 2)
    % same as above
    mat_out = repmat(mat_in, 1, 2) ; 
    mat_out(:,1:nc) = (mat_out(:,1:nc) == 1) | (mat_out(:,1:nc) == 3) ; 
    mat_out(:,(nc+1):end) = (mat_out(:,(nc+1):end) == 2) | ...
        (mat_out(:,(nc+1):end) == 3) ; 
else
    fprintf('Invalid target dimension: %d \n',target_dim)
    keyboard
end


end