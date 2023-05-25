% -------------------------------------------------------------------------
% function to normalize matrix of anatomical overlap
% -------------------------------------------------------------------------
function overlap_data_norm = ...
    normalize_overlap_data(overlap_data, vox_num_data1, vox_num_data2, normType)
% ---------------------
%% inputs
if ~exist('vox_num_data1','var') || isempty(vox_num_data1)
    vox_num_data1 = [] ; % normally a vector conainting voxel count for each stack
    normType = 'max' ;
end
if ~exist('vox_num_data2','var') || isempty(vox_num_data2)
    vox_num_data2 = [] ; % normally a vector conainting voxel count for each stack
    normType = 'max' ;
end
if ~exist('normType','var') || isempty(normType)
    normType = 'max' ; % method of normalization. can be 'max', 'sym_size', 'asym_size'
end

% see if input is table or array
if istable(overlap_data)
    tableFlag = true;
else
    tableFlag = false ;
end

% --------------------------------------------------------
%% remove any deleted entries from overlap table/matrix
if tableFlag 
    % get row and column names
    rowNames = overlap_data.Properties.RowNames ; 
    colNames = overlap_data.Properties.VariableNames ; 
    
    % find any rows or columns with "deleted" in the name
    deleted_row_idx = cellfun(@(y) contains(y, 'delete', 'IgnoreCase',1),...
        rowNames) ; 
    deleted_col_idx = cellfun(@(y) contains(y, 'delete', 'IgnoreCase',1),...
        colNames) ; 
    
    % keep only rows/cols without the "deleted" label
    overlap_data = overlap_data(~deleted_row_idx, ~deleted_col_idx) ; 
end

% ---------------------------------
%% initialize normalized output
[N_stacks1, N_stacks2] = size(overlap_data) ;
% switch numbers on voxel size vector? I think my overlap matrices might
% not have consistent order
if (size(vox_num_data1,1) == N_stacks2) && ...
        (size(vox_num_data2,1) == N_stacks1)
    temp = vox_num_data1 ;
    vox_num_data1 = vox_num_data2 ;
    vox_num_data2 = temp ;
end

% initialize
if tableFlag
    overlap_norm_init = nan(size(overlap_data)) ;
    overlap_data_norm = array2table(overlap_norm_init) ;
    overlap_data_norm.Properties.RowNames = ...
        overlap_data.Properties.RowNames ;
    overlap_data_norm.Properties.VariableNames = ...
        overlap_data.Properties.VariableNames ;
else
    overlap_data_norm = nan(size(overlap_data)) ;
end
% --------------------------------------------------------
%% deal with overlap mat diagonal and symmetrize if need be
if tableFlag
    % this is dumb, but if it's a table it needs curly brackets...
    if (N_stacks1 == N_stacks2)
        for i = 1:N_stacks1
            overlap_data{i,i} = nan ;
        end
        for ii = 1:N_stacks1
            for jj = 1:ii
                overlap_data{ii,jj} = overlap_data{jj,ii} ;
                overlap_data{jj,ii} = overlap_data{jj,ii} ;
            end
        end
    end
    
else
    % ...if it's a matrix, parentheses
    if (N_stacks1 == N_stacks2)
        for i = 1:N_stacks1
            overlap_data(i,i) = nan ;
        end
        for ii = 1:N_stacks1
            for jj = 1:ii
                overlap_data(ii,jj) = overlap_data(jj,ii) ;
                overlap_data(jj,ii) = overlap_data(jj,ii) ;
            end
        end
        
    end
end
% ------------------------------------------------------
%% allow for a few different types of normalization
switch normType
    case 'max'
        % normalize so maximum overlap is 1. does not take into account
        % overall size of expression pattern
        if tableFlag
            max_overlap = nanmax(overlap_data{:,:},[],'all') ;
            overlap_data_norm{:,:} = overlap_data{:,:}./max_overlap ;
        else
            max_overlap = nanmax(overlap_data(:)) ;
            overlap_data_norm = overlap_data./max_overlap ;
        end
    case 'sym_size'
        % divide the voxel overlap number by the average of the voxel
        % counts for the two stacks being compared. this will lead to a
        % symmetric overlap matrix
        % -----------------------------------------
        % loop over rows
        for r = 1:N_stacks1
            % -------------------
            % vox num for rows:
            % -------------------
            if istable(vox_num_data1) && tableFlag
                % if both overlap data and voxel count are in table
                % format, we can be careful and check for matching line
                % names
                name1_curr = overlap_data.Properties.RowNames{r}  ;
                match_idx1 = cellfun(@(y) strcmp(y, name1_curr), ...
                    vox_num_data1.lineName) ;
                % cludgy fix for hyphen vs underscore bullshit
                if (sum(match_idx1) == 0)
                    match_idx1 = cellfun(@(y) strcmp(strrep(y,'-','_'),...
                        name1_curr), vox_num_data1.lineName) ;
                end
                % if we still can't find a match, there's a problem
                if (sum(match_idx1) ~= 1)
                    fprintf('Could not find a match for line: %s \n', ...
                        name1_curr)
                    keyboard
                end
                % pull out matching entry
                vox_size1 = vox_num_data1.voxCount(match_idx1) ;
            elseif istable(vox_num_data1) && ~tableFlag
                % if overlap data is not in table form, just need to
                % read from correct field
                vox_size1 = vox_num_data1.voxCount(r) ;
            else
                % if everything is arrays, standard read
                vox_size1 = vox_num_data1(r) ;
            end
            
            % ----------------------------------------------
            % loop over columns
            for c = 1:N_stacks2
                % ------------------
                % vox num for cols:
                % ------------------
                if istable(vox_num_data2) && tableFlag
                    % if both overlap data and voxel count are in table
                    % format, we can be careful and check for matching line
                    % names
                    name2_curr = overlap_data.Properties.VariableNames{c} ;
                    match_idx2 = cellfun(@(y) strcmp(y, name2_curr), ...
                        vox_num_data2.lineName) ;
                    % cludgy fix for hyphen vs underscore bullshit
                    if (sum(match_idx2) == 0)
                        match_idx2 = cellfun(@(y) strcmp(strrep(y,'-','_'),...
                            name2_curr), vox_num_data2.lineName) ;
                    end
                    % if we still can't find a match, there's a problem
                    if (sum(match_idx2) ~= 1)
                        fprintf('Could not find a match for line: %s \n', ...
                            name2_curr)
                        keyboard
                    end
                    % pull out matching entry
                    vox_size2 = vox_num_data2.voxCount(match_idx2) ;
                elseif istable(vox_num_data2) && ~tableFlag
                    % if overlap data is not in table form, just need to
                    % read from correct field
                    vox_size2 = vox_num_data2.voxCount(c) ;
                else
                    % if everything is arrays, standard read
                    vox_size2 = vox_num_data2(c) ;
                end
                
                % ------------------------------------------------------
                % use these results to get mean voxel size and use to
                % normalize
                mean_size = nanmean([vox_size1, vox_size2]) ;
                overlap_data_norm{r, c} = overlap_data{r, c}./mean_size ;
            end
        end
    case 'asym_size'
        % divide the voxel overlap number by the voxel count of the stack
        % corresponding to the row. This will be asymmetric
        % ------------------------------------
        % loop over rows
        for r = 1:N_stacks1
            % ----------------------
            % vox num for rows
            % ----------------------
            if istable(vox_num_data1) && tableFlag
                % if both overlap data and voxel count are in table
                % format, we can be careful and check for matching line
                % names
                name1_curr = overlap_data.Properties.RowNames{r}  ;
                match_idx1 = cellfun(@(y) strcmp(y, name1_curr), ...
                    vox_num_data1.lineName) ;
                % cludgy fix for hyphen vs underscore bullshit
                if (sum(match_idx1) == 0)
                    match_idx1 = cellfun(@(y) strcmp(strrep(y,'-','_'),...
                        name1_curr), vox_num_data1.lineName) ;
                end
                % if we still can't find a match, there's a problem
                if (sum(match_idx1) ~= 1)
                    fprintf('Could not find a match for line: %s \n', ...
                        name1_curr)
                    keyboard
                end
                % pull out matching entry
                vox_size1 = vox_num_data1.voxCount(match_idx1) ;
            elseif istable(vox_num_data1) && ~tableFlag
                % if overlap data is not in table form, just need to
                % read from correct field
                vox_size1 = vox_num_data1.voxCount(r) ;
            else
                % if everything is arrays, standard read
                vox_size1 = vox_num_data1(r) ;
            end
            % -------------------------------
            % loop over columns
            for c = 1:N_stacks2
               
                % divide by just size corresponding to row image
                overlap_data_norm{r, c} = overlap_data{r, c}./vox_size1 ;
            end
        end
    otherwise
        disp('Invalid normalization type')
        overlap_data_norm = [] ;
        return
end

end