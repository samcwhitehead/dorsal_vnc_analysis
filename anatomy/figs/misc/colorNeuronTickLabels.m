% -------------------------------------------------------------------------
% function to color neuron name tick labels on a given axis
% -------------------------------------------------------------------------
function ax = colorNeuronTickLabels(ax, neuronColorStruct)
% ----------------------------------------------------------------------
% make sure axis tick labels are interpreted with 'tex' -- this lets us
% color them
ax.TickLabelInterpreter = 'tex' ;

% ---------------------------------------------------------------------
% start with x ticks. loop over labels. if neuron label is in our
% struct, change the tick label to that color
x_tick_labels = get(ax, 'XTickLabel') ;
for xn = 1:length(x_tick_labels)
    % current tick label
    tick_label_curr = strrep(x_tick_labels{xn},'-','_') ;
    
    % check that we have a color match
    if isfield(neuronColorStruct, tick_label_curr)
        % if we find a matching color, read it out and replace tick label
        neuronColor = neuronColorStruct.(tick_label_curr) ;
        tick_label_new = sprintf('{{\\color[rgb]{%f, %f, %f} %s}}', ...
            neuronColor(1), neuronColor(2), neuronColor(3), ...
            tick_label_curr) ;
    else
        % if we don't find a match, just get it to print black?
        tick_label_new = sprintf('{{\\color[rgb]{%f, %f, %f} %s}}', ...
            0.0, 0.0, 0.0, tick_label_curr) ;
    end
    
    % add label back to array; make sure to switch out hyphen and
    % underscore
    x_tick_labels{xn} = strrep(tick_label_new,'_','-') ;
    
end

% set axis to have new x tick labels
set(ax,'XTickLabel',x_tick_labels) ;

% --------------------------------------------------------------------
% now do the same with y labels
y_tick_labels = get(ax, 'YTickLabel') ;
for yn = 1:length(y_tick_labels)
    % current tick label
    tick_label_curr = strrep(y_tick_labels{yn},'-','_') ;
    
    % check that we have a color match
    if isfield(neuronColorStruct, tick_label_curr)
        % if we find a matching color, read it out and replace tick label
        neuronColor = neuronColorStruct.(tick_label_curr) ;
        tick_label_new = sprintf('{{\\color[rgb]{%f, %f, %f} %s}}', ...
            neuronColor(1), neuronColor(2), neuronColor(3), ...
            tick_label_curr) ;
        
    else
        % if we don't find a match, just get it to print black?
        tick_label_new = sprintf('{{\\color[rgb]{%f, %f, %f} %s}}', ...
            0.0, 0.0, 0.0, tick_label_curr) ;
    end
    
    % add y tick label back to array; switch hyphen and underscore
    y_tick_labels{yn} = strrep(tick_label_new,'_','-') ;
end

% set axis to have new y tick labels
set(ax,'YTickLabel',y_tick_labels) ;

end