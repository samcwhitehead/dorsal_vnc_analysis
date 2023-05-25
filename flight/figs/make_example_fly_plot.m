%--------------------------------------------------------------------------
% script to generate a plot of a single fly (vs a control fly) to
% illustrate kinefly angles
%--------------------------------------------------------------------------
[mfilePath, ~, ~] = fileparts(mfilename('fullpath')) ; 
parentDirectory = fileparts(mfilePath) ;
dataPath = fullfile(parentDirectory, 'data') ;
savePath = fullfile(parentDirectory,'figs','output') ;

saveFlag = false ;
rawTracesFlag = false ;
smoothFlag = true ; 
excludeTrialsFlag = true ;

%---------------------------
% which fly?
driver = 'SS37246' ;
flyInd = 101 ; %154 ; % 101 ; 

%----------------------------
% make plot
h_main = plot_single_fly(dataPath, driver, flyInd, rawTracesFlag, ...
    smoothFlag, excludeTrialsFlag) ;
set(h_main, 'Renderer','painters')
%--------------------
% save results
if saveFlag
    savefig(h_main, fullfile(savePath, ...
        [driver '_' num2str(flyInd,'%03d') '.fig']))
    print(h_main, fullfile(savePath, ...
        [driver '_' num2str(flyInd,'%03d') '.svg']),'-dsvg')
    print(h_main, fullfile(savePath, ...
        [driver '_' num2str(flyInd,'%03d') '.png']),'-dpng','-r300')
end
