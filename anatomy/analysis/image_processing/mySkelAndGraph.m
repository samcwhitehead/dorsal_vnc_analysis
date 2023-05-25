% -------------------------------------------------------------------------
% function to convert 3D binary volume into a graph representation using
% the skel2graph3d package from matlab file exchange. this function just
% mimics the script "Test_Skel2Graph3D.m" from that same package.
%
% NB: this function takes a binary volume as input, and performs
% skeletonization within
%
%{
dataPath = 'D:\Fly Imaging\Erica Interneuron Stacks\IN\binarized_non_sym\' ; 
fn = 'NH1_8Balphat3_3093_SS47214-20180202_25_H1_bw.mat' ; 
bwMat = importdata(fullfile(dataPath, fn)) ; 

[skel_struct, h_main] = mySkelAndGraph(bwMat,0, true) ;
%}
% -------------------------------------------------------------------------
function [skel_struct, h_main] = mySkelAndGraph(bwMat,THR, plotFlag)
% ----------------------
%% inputs
if ~exist('THR','var') || isempty(THR)
   THR = 0 ; % minimum branch length for skeleton, used to filter artifacts 
end
if ~exist('plotFlag','var') || isempty(plotFlag)
   plotFlag = false ; % make plot of results?
end

% ----------------------------------------------
%% first round (skeletonize and get graph rep)
skel_init = Skeleton3D(bwMat);

w = size(skel_init,1);
l = size(skel_init,2);
h = size(skel_init,3);

% initial step: condense, convert to voxels and back, detect cells
[~,node_init,link_init] = Skel2Graph3D(skel_init,THR);

% total length of network
wl = sum(cellfun('length',{node_init.links}));

skel_final = Graph2Skel3D(node_init,link_init,w,l,h);
[~,node_final,link_final] = Skel2Graph3D(skel_final,THR);

% calculate new total length of network
wl_new = sum(cellfun('length',{node_final.links}));

% ------------------------------------------------------------------------
%% iterate the same steps until network length changed by less than 0.5%
while (wl_new~=wl)

    wl = wl_new;   
    
     skel_final = Graph2Skel3D(node_final,link_final,w,l,h);
     [~,node_final,link_final] = Skel2Graph3D(skel_final,THR);

     wl_new = sum(cellfun('length',{node_final.links}));
       
end

% ----------------------------------------------------------------------
%% store info in one structure 
skel_struct = struct() ;
skel_struct.skel = skel_final ;
skel_struct.node = node_final ;
skel_struct.link = link_final ;

% ----------------------------------------------------------------------
%% display result?
if plotFlag
    h_main = figure();
    hold on;
    for i=1:length(node_final)
        x1 = node_final(i).comx;
        y1 = node_final(i).comy;
        z1 = node_final(i).comz;
        
        if(node_final(i).ep==1)
            ncol = 'c';
        else
            ncol = 'y';
        end
        
        for j=1:length(node_final(i).links)    % draw all connections of each node
            if(node_final(node_final(i).conn(j)).ep==1)
                col='k'; % branches are black
            else
                col='r'; % links are red
            end
            if(node_final(i).ep==1)
                col='k';
            end
            
            
            % draw edges as lines using voxel positions
            for k=1:length(link_final(node_final(i).links(j)).point)-1
                [x3,y3,z3] = ind2sub([w,l,h], ...
                    link_final(node_final(i).links(j)).point(k));
                [x2,y2,z2] = ind2sub([w,l,h], ...
                    link_final(node_final(i).links(j)).point(k+1));
                line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
            end
        end
        
        % draw all nodes as yellow circles
        plot3(y1,x1,z1,'o','Markersize',9,...
            'MarkerFaceColor',ncol,...
            'Color','k');
    end
    axis image;axis off;
    set(gcf,'Color','white');
    drawnow;
    view(-17,46);
else
    h_main = [] ; 
end

end