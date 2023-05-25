%--------------------------------------------------------------------------
% functions to calculate distances between connectivity manual annotations.
% the annotations take the form of a vector where each dimension is a
% different brain/VNC region. The possible values for each entry are:
%       0 = no connection
%       1 = input
%       2 = output
%       3 = both input and output
% Going to start with a modified Hamming distance
%
%--------------------------------------------------------------------------
function D = outputOnly_distfun(zi, zj)

D = zeros(size(zj,1),1) ;

for k = 1:size(zj,1)
    for m = 1:length(zi)
        x = zi(m) ;
        y = zj(k, m) ;
        if (x == 2) && (y == 2)
            val_add = 0 ;
        elseif ((x == 3) && (y == 2)) || ...
                ((x == 2) && (y == 3))
            val_add = 0.5 ;  % should this just be zero?
        else
            val_add = 1 ;
        end
        
        D(k) = D(k) + val_add ;
    end
    
    D(k) = D(k) / length(zi) ; % normalize by total length
end

end