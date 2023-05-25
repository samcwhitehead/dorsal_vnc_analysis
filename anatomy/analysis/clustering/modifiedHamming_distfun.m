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
function D = modifiedHamming_distfun(zi, zj)

D = zeros(size(zj,1),1) ;

for k = 1:size(zj,1)
    val_sum = 0 ; 
    for m = 1:length(zi)
        x = zi(m) ;
        y = zj(k, m) ;
        if (x == y)
            val_add = 1 ;
        elseif ((x == 3) && (y > 0)) || ...
                ((x > 0) && (y == 3))
            val_add = 0.5 ;  % should this just be zero?
        else
            val_add = 0 ;
        end
        
        val_sum = val_sum + val_add ;
    end
    
    D(k) = 1 -  (val_sum / length(zi)) ; % normalize by total length
end

end