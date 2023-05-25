% -------------------------------------------------------------------------
% quick/dirty function to parse motor neuron names in new format (e.g.
% DLMNab to 'DLMNa/b')
% -------------------------------------------------------------------------
function mn_names_out = parse_mn_names(mn_names_in)
% ------------------------------
%% define translation struct
trans_struct = struct() ;

trans_struct.DLMNab = 'DLMNa/b' ;
trans_struct.DLMNc = 'DVMNc' ;
trans_struct.DLMNd = 'DVMNd' ;
trans_struct.DLMNe = 'DVMNe' ;
trans_struct.DLMNf = 'DVMNf' ;

trans_struct.DVMN1a = 'DVMN1a' ;
trans_struct.DVMN1b = 'DVMN1b' ;
trans_struct.DVMN1b2a = 'DVMN1b,2a' ;
trans_struct.DVMN2b = 'DVMN2b' ;
trans_struct.DVMN2c = 'DVMN2c' ;
trans_struct.DVMN3a = 'DVMN3a' ;
trans_struct.DVMN3b = 'DVMN3b' ;

trans_struct.b1 = 'b1 MN' ;
trans_struct.b2 = 'b2 MN' ;

trans_struct.i1 = 'i1 MN' ;
trans_struct.i2 = 'i2 MN' ;

trans_struct.iii1 = 'iii1 MN' ;
trans_struct.iii3 = 'iii3 MN' ;

trans_struct.hg1 = 'hg1 MN' ;
trans_struct.hg2 = 'hg2 MN' ;
trans_struct.hg3 = 'hg3 MN' ;

trans_struct.ps1 = 'ps1 MN' ;

trans_struct.tp1 = 'tp1 MN' ;
trans_struct.tp2 = 'tp2 MN' ;
trans_struct.tpN = 'tpN MN' ;

trans_struct.tt = 'tt MN' ;

% --------------------------------------
%% use struct to get parsed neuron name
if iscell(mn_names_in)
    % loop through neuron names if we are given multiple names
    mn_names_out = mn_names_in ;
    for k = 1:length(mn_names_in)
        mn_names_out{k} = trans_struct.(mn_names_in{k}) ;
    end
else
    % if it's just one entry, output string
    mn_names_out = trans_struct.(mn_names_in) ;
end
end