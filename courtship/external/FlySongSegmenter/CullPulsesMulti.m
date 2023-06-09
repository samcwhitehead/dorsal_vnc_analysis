function CullPulsesMulti(folder,pulseInfo_name,Lik_name,Lik_type,minRange,maxRange)
%USAGE CullPulsesMulti(folder,pulseInfo_name,Lik_type,range)
%e.g.
%Lik_name = 'Lik_pulseInfo' or 'Lik_pulseInfo2'
%Lik_type = 'LLR_fh'
%minRange = 0
%maxRange = 'max'

pI_name = char(pulseInfo_name);
Lik_name = char(Lik_name);
Lik_type = char(Lik_type);

if strcmp(folder(end),'/') == 0
    folder = [folder '/'];
end
dir_list = dir(folder);
file_num = length(dir_list);

for y = 1:file_num
    file = dir_list(y).name; %pull out the file name
    [~,root,ext] = fileparts(file);
    path_file = [folder file];
    TG = strcmp(ext,'.mat');
    
    if TG == 1
        
        fprintf([file '\n']);

        W = who('-file',path_file);
        varstruc =struct;
        load(path_file);
        
        for ii = 1:numel(W)
            varstruc.(W{ii}) = eval(W{ii});
        end
        
        pIData = load(path_file,pI_name);
        pI_data = pIData.(pI_name);
       
        LikData = load(path_file,Lik_name);
        a = fieldnames(LikData);
        LikData = getfield(LikData,char(a));
        Lik_data = LikData.(Lik_type);
        
        
        if nargin < 6
            maxRange = max(Lik_data) + 1;
        end
        
        culled_pulseInfo = ModelCullPulses(pI_data,Lik_data,[minRange maxRange]);
        
        if strcmp(pI_name,'pulseInfo') == 1

            varstruc.culled_pulseInfo = culled_pulseInfo;
            varstruc.culled_pulseInfo.variables.pulseInfo_ver = pI_name;
            varstruc.culled_pulseInfo.variables.Lik_ver = Lik_name;
            varstruc.culled_pulseInfo.variables.range = [minRange maxRange];
            varstruc.culled_pulseInfo.variables.date = date;
            varstruc.culled_pulseInfo.variables.time = clock;
                        
        elseif strcmp(pI_name,'pulseInfo2') == 1
            varstruc.culled_pulseInfo2 = culled_pulseInfo;
            varstruc.culled_pulseInfo2.variables.pulseInfo_ver = pI_name;
            varstruc.culled_pulseInfo2.variables.Lik_ver = Lik_name;
            varstruc.culled_pulseInfo2.variables.range = [minRange maxRange];
            varstruc.culled_pulseInfo2.variables.date = date;
            varstruc.culled_pulseInfo2.variables.time = clock;

        end
 
        
        varstruc.culled_pulseInfo = culled_pulseInfo;   
        varstruc.culled_pulseInfo.variables.pulseInfo_ver = pI_name;
        varstruc.culled_pulseInfo.variables.Lik_ver = Lik_name;
        varstruc.culled_pulseInfo.variables.range = [minRange maxRange];
        varstruc.culled_pulseInfo.variables.date = date;
        varstruc.culled_pulseInfo.variables.time = clock;
        
        save(path_file,'-struct','varstruc','-mat')%save all variables in original file
       
    end
end


