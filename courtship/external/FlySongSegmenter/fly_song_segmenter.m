function fly_song_segmenter(song_path, varargin)
    % Convert Unix-style command line arguments to MATLAB-style.
    varargin = regexprep(varargin, '^-p$', 'params_path');
    varargin = regexprep(varargin, '^-s$', 'sample_range');
    varargin = regexprep(varargin, '^-c$', 'channel_num');
    
    p = inputParser;
    p.addRequired('song_path', @ischar);
    p.addParamValue('params_path', '', @ischar);
    p.addParamValue('sample_range', '', @(x)isempty(x) || ~isempty(regexp(x, '^[0-9]*:[0-9]*$', 'once')));
    p.addParamValue('channel_num', '', @(x)isempty(x) || ~isempty(regexp(x, '^[0-9]*$', 'once')));
    try
        p.parse(song_path, varargin{:});
    catch ME
        disp(ME.message);
        disp('Usage:');
        disp('    fly_song_segmenter song_path [-p params_path] [-s sample_range] [-c channel_num]');
        disp('e.g.');
        disp('    fly_song_segmenter song.wav -p params.m -s 1000:2000 -c 7');
        return
    end
    
    if isempty(p.Results.sample_range)
        sample_range = [];
    else
        tokens = regexp(p.Results.sample_range, '^([0-9]*):([0-9]*)$', 'tokens');
        sample_range = [str2double(tokens{1}{1}) str2double(tokens{1}{2})];
    end
    
    [~, ~, ext] = fileparts(p.Results.song_path);
    if strcmp(ext, '.daq')
        FlySongSegmenterDAQ(p.Results.song_path, str2num(p.Results.channel_num), sample_range, p.Results.params_path);
    elseif strcmp(ext, '.wav')
        FlySongSegmenterWAV(p.Results.song_path, str2num(p.Results.channel_num), sample_range, p.Results.params_path);
%        if(~isempty(p.Results.channel_num))
%          error('-c channel_num only valid with .daq files');
%        end
%        if(isempty(sample_range))
%          [data, Fs] = wavread(p.Results.song_path);
%        else
%          [data, Fs] = wavread(p.Results.song_path,sample_range);
%        end
%        process_song_data(p.Results.song_path, data, Fs, p.Results.params_path);
%        %[data, Sines, Pulses, Params]=FlySongSegmenter(data, [], p.Results.params_path, Fs);
%        %save([p.Results.song_path(1:end-4) '.mat'], 'data','Sines','Pulses','Params','-v7.3');
    elseif strcmp(ext, '.au')
        if(~isempty(p.Results.channel_num))
          error('-c channel_num only valid with .daq files');
        end
        if(isempty(sample_range))
          [data, Fs] = auread(p.Results.song_path);
        else
          [data, Fs] = auread(p.Results.song_path,sample_range);
        end
        process_song_data(p.Results.song_path, data, Fs, p.Results.params_path);
        %[data, Sines, Pulses, Params]=FlySongSegmenter(data, [], p.Results.params_path, Fs);
        %save([p.Results.song_path(1:end-4) '.mat'], 'data','Sines','Pulses','Params','-v7.3');
    else
        error('Unknown song file format: %s', ext);
    end
end


function process_song_data(song_path, song_data, Fs, params_path)
    [parentDir, song_name, ~] = fileparts(song_path);
    outfile  = fullfile(parentDir, [song_name '.mat']);
    if ~exist(outfile, 'file')
%        if ~isempty(sample_range)
%            song_data = song_data(sample_range(1):sample_range(2));
%        end
        [Data, Sines, Pulses, Params] = FlySongSegmenter(song_data, [], params_path); %#ok<NASGU,ASGLU>
        save(outfile, 'Data','Sines','Pulses','Params','-v7.3')
        clear song_data Data Sines Pulses Params;
    else
        fprintf('File %s exists. Skipping.\n', outfile)
    end
end
