function pathOut = setPathParadigm
    [curr_path, ~, ~] = fileparts(mfilename('fullpath'));
    pathOut.home = curr_path;
    
    pathOut.info = fullfile(curr_path, 'scripts');
    
    pathOut.log = fullfile(curr_path, 'log');
    if (~exist(pathOut.log, 'dir'));
        mkdir(pathOut.log);
    end
    pathOut.results = fullfile(curr_path, 'results');
    if (~exist(pathOut.results, 'dir'));
        mkdir(pathOut.results);
    end
    
    pathOut.info = fullfile(curr_path, 'info');
    if (~exist(pathOut.info, 'dir'));
        mkdir(pathOut.info);
    end
end