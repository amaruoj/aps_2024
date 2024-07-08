function meanfield(folder)
% input:
% folder:   directory to the folder containing the .pbin and .pcd files

% output:
% meanfield for nozzle at given folder

% prep data once so we don't need to parse every time
pbin = dir(fullfile(folder,'*.pbin'));
pcd = dir(fullfile(folder,'*.pcd'));

% cylindrical data coordinate definition
nx = 751;
ntheta = 128;
nr = 151;
nvars = 5;

% prep chunk processing
nfiles = length(pcd);
chunk_size = 2000;  % each file is 0.1 acoustic time units, 2000 timesteps = 200
nchunks = ceil(nfiles / chunk_size);
global_data = zeros(ntheta, nr, nx, nvars); % init global mean for all chunks
test_name = folder(strfind(folder, 'pNozzle'):end);
test_name = test_name(9:12);

% chunk processing loop
for k = 0:nchunks - 1
    disp(['There are ', num2str(nfiles-(k*chunk_size)), ' files left to process in total!'])
    disp([num2str(nchunks-k), ' chunks remaining.'])
    remaining = nfiles - (k * chunk_size);
    if remaining < chunk_size   % handles edge case, less than perfect chunks
        nt = remaining;
    else
        nt = chunk_size;
    end

    chunk_data = zeros(ntheta, nr, nx, nvars);  % init running mean for chunk

    % data processing loop
    disp('reading from readBox.m...')
    for i = 1:nt
        % message every 100 timesteps to have an idea of where things go
        % wrong if there is a bug
        if ~mod(i,100)
            disp([num2str(i),' files processed!'])
        end
        % use readBox.m
        idx = i + k*chunk_size;
        currPcd = pcd(idx);
        [~, data] = readBox(fullfile(folder, pbin.name),fullfile(folder, currPcd.name));
        data = data';

        % reshaping data
        % data = permute(data, [2,3,1]);
        data = reshape(data, ntheta, nr, nx, nvars);

        % update running mean of this chunk
        chunk_data = chunk_data + (data ./ nt);
    end

    % calculate mean of chunk and save as a .mat
    tic
    disp('saving mean chunk matrix')
    filename = append('meanchunk_',test_name,'_',num2str(k+1));
    dirname = append('matrices_',test_name);
    out_dir = fullfile('..',dirname,'mean_data');
    if ~exist(out_dir,'dir')
        mkdir(out_dir);
    end
    vol_data = chunk_data;
    save(fullfile(out_dir, filename),'vol_data','-v7.3');
    toc

    % update global mean
    weight = nt / nfiles;   % not all chunks are the same size, some have more weight than others
    global_data = global_data + (weight .* chunk_data);
    
    % manage memory
    clear vol_data data chunk_data
end
% save .mat file of global mean
tic
disp('all done!! averaging all the chunks together...')
vol_data = global_data;
filename = append('meanfield_',test_name);
dirname = append('matrices_',test_name);
out_dir = fullfile('..',dirname,'mean_data');
save(fullfile(out_dir, filename),'vol_data','-v7.3');
toc
disp('finished! (*¯︶¯*)')
end