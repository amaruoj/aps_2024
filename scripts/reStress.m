function reStress(folder)
% input:
% folder:   directory to the folder containing the .pbin and .pcd files

% output:
% saved matrices showing the reynolds stress of each jet

% prep data once so we don't need to parse every time
pbin = dir(fullfile(folder,'*.pbin'));
pcd = dir(fullfile(folder,'*.pcd'));

% cylindrical data coordinate definition
nx = 751;
ntheta = 128;
nr = 151;

% prep data processing loop
vars = {'U','V','W'};
nvars = length(vars);
nfiles = length(pcd);
test_name = folder(strfind(folder, 'pNozzle'):end);
test_name = test_name(9:12);
stress = zeros(ntheta,nr,nx,6);

% cylindrical data coordinate definition
nx = 751; ntheta = 128; nr = 151;
x = linspace(0,30, nx)'; r = linspace(0,6, nr)';
theta = linspace(0, 2*pi-2*pi/128, ntheta)';
[~, THETA, ~] = meshgrid(r, theta, x);

% step 1: retrieve meanfield data
meanfield = load(fullfile('..',append('matrices_',test_name), ...
    'mean_data',append('meanfield_',test_name))).vol_data;
v = meanfield(:,:,:,2); w = meanfield(:,:,:,3);   % axial experiences no change
meanfield(:,:,:,2) = v .* cos(THETA) + w .* sin(THETA);  % radial
meanfield(:,:,:,3) = -v .* sin(THETA) + w .* cos(THETA); % azimuthal

% data processing loop, no chunks for this analysis
tic
disp('reading from readBox.m...')
for i = 1:nfiles
    if ~mod(i,100)
        disp([num2str(i),' files processed!'])
    end
    % use readBox.m
    currPcd = pcd(i);
    [~, data] = readBox(fullfile(folder, pbin.name),fullfile(folder, currPcd.name));
    data = data';

    % reshaping data
    data = reshape(data, ntheta, nr, nx, 5);
    v = data(:,:,:,2); w = data(:,:,:,3);   % axial experiences no change
    data(:,:,:,2) = v .* cos(THETA) + w .* sin(THETA);  % radial
    data(:,:,:,3) = -v .* sin(THETA) + w .* cos(THETA); % azimuthal


    % step 2: calculate fluctuation
    fluc = data - meanfield;

    % steps 3 & 4: calculate combined fluctuations and find the mean
    idx = 1;
    for j = 1:nvars
        for k = j:nvars % so as to not repeat calculations that were already made
            fluc2 = fluc(:,:,:,j) .* fluc(:,:,:,k);
            stress(:,:,:,idx) = stress(:,:,:,idx) + (fluc2 ./ nfiles);
            idx = idx + 1;
        end
    end
end
toc

% save .mat file of reynolds stress
tic
disp('all done!! saving reynolds stress...')
filename = append('reynolds_stress_',test_name);
dirname = append('matrices_',test_name);
out_dir = fullfile('..',dirname,'stress');
if ~exist(out_dir,'dir')
    mkdir(out_dir);
end
save(fullfile(out_dir,filename),'stress','-v7.3');
toc
disp('finished! (*¯︶¯*)')
end