function azimuthalFFT(folder)
% input:
% folder:   directory to the folder containing pbin and pcd files
% data

% output:
% unsure what i need to output here but probably matrices with fft data and
% plots of the fourier decomposition in the azimuthal direction (modes 0-2)

% cylindrical data coordinate definition
nx = 751;
ntheta = 128;
nr = 151;
vars = {'U','V','W','Pressure','Density'};
nvars = length(vars);
x = linspace(0,30, nx)';
r = linspace(0,6, nr)';
theta = linspace(0,2*pi,ntheta)';
m = (-ntheta/2:1:ntheta/2-1)';
[R,THETA] = meshgrid(r,theta);
X = R .* cos(THETA);
Y = R .* sin(THETA);

% pull first datapoint from dataset with readBox to make fluctuation
% calculations
pbin = dir(fullfile(folder,'*.pbin')).name;
pcd = dir(fullfile(folder,'*.pcd'));
pcd = pcd(1).name;
[~, data] = readBox(fullfile(folder,pbin), fullfile(folder,pcd));
data = data';
data = reshape(data,ntheta,nr,nx,nvars);    % reshape data
data = permute(data, [3,2,1,4]);

% pull mean data over all timesteps for given test
test_name = folder(strfind(folder, 'pNozzle'):end);
test_name = test_name(9:12);
mean_data = load(fullfile('..',append('matrices_',test_name),'mean_data', ...
    append('meanfield_',test_name))).vol_data;
%mean_data = load(fullfile('..','vm_mats',append('meanfield_',test_name))).vol_data;
mean_data = permute(mean_data, [3,2,1,4]);

% apply Uj normalization
[Uj, pj, rhoj] = normData(test_name);
data(:,:,:,1:3) = data(:,:,:,1:3)./Uj;
data(:,:,:,4) = data(:,:,:,4)./pj;
data(:,:,:,5) = data(:,:,:,5)./rhoj;
mean_data(:,:,:,1:3) = mean_data(:,:,:,1:3)./Uj;
mean_data(:,:,:,4) = mean_data(:,:,:,4)./pj;
mean_data(:,:,:,5) = mean_data(:,:,:,5)./rhoj;

% calculate fluctuation averaged over theta
theta_mean_data = mean(mean_data, 3);
fluc = data - theta_mean_data;

% compute fourier transform
fluc_hat = fft(fluc, ntheta, 3);
fluc_hat = fftshift(fluc_hat, 3);

% find fourier coefficients
fluc_0 = permute(fluc_hat(:,:,ntheta/2+1,:), [1,2,4,3]);
fluc_a1 = permute(real(2.*fluc_hat(:,:,ntheta/2+2,:)), [1,2,4,3]);
fluc_b1 = permute(imag(2.*fluc_hat(:,:,ntheta/2+2,:)), [1,2,4,3]);
fluc_a2 = permute(real(2.*fluc_hat(:,:,ntheta/2+3,:)), [1,2,4,3]);
fluc_b2 = permute(imag(2.*fluc_hat(:,:,ntheta/2+3,:)), [1,2,4,3]);
modes = {fluc_0, fluc_a1, fluc_b1, fluc_a2, fluc_b2};
modenames = {'m = 0', 'm = 1a', 'm = 1b', 'm = 2a', 'm = 2b'};

% plot fluctuation for all vars
for i = 1:nvars
    currVar = vars{i};
    figure
    sgtitle([currVar, ' fluctuations in the ', test_name, ' case'])
    subplot(6,1,1);
    contourf(x,r,fluc(:,:,1,i)','edgecolor','none');
    colorbar
    ylabel('$y/D_e$','interpreter','latex');
    title('2D plane')
    ax = gca; ax.XLim = [0 10]; ax.YLim = [0 2];
    set(gcf, 'position', [100,100,1000,750]);
    for j = 1:length(modes)
        currMode = modes{j};
        subplot(6,1,j+1);
        contourf(x,r,real(currMode(:,:,i))','edgecolor','none');
        colorbar
        ylabel('$y/D_e$','interpreter','latex');
        title(modenames{j});
        ax = gca; ax.XLim = [0 10]; ax.YLim = [0 2];
        set(gcf, 'position', [100,100,1000,750]);
        if j == length(modes)
            xlabel('$x/D_e$','interpreter','latex');
        end
    end

    % save plot
    out_dir = fullfile('..','figs',test_name);
    figName = append('azimuthalFFT_fluc_',vars{i},'_',test_name,'.fig');
    pngName = append('azimuthalFFT_fluc_',vars{i},'_',test_name,'.png');
    saveas(gcf,fullfile(out_dir,figName));
    saveas(gcf,fullfile(out_dir,pngName));
end

% load stress tensor, already in shape [nx,nr,ntheata,nvars]
test_name = folder(strfind(folder, 'M0'):strfind(folder, 'M0')+3);
stress = load(fullfile('..',append('matrices_',test_name), ...
    'stress',append('reynolds_stress_',test_name))).stress;
stress = stress ./ (Uj^2);
stress = permute(stress, [3,2,1,4]);
vars = {'UU','UV','UW','VV','VW','WW'};
nvars = length(vars);

% loop over components of the stress tensor
for i = 1:nvars
    % compute fourier transform
    q_hat = fft(stress, ntheta, 3);
    q_hat = fftshift(q_hat, 3);

    % find fourier coefficients
    q_0 = permute(q_hat(:,:,ntheta/2+1,:), [1,2,4,3]);
    q_a1 = permute(real(2.*q_hat(:,:,ntheta/2+2,:)), [1,2,4,3]);
    q_b1 = permute(imag(2.*q_hat(:,:,ntheta/2+2,:)), [1,2,4,3]);
    q_a2 = permute(real(2.*q_hat(:,:,ntheta/2+3,:)), [1,2,4,3]);
    q_b2 = permute(imag(2.*q_hat(:,:,ntheta/2+3,:)), [1,2,4,3]);
    modes = {q_0, q_a1, q_b1, q_a2, q_b2};

    % plot modes
    currVar = vars{i};
    figure
    sgtitle([currVar, ' Component of the Reynolds Stress Tensor in the ', ...
        test_name, ' case'])
    subplot(6,1,1);
    contourf(x,r,stress(:,:,1,i)','edgecolor','none');
    colorbar
    ylabel('$y/D_e$','interpreter','latex');
    title('2D plane')
    ax = gca; ax.XLim = [0 10]; ax.YLim = [0 2];
    set(gcf, 'position', [100,100,1000,750]);
    for j = 1:length(modes)
        currMode = modes{j};
        subplot(6,1,j+1);
        contourf(x,r,real(currMode(:,:,i))','edgecolor','none');
        colorbar
        ylabel('$y/D_e$','interpreter','latex');
        title(modenames{j});
        ax = gca; ax.XLim = [0 10]; ax.YLim = [0 2];
        set(gcf, 'position', [100,100,1000,750]);
        if j == length(modes)
            xlabel('$x/D_e$','interpreter','latex');
        end
    end

    % save plot
    out_dir = fullfile('..','figs',test_name);
    figName = append('azimuthalFFT_stress_',vars{i},'_',test_name,'.fig');
    pngName = append('azimuthalFFT_stress_',vars{i},'_',test_name,'.png');
    saveas(gcf,fullfile(out_dir,figName));
    saveas(gcf,fullfile(out_dir,pngName));
end
end