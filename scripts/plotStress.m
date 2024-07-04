function plotStress(folder)
% input:
% folder: path to folder with .mat files containing reynolds stress data from
% post-processing chunks

% output:
% .png and .fig files showing the stress contour at a central plane and
% line profiles for the stresses with the normalized radius and stuff

% prep data for looping
stress = dir(fullfile(folder,'reynolds_stress*'));
vars = {'U','V','W'};
nvars = length(vars);
test_name = folder(strfind(folder, 'M0'):end);

% cylindrical data coordinate definition
x = linspace(0,30, 751)';
r = linspace(0,6, 151)';
% theta = linspace(0,2*pi,128)';
nx = 751;
% ntheta = 128;
nr = 151;

% plot contours
disp('generating contour plots for each component of the stress tensor...')
stress_plot = permute(stress, [2,3,1,4]);   % now in order [x,r,theta,<ij>]
for i = 1:nvars
    var1 = vars{i};
    figure
    for j = i:nvars % so as to not repeat plots that were already made
        var2 = vars{j};
        subplot(nvars,1,j);
        idx = ((i-1)*3) + j;
        sgtitle('Components of the Reynolds Stress Tensor (No Normalization)')
        % plot contour in direction <i,j>
        contourf(x,r,stress_plot(:,:,1,idx),'edgecolor','none');
        hold on
        contourf(x,r.*-1,stress_plot(:,:,65,idx),'edgecolor','none');
        colorbar;
        axis equal

        % title shenanigans
        xlabel("X/D_e, X-Distance from Nozzle Exit");
        ylabel("Y/D_e, Y-Distance from Nozzle Exit");
        title([var1,var2])
    end

    % save fig
    tic
    disp(['saving contour plot ',var1,var2])
    figName = append('stresscontour_',var1,var2,'_',test_name,'.fig');
    pngName = append('stresscontour_',var1,var2,'_',test_name,'.png');
    out_dir = fullfile('..','figs',test_name);
    if ~exist(out_dir,'dir')
        mkdir(out_dir);
    end
    saveas(gcf,fullfile(out_dir,figName));
    saveas(gcf,fullfile(out_dir,pngName));
    disp(['done! saved as ',figName,' AND .png! ♪(´▽｀)'])
    toc
end

% plot lines, much more involved :C
% normalize radius
disp('generating line plots for stress tensor profiles...')
meanfield = dir(fullfile('..',append('matrices_',test_name), ...
    'mean_data','meanfield*')).vol_data;
meanfield = permute(meanfield, [2,3,1,4]);  % now in order [x,r,theta,var]
centerU = meanfield(:,1,1,1);   % meanfield centerline velocity in the u-direction
halfU = centerU ./ 2;    % half of meanfield centerline velocity
meanU = meanfield(:,:,1,1); % meanfield in the U direction at theta = 0
[~,halfR] = find(abs(meanU - halfU) < 0.001);   % r_1/2 for normalization
normR = zeros(nx,nr);
for i = 1:nx
    normR(i,:) = r ./ halfR(i);
end

% normalize stress tensor
normStress = stress_plot ./ (centerU.^2);

% plot normalized tensor versus normalized radius for each combination of
% components at different x-stations (same from meanplot)
xidx = [1,26,126,251,376,501,626,751];
xStations = [0,1,5,10,15,20,25,30];
for i = 1:length(xStations)
    % pull data to plot
    currStress = normStress(xidx(i),:,1,:);
    currR = normR(xidx(i),:);

    % plot
    figure
    plot(currR, currStress);
    title('Profiles of Reynolds Stresses')
    subtitle(['X/D_e = ', num2str(xStations(i))]);
    xlabel('$r/r_{1/2}$','interpreter','latex');
    ylabel('$\frac{\langle u_i u_j\rangle}{U^2_0}$','interpreter','latex','Rotation',0);
    legend('$\langle u^2\rangle$','$\langle uv\rangle$','$\langle uw\rangle$','$\langle v^2\rangle$', ...
        '$\langle vw\rangle$','$\langle w^2\rangle$','interpreter','latex');

    % save fig
    tic
    disp(['saving stress profile ',num2str(i),' of ',num2str(length(xStations))]);
    figName = append('stressprofile_X',num2str(xStations(i)),'_',test_name,'.fig');
    pngName = append('stressprofile_X',num2str(xStations(i)),'_',test_name,'.png');
    saveas(gcf,fullfile(out_dir,figName));
    saveas(gcf,fullfile(out_dir,pngName));
    disp(['done! saved as ',figName,' AND .png! ♪(´▽｀)'])
    toc
end
end