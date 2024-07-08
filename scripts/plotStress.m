function plotStress(folder)
% input:
% folder: path to folder with .mat files containing reynolds stress data from
% post-processing chunks

% output:
% .png and .fig files showing the stress contour at a central plane and
% line profiles for the stresses with the normalized radius and stuff

% prep data for looping
test_name = folder(strfind(folder, 'M0'):strfind(folder, 'M0')+3);
stress = load(fullfile(folder,append('reynolds_stress_',test_name))).stress;
vars = {'U','V','W'};
nvars = length(vars);

% cylindrical data coordinate definition
x = linspace(0,30, 751)';
r = linspace(0,6, 151)';
% theta = linspace(0,2*pi,128)';
nx = 751;
% ntheta = 128;
nr = 151;

% plot contours
disp('generating contour plots for each component of the stress tensor...')
stress_plot = permute(stress, [3,2,1,4]);   % now in order [x,r,theta,<ij>]
idx = 1;
mat_idx = 1;
fig_idx = 1;
for i = 1:nvars
    var1 = vars{i};
    if idx == 1 || idx == 4
        figure
        sgtitle('Components of the Reynolds Stress Tensor (No Normalization)')
        idx = 1;
    end
    for j = i:nvars % so as to not repeat plots that were already made
        var2 = vars{j};
        subplot(nvars,1,idx);
        % plot contour in direction <i,j>
        contourf(x,r,stress_plot(:,:,1,mat_idx)','edgecolor','none');
        hold on
        contourf(x,r.*-1,stress_plot(:,:,67,mat_idx)','edgecolor','none');
        colorbar;
        % axis equal

        % title shenanigans
        xlabel("X/D_e, X-Distance from Nozzle Exit");
        ylabel("Y/D_e");
        title([var1,var2,' Component'])
        idx = idx + 1;
        mat_idx = mat_idx + 1;
    end

    % save fig
    if idx == 4
        tic
        disp(['saving contour plot ',num2str(fig_idx)])
        figName = append('stresscontour_',num2str(fig_idx),'_',test_name,'.fig');
        pngName = append('stresscontour_',num2str(fig_idx),'_',test_name,'.png');
        out_dir = fullfile('..','figs',test_name);
        if ~exist(out_dir,'dir')
            mkdir(out_dir);
        end
        saveas(gcf,fullfile(out_dir,figName));
        saveas(gcf,fullfile(out_dir,pngName));
        fig_idx = fig_idx + 1;
        disp(['done! saved as ',figName,' AND .png! ♪(´▽｀)'])
        toc
    end
end

% plot lines, much more involved :C
% normalize radius
disp('generating line plots for stress tensor profiles...')
meanfield = load(fullfile('..',append('matrices_',test_name), ...
    'mean_data',append('meanfield_',test_name))).vol_data;
meanfield = permute(meanfield, [3,2,1,4]);  % now in order [x,r,theta,var]
centerU = meanfield(:,1,1,1);   % meanfield centerline velocity in the u-direction
halfU = centerU ./ 2;    % half of meanfield centerline velocity
meanU = meanfield(:,:,1,1); % meanfield in the U direction at theta = 0
halfR = zeros(nx,1); halfU_fit = zeros(nx,1); normR = zeros(nx,nr);
for i = 1:nx
    if i <= 70; type = 'gauss2';
    else; type = 'poly3'; end
    f = fit(r,meanU(i,:)',type);
    fun = @(r) halfU(i) - f(r);
    halfR(i) = fzero(fun,[0 5]);
    halfU_fit(i) = f(halfR(i));
    normR(i,:) = r ./ halfR(i);
end

% plot jet half width and half velocity for debugging purposes
figure
fig1 = plot(x,halfR,'o');
title("Jet Half-Width versus X-Distance");
xlabel("X/D_e");
ylabel("r_{1/2}",'rotation',0);

figure
fig2 = plot(x,halfU,'color','magenta','linewidth',2);
hold on
plot(x,halfU_fit,'--','color','black','linewidth',2)
plot(x,centerU,'color','blue')
title("Jet Velocity Fit Test")
xlabel("X/D_e")
ylabel("Jet Velocity")
legend("true U_{1/2}","fitted U_{1/2}","Centerline U",'color','blue','linewidth',2)

% save plots
figName1 = "fit_test_r.fig";
pngName1 = "fit_test_r.png";
figName2 = "fit_test_u.fig";
pngName2 = "fit_test_u.png";
saveas(fig1,fullfile(out_dir,figName1));
saveas(fig1,fullfile(out_dir,pngName1));
saveas(fig2,fullfile(out_dir,figName2));
saveas(fig2,fullfile(out_dir,pngName2));

% normalize stress tensor
normStress = stress_plot ./ (centerU.^2);

% plot normalized tensor versus normalized radius for each combination of
% components at different x-stations (same from meanplot)
xidx = [1,26,126,251,376,501,626,751];
xStations = [0,1,5,10,15,20,25,30];
for i = 1:length(xStations)
    % pull data to plot
    currStress = permute(normStress(xidx(i),:,1,:), [2,4,1,3]);
    currR = normR(xidx(i),:)';

    % plot
    figure
    plot(currR, currStress);
    title('Profiles of Reynolds Stresses')
    subtitle(['X/D_e = ', num2str(xStations(i))]);
    xlabel('$r/r_{1/2}$','interpreter','latex');
    ylabel('$\frac{\langle u_i u_j\rangle}{U^2_0}$','interpreter','latex','Rotation',0);
    legend('$\langle uu\rangle$','$\langle uv\rangle$','$\langle uw\rangle$','$\langle vv\rangle$', ...
        '$\langle vw\rangle$','$\langle ww\rangle$','interpreter','latex');

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