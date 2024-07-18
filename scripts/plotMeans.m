function plotMeans(folder)
% input:
% folder: path to folder with .mat files containing mean field data from
% post-processing chunks (now including the final, complete mean field!)

% output:
% .png and .fig files showing the central plane and cross section at each
% chunk and the final, complete meanfield

% prep data for looping
mats = dir(fullfile(folder, 'mean_data', 'mean*'));
vars = {'u','v','w','Pressure','Density'};
nvars = length(vars);
chunk_size = 2000;
test_name = folder(strfind(folder, 'M0'):end);
[Uj, pj, rhoj] = normData(test_name);
radiusLabel = "R/D_e, Radial Distance from Nozzle Exit";

% cylindrical data coordinate definition
x = linspace(0,30, 751)';
r = linspace(0,6, 151)';
theta = linspace(0,2*pi,128)';

% prep meshgrid for cross section plotting
[R,Theta] = meshgrid(r,theta);
X = R .* cos(Theta);
Y = R .* sin(Theta);

% loop over all mats in folder
for i = 1:length(mats)
    % load matrix data from saved mean files
    data = load(fullfile(folder,'mean_data', mats(i).name)).vol_data;
    data(:,:,:,1:3) = data(:,:,:,1:3)./Uj;
    data(:,:,:,4) = data(:,:,:,4)./pj;
    data(:,:,:,5) = data(:,:,:,5)./rhoj;

    % set time start for plot
    time_start = chunk_size * (i-1);
    time_end = time_start + chunk_size;
    if i == (length(mats)-1)
        time_end = 'end';
    end

    % loop over nvars = 5
    for j = 1:nvars
        % extract data for var j
        var = vars{j};
        mean_plot = permute(data(:,:,:,j), [2,3,1]);

        % prep title names to avoid more mess
        velTitle = append('Mean Velocity Contour for the ',var, ...
            ' Component of the Velocity in the ',test_name,' case');
        otherTitle = append('Mean ',var,' Countour in the ',test_name,' case');
        subTitle = append('Acoustic Time Steps ',num2str(time_start), '-',num2str(time_end));

        % special title for the final matrix which is always the full mean
        % field
        if strfind(mats(i).name, 'meanfield')
            velTitle = append('Mean Velocity Contour for the ',var, ...
                ' Component of the Velocity over All Timesteps');
            otherTitle = append('Mean ',var,' Contour over All Timesteps');
            subTitle = "";
        end

        % plot at central plane (theta = 0, pi)
        figure
        sgtitle("");
        subplot(2,1,1)
        contourf(x,r,mean_plot(:,:,1),'edgecolor','none');
        hold on
        contourf(x,r.*-1,mean_plot(:,:,67),'edgecolor','none');
        colorbar;
        if var == 'u'; caxis([0 1]); end
        if var == 'u'; caxis([0 1]); end
	      if var == 'v'; caxis([-0.03 0.03]); end
        if var == 'w'; caxis([-0.03 0.03]); end
        if contains(test_name, "M0p8")
		      if var == "Pressure"; caxis([0.99 1.001]); end
		      if var == "Density"; caxis([0.88 1.03]); end
	      elseif contains(test_name, "M0p9")
	      	if var == "Pressure"; caxis([0.98 1.001]); end
	      	if var == "Density"; caxis([0.86 1.001]); end
      	else
        	if var == "Pressure"; caxis([0.99 1.001]); end
        	if var == "Density"; caxis([0.95 1]); end
      	end
        axis equal;

        % title shenanigans for central plane plot
        if j < 4
            title([velTitle, ' at \theta = 0']);
        else
            title([otherTitle, ' at \theta = 0']);
        end
        subtitle(subTitle);
        xlabel("X/D_e, X-Distance from Nozzle Exit");
        ylabel("Y/D_e, Y-Distance from Nozzle Exit");

        % plot at central plane (theta = pi/2, 3pi/2)
        subplot(2,1,2);
        contourf(x,r,mean_plot(:,:,33),'edgecolor','none');
        hold on
        contourf(x,r.*-1,mean_plot(:,:,97),'edgecolor','none');
        colorbar;
        if var == 'u'; caxis([0 1]); end
        if var == 'u'; caxis([0 1]); end
	      if var == 'v'; caxis([-0.03 0.03]); end
        if var == 'w'; caxis([-0.03 0.03]); end
        if contains(test_name, "M0p8")
		      if var == "Pressure"; caxis([0.99 1.001]); end
		      if var == "Density"; caxis([0.88 1.03]); end
	      elseif contains(test_name, "M0p9")
	      	if var == "Pressure"; caxis([0.98 1.001]); end
	      	if var == "Density"; caxis([0.86 1.001]); end
      	else
        	if var == "Pressure"; caxis([0.99 1.001]); end
        	if var == "Density"; caxis([0.95 1]); end
      	end
        axis equal;

        % title shenanigans for central plane plot
        if j < 4
            title([velTitle, ' at \theta = \pi/2']);
        else
            title([otherTitle, ' at \theta = \pi/2']);
        end
        subtitle(subTitle);
        xlabel("X/D_e, X-Distance from Nozzle Exit");
        ylabel("Y/D_e, Y-Distance from Nozzle Exit");
        set(gcf, 'position', [100,100,500,500]);

        % save plots
        disp('saving central plane plot...')
        figNameZ = append('meancontour_',var,'_',test_name,'_',num2str(i),'.fig');
        pngNameZ = append('meancontour_',var,'_',test_name,'_',num2str(i),'.png');
        if i == length(mats)
            figNameZ = append('meancontour_',var,'_',test_name,'_all.fig');
            pngNameZ = append('meancontour_',var,'_',test_name,'_all.png');
        end
        out_dir = fullfile('..','figs',test_name);
        if ~exist(out_dir,'dir')
            mkdir(out_dir);
        end
        saveas(gcf,fullfile(out_dir,figNameZ));
        saveas(gcf,fullfile(out_dir,pngNameZ));
        disp(['done! saved as ',figNameZ,' AND .png! ♪(´▽｀)'])

        % find cross section at different x-stations
        X0 = permute(mean_plot(:,1,:), [1,3,2])';
        X1 = permute(mean_plot(:,26,:), [1,3,2])';
        X5 = permute(mean_plot(:,126,:), [1,3,2])';
        X10 = permute(mean_plot(:,251,:), [1,3,2])';
        X15 = permute(mean_plot(:,376,:), [1,3,2])';
        X20 = permute(mean_plot(:,501,:), [1,3,2])';
        X25 = permute(mean_plot(:,626,:), [1,3,2])';
        X30 = permute(mean_plot(:,751,:), [1,3,2])';
        crossSections = {X0 X1 X5 X10 X15 X20 X25 X30};
        xStations = [0,1,5,10,15,20,25,30];

        % plot at each station
        figure
        subplot(2,2,1);
        if j < 4
            sgtitle(velTitle);
        else
            sgtitle(otherTitle);
        end
        idx = 1;
        for k = 1:length(xStations)
            curr = crossSections{k};

            % split plots into two sections, restart plotting for second
            % half and save first half
            if k == 5
                % save plot
                disp('saving cross section plots...')
                figNameR = append('cross_meancontour_X0-10','_',var,'_',test_name,'_',num2str(i),'.fig');
                pngNameR = append('cross_meancontour_X0-10','_',var,'_',test_name,'_',num2str(i),'.png');
                if i == length(mats)
                    figNameR = append('cross_meancontour_X0-10','_',var,'_',test_name,'_all.fig');
                    pngNameR = append('cross_meancontour_X0-10','_',var,'_',test_name,'_all.png');
                end
                saveas(gcf,fullfile(out_dir,figNameR));
                saveas(gcf,fullfile(out_dir,pngNameR));
                disp(['done! saved as ',figNameR,' AND .png! ♪(´▽｀)'])

                % restart plotting for second half by initializing again
                figure
                subplot(2,2,1);
                idx = 1;
                if j < 4
                    sgtitle(velTitle);
                else
                    sgtitle(otherTitle);
                end
            end
            subplot(2,2,idx);
            idx = idx + 1;
            contourf(X,Y,curr,'edgecolor','none');
            colorbar;
            if var == 'u'; caxis([0 1]); end
            if var == 'u'; caxis([0 1]); end
            if var == 'v'; caxis([-0.03 0.03]); end
            if var == 'w'; caxis([-0.03 0.03]); end
            if contains(test_name, "M0p8")
              if var == "Pressure"; caxis([0.99 1.001]); end
              if var == "Density"; caxis([0.88 1.03]); end
            elseif contains(test_name, "M0p9")
              if var == "Pressure"; caxis([0.98 1.001]); end
              if var == "Density"; caxis([0.86 1.001]); end
            else
              if var == "Pressure"; caxis([0.99 1.001]); end
              if var == "Density"; caxis([0.95 1]); end
            end
            if k < 4
                ax = gca;
                ax.XLim = [-1 1];
                ax.YLim = [-1 1];
            elseif k == 4
                ax = gca;
                ax.XLim = [-2 2];
                ax.YLim = [-2 2];
            elseif k == 5
                ax = gca;
                ax.XLim = [-3 3];
                ax.YLim = [-3 3];
            elseif k == 6
                ax = gca;
                ax.XLim = [-4 4];
                ax.YLim = [-4 4];
            end

            % title shenanigans for cross section plots
            title(['X/D_e = ',num2str(xStations(k))]);
            subtitle(subTitle);
            xlabel(radiusLabel);
            ylabel("R/D_e");
        end

        % save plot at the very end!
        disp('saving cross section plots...')
        figNameR = append('cross_meancontour_X15-30','_',var,'_',test_name,'_',num2str(i),'.fig');
        pngNameR = append('cross_meancontour_X15-30','_',var,'_',test_name,'_',num2str(i),'.png');
        if i == length(mats)
            figNameR = append('cross_meancontour_X15-30','_',var,'_',test_name,'.fig');
            pngNameR = append('cross_meancontour_X15-30','_',var,'_',test_name,'.png');
        end
        saveas(gcf,fullfile(out_dir,figNameR));
        saveas(gcf,fullfile(out_dir,pngNameR));
        disp(['done! saved as ',figNameR,' AND .png! ♪(´▽｀)'])
    end
end
