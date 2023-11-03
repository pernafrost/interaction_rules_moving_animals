
% PLOT_TURNING_ACCELERATION Plots the rules of interaction as implemented
% in the simulation
%
%   [] = plot_turning_acceleration(particleType, typeNumber, fps)
%
% INPUT
%   particleType: contains all the interaction parameters of different
%   classes of particles
%   fps: frames per second
%
% OUTPUT
%   None
%
% DESCRIPTION
% The function plots the interaction rules between particles as implemented
% in the simulation
%
% NOTE
% Please notice that the schooling interactions are defined inside this
% function; check that these are the same functions implemented in
% schooling_interactions.m

% Copyright: Andrea Perna, perna@math.uu.se
% Department of Mathematics, Uppsala University
% Box 480, 75106 Uppsala, Sweden
% http://perna.fr

% Latest update 2013-03-05

function plot_turning_acceleration(particleType, fps, varargin)

if (size(varargin, 2) >= 1) && (isfinite(varargin{1}))% && (~isempty(varargin{1})
    defaultPhi = varargin{1};
else
    defaultPhi = 0;
end

if (size(varargin, 2) >= 2) && (isfinite(varargin{2}))% && (~isempty(varargin{2})
    defaultDist = varargin{2};
else
    defaultDist = 50;
end

if (size(varargin, 2) >= 3) && (isfinite(varargin{3}))% && (~isempty(varargin{3})
    defaultSpeed = varargin{3};
else
    % here the default speed of the neighbour is equal to the default speed of the particle
    defaultSpeed = particleType.preferredSpeed(1);
end

if (size(varargin, 2) >= 4) && (isfinite(varargin{4}))% && (~isempty(varargin{3})
    plotInXY = varargin{4};
else
    % decide whether we want the plot in deltaX, deltaY or in rho, theta
    plotInXY = 0; % by default plot in rho, theta
end


% if 0 it only plots some phi bins
useAllPhiBins = 0;


% font size for the figures
useFontSize = 18;


% create the colormap
x=(0:1:255)/255;

R = x/0.32 - 0.78125;
R(R<0) = 0;
R(R>1) = 1;

G = 2*x - 0.84;
G(G<0) = 0;
G(G>1) = 1;


B(find(x<=0.25)) = 4*x(x<=0.25);
B(intersect(find(x>0.25), find(x<=0.42))) = 1;
B(intersect(find(x>0.42), find(x<=0.92))) = -2*x(x>0.42 & x<=0.92) + 1.84;
B(find(x>0.92)) = x(x>0.92)/0.08 - 11.5;

myColourMap=[R', G', B'];



%% create two particles that will be used to simulate interactions
focalParticle.x = 0;
focalParticle.y = 0;
focalParticle.phi = 0;
focalParticle.speed = particleType.preferredSpeed(1);
focalParticle.species = 1;

neighbour.distance = defaultDist;
neighbour.phi = defaultPhi;
neighbour.speed = defaultSpeed;
neighbour.species = 1; % the neighbour is of the same species right now



if (plotInXY == 1)
    
    
    %% bins for the figure
    deltaXBins = -10:0.5:10;
    deltaYBins = -10:0.5:10;
    centerDeltaXBins = deltaXBins(1:end-1) + diff(deltaXBins)/2;
    centerDeltaYBins = deltaYBins(1:end-1) + diff(deltaYBins)/2;
    
    thetaBins = linspace(-pi, pi, 25); % radians (13 bins make 30 degrees each)
    centerThetaBins = thetaBins(1:end-1) + diff(thetaBins)/2;
    
    if useAllPhiBins
        phiBins = linspace(-pi, pi, 25); % radians (13 bins make 30 degrees each)
        centerPhiBins = phiBins(1:end-1) + diff(phiBins)/2;
    else
        phiBins = linspace(-pi/20, pi/20, 11);
        centerPhiBins = phiBins(1:end-1) + diff(phiBins)/2;
    end
    
    
    % initialize empty arrays
    turningDeltaXDeltaY = NaN(length(deltaXBins), length(deltaYBins));
    accelerationDeltaXDeltaY = NaN(length(deltaXBins), length(deltaYBins));
    
    turningPhiTheta = NaN(length(phiBins), length(thetaBins));
    accelerationPhiTheta = NaN(length(phiBins), length(thetaBins));
    
    
    
    for countTheta=1:length(thetaBins)-1
        countTheta
        neighbour.theta = centerThetaBins(countTheta);
        
        for countPhi=1:length(phiBins)-1
            
            neighbour.distance = defaultDist;
            neighbour.phi = centerPhiBins(countPhi);
            
            
            P = nan(2,5);
            P(1,:) = [focalParticle.x, focalParticle.y, focalParticle.phi, focalParticle.speed, focalParticle.species];
            P(2,:) = [neighbour.distance * cos(neighbour.theta), neighbour.distance * sin(neighbour.theta), neighbour.phi, neighbour.speed, neighbour.species];
            
            Pnew = schooling_interactions(P, particleType, fps);
            clear schooling_interactions; % I have a persistent variable inside the function that I want to clear here
            [turning, Vtrue] = cart2pol(Pnew(1,1), Pnew(1,2));
            turningPhiTheta(countPhi, countTheta) = turning;
            accelerationPhiTheta(countPhi, countTheta) = (Vtrue - focalParticle.speed)/fps;
            
        end
    end
    
    
    
    for countDeltaX = 1:length(deltaXBins)-1
        countDeltaX
        for countDeltaY = 1:length(deltaYBins)-1
            
            % this transformation from cart to pol coordinates is not
            % useful, but I keep it here to maintain the same structure as
            % in the other 2D histograms
            [neighbour.theta, neighbour.distance] = cart2pol(centerDeltaXBins(countDeltaX), centerDeltaYBins(countDeltaY));
            neighbour.phi = defaultPhi;
            
            
            P = nan(2,5);
            P(1,:) = [focalParticle.x, focalParticle.y, focalParticle.phi, focalParticle.speed, focalParticle.species];
            P(2,:) = [neighbour.distance * cos(neighbour.theta), neighbour.distance * sin(neighbour.theta), neighbour.phi, neighbour.speed, neighbour.species];
            
            Pnew = schooling_interactions(P, particleType, fps);
            clear schooling_interactions; % I have a persistent variable inside the function that I want to clear here
            [turning, Vtrue] = cart2pol(Pnew(1,1), Pnew(1,2));
            
            turningDeltaXDeltaY(countDeltaX, countDeltaY) = turning;
            accelerationDeltaXDeltaY(countDeltaX, countDeltaY) = (Vtrue - focalParticle.speed)/fps;
        end
    end
    
    
    
    
    
    %% figure for turning angle
    valueForColourScale = (max(max(abs(turningDeltaXDeltaY))));
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(deltaXBins,deltaYBins,turningDeltaXDeltaY); caxis([-valueForColourScale, valueForColourScale]);
    hold on;
    fill([0, -0.25, 0, 0.25],[0.5, -0.4, 0, -0.4], 'w', 'LineWidth', 1, 'LineStyle', 'none');
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel('left-right distance to neighbour (cm)'); ylabel('front-back distance to neighbour (cm)');
    title(sprintf('Neighbour orientation phi=%0.1f', defaultPhi), 'FontName', 'Arial');
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'turning angle', 'FontSize', useFontSize);
    colormap(myColourMap);
    
    figureFileName = 'implemented_turning_angle_vs_deltaX_deltaY.eps';
    % print(gcf, '-depsc2', figureFileName);
    
    
    
    
    
    
    
    
    
    %% figure for acceleration
    valueForColourScale = max(max(abs(accelerationDeltaXDeltaY)));
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(deltaXBins,deltaYBins,accelerationDeltaXDeltaY); caxis([-valueForColourScale, valueForColourScale]);
    hold on;
    fill([0, -0.25, 0, 0.25],[0.5, -0.4, 0, -0.4], 'w', 'LineWidth', 1, 'LineStyle', 'none');
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
         'FontName', 'Arial');
    
    
    xlabel('left-right distance to neighbour (cm)'); ylabel('front-back distance to neighbour (cm)');
    title(sprintf('Neighbour orientation phi=%0.1f', defaultPhi), 'FontName', 'Arial');
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'acceleration (cm/s^2)', 'FontSize', useFontSize);
    colormap(myColourMap);
    
    figureFileName = 'implemented_acceleration_vs_R_theta.eps';
    % print(gcf, '-depsc2', figureFileName);
    
    
    
    
    
    
    
    
    
    
    
    
    % figure for turning angle phi theta
    valueForColourScale = (max(max(abs(turningPhiTheta))));
    figure,
    set(gcf, 'Position', [1 1 1000 640]);
    imagesc(turningPhiTheta(1:end-1, 1:end-1) , [-valueForColourScale, valueForColourScale]);
    
    
    if useAllPhiBins
        set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
            0.5:4:length(phiBins), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'YTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'FontName', 'symbol');
    else
        set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
            linspace(0.5, length(phiBins)-0.5, 7), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'YTickLabel', {'-p/20', '-p/30', '-p/60', '0', 'p/60', 'p/30', 'p/20'}, 'FontName', 'symbol');
    end
    
    
    xlabel('\theta (rad)', 'FontName', 'arial'); ylabel('\phi (rad)', 'FontName', 'arial');
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'turning (rad / time step)', 'FontSize', useFontSize);
    axis xy equal tight;
    text(2.5, length(thetaBins) - 2.5, 'B', 'FontName', 'Arial', 'FontSize', 56);
    
    title(sprintf('Neighbour distance r=%0.1f', defaultDist), 'FontName', 'Arial');
    colormap(myColourMap);
    
    figureFileName = 'implemented_turning_angle_vs_phi_theta.eps';
    % print(gcf, '-depsc2', figureFileName);
    
    
    
    
    
    
    
    
    
    
    % figure for acceleration phi theta
    valueForColourScale = (max(max(abs(accelerationPhiTheta))));
    figure,
    set(gcf, 'Position', [1 1 1000 640]);
    imagesc(accelerationPhiTheta(1:end-1, 1:end-1) , [-valueForColourScale, valueForColourScale]);
    
    if useAllPhiBins
        set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
            0.5:4:length(phiBins), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'YTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'FontName', 'symbol');
    else
        set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
            linspace(0.5, length(phiBins)-0.5, 7), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'YTickLabel', {'-p/20', '-p/30', '-p/60', '0', 'p/60', 'p/30', 'p/20'}, 'FontName', 'symbol');
    end
    
    xlabel('\theta (rad)', 'FontName', 'arial'); ylabel('\phi (rad)', 'FontName', 'arial');
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'acceleration (cm/s^2)', 'FontSize', useFontSize);
    axis xy equal tight;
    title(sprintf('Neighbour distance r=%0.1f', defaultDist), 'FontName', 'Arial');
    colormap(myColourMap);
    
    figureFileName = 'implemented_acceleration_vs_phi_theta.eps';
    % print(gcf, '-depsc2', figureFileName);
    
    
    
else
    %% bins in distance and theta
    
    
    % bins for the figure
    distanceBins = 0:0.5:10;
    thetaBins = linspace(-pi, pi, 25); % radians (13 bins make 30 degrees each)
    centerThetaBins = thetaBins(1:end-1) + diff(thetaBins)/2;
    centerDistanceBins = distanceBins(1:end-1) + diff(distanceBins)/2;
    
    
    if useAllPhiBins
        phiBins = linspace(-pi, pi, 25); % radians (13 bins make 30 degrees each)
        centerPhiBins = phiBins(1:end-1) + diff(phiBins)/2;
    else
        phiBins = linspace(-pi/20, pi/20, 11);
        centerPhiBins = phiBins(1:end-1) + diff(phiBins)/2;
    end
    
    
    % initialize empty arrays
    turningDTheta = NaN(length(distanceBins), length(thetaBins));
    accelerationDTheta = NaN(length(distanceBins), length(thetaBins));
    
    turningPhiTheta = NaN(length(phiBins), length(thetaBins));
    accelerationPhiTheta = NaN(length(phiBins), length(thetaBins));
    
    
    
    for countTheta=1:length(thetaBins)-1
        countTheta
        neighbour.theta = centerThetaBins(countTheta);
        
        
        for countD=1:length(distanceBins)-1
            
            
            neighbour.distance = centerDistanceBins(countD);
            neighbour.phi = defaultPhi;
            
            
            P = nan(2,5);
            P(1,:) = [focalParticle.x, focalParticle.y, focalParticle.phi, focalParticle.speed, focalParticle.species];
            P(2,:) = [neighbour.distance * cos(neighbour.theta), neighbour.distance * sin(neighbour.theta), neighbour.phi, neighbour.speed, neighbour.species];
            
            Pnew = schooling_interactions(P, particleType, fps);
            clear schooling_interactions; % I have a persistent variable inside the function that I want to clear here
            [turning, Vtrue] = cart2pol(Pnew(1,1), Pnew(1,2));
            
            turningDTheta(countD, countTheta) = turning;
            accelerationDTheta(countD, countTheta) = (Vtrue - focalParticle.speed)/fps;
        end
        
        for countPhi=1:length(phiBins)-1
            
            neighbour.distance = defaultDist;
            neighbour.phi = centerPhiBins(countPhi);
            
            
            P = nan(2,5);
            P(1,:) = [focalParticle.x, focalParticle.y, focalParticle.phi, focalParticle.speed, focalParticle.species];
            P(2,:) = [neighbour.distance * cos(neighbour.theta), neighbour.distance * sin(neighbour.theta), neighbour.phi, neighbour.speed, neighbour.species];
            
            Pnew = schooling_interactions(P, particleType, fps);
            clear schooling_interactions; % I have a persistent variable inside the function that I want to clear here
            [turning, Vtrue] = cart2pol(Pnew(1,1), Pnew(1,2));
            turningPhiTheta(countPhi, countTheta) = turning;
            accelerationPhiTheta(countPhi, countTheta) = (Vtrue - focalParticle.speed)/fps;
            
        end
    end
    
    
    X = distanceBins'*cos(pi/2 -thetaBins);
    Y = distanceBins'*sin(pi/2 -thetaBins);
    
    
    
    
    
    
    
    
    
    %% figure for turning angle
    valueForColourScale = (max(max(abs(turningDTheta))));
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(X,Y,turningDTheta); caxis([-valueForColourScale, valueForColourScale]);
    hold on;
    fill([0, -0.25, 0, 0.25],[0.5, -0.4, 0, -0.4], 'w', 'LineWidth', 1, 'LineStyle', 'none');
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'XLim', [-max(distanceBins)*12/10 max(distanceBins)*12/10], 'YLim', [-max(distanceBins)*12/10, max(distanceBins)*12/10], ... 'XTickLabel', -max(distanceBins):1:max(distanceBins), 'YTickLabel',  -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    text(0,max(distanceBins)*11/10,'\vartheta = 0', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
    text(max(distanceBins)*11/10,0,'\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
    text(-max(distanceBins)*11/10,0,'-\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
    text(0,-max(distanceBins)*11/10,'\pi', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
    
    % text(-max(distanceBins)*11/10, max(distanceBins)*11/10, 'A', 'FontName', 'Arial', 'FontSize', 50);
    %
    % if defaultPhi == 0
    % text(-max(distanceBins)*11/10, +max(distanceBins)*11/10, {['\phi=0\pi']}, 'FontName', 'Arial', 'FontSize', useFontSize, 'Color', 'k');
    % else
    % text(-max(distanceBins)*11/10, +max(distanceBins)*11/10, {['\phi =', num2str(rats(defaultPhi/pi)), '\pi']}, 'FontName', 'Arial', 'FontSize', useFontSize, 'Color', 'k');
    % end
    
    
    xlabel('distance from neighbour r (cm)'); ylabel('distance from neighbour r (cm)');
    title(sprintf('Neighbour orientation phi=%0.1f', defaultPhi), 'FontName', 'Arial');
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'turning angle', 'FontSize', useFontSize);
    colormap(myColourMap);
    
    figureFileName = 'implemented_turning_angle_vs_R_theta.eps';
    % print(gcf, '-depsc2', figureFileName);
    
    
    
    
    
    
    
    
    
    %% figure for acceleration
    valueForColourScale = max(max(abs(accelerationDTheta)));
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(X,Y,accelerationDTheta); caxis([-valueForColourScale, valueForColourScale]);
    hold on;
    fill([0, -0.25, 0, 0.25],[0.5, -0.4, 0, -0.4], 'w', 'LineWidth', 1, 'LineStyle', 'none');
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'XLim', [-max(distanceBins)*12/10 max(distanceBins)*12/10], 'YLim', [-max(distanceBins)*12/10, max(distanceBins)*12/10], ... 'XTickLabel', -max(distanceBins):1:max(distanceBins), 'YTickLabel',  -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    text(0,max(distanceBins)*11/10,'\vartheta = 0', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
    text(max(distanceBins)*11/10,0,'\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
    text(-max(distanceBins)*11/10,0,'-\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
    text(0,-max(distanceBins)*11/10,'\pi', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
    
    % text(-max(distanceBins)*11/10, max(distanceBins)*11/10, 'C', 'FontName', 'Arial', 'FontSize', 50);
    
    %
    % if defaultPhi == 0
    % text(-max(distanceBins)*11/10, +max(distanceBins)*11/10, {'\phi=0\pi'}, 'FontName', 'Arial', 'FontSize', useFontSize, 'Color', 'k');
    % else
    % text(-max(distanceBins)*11/10, +max(distanceBins)*11/10, {'\phi =', num2str(rats(defaultPhi/pi)), '\pi'}, 'FontName', 'Arial', 'FontSize', useFontSize, 'Color', 'k');
    % end
    
    
    xlabel('distance from neighbour r (cm)'); ylabel('distance from neighbour r (cm)');
    title(sprintf('Neighbour orientation phi=%0.1f', defaultPhi), 'FontName', 'Arial');
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'acceleration (cm/s^2)', 'FontSize', useFontSize);
    colormap(myColourMap);
    
    figureFileName = 'implemented_acceleration_vs_R_theta.eps';
    % print(gcf, '-depsc2', figureFileName);
    
    
    
    
    
    
    
    
    
    
    
    
    % figure for turning angle phi theta
    valueForColourScale = (max(max(abs(turningPhiTheta))));
    figure,
    set(gcf, 'Position', [1 1 1000 640]);
    imagesc(turningPhiTheta(1:end-1, 1:end-1) , [-valueForColourScale, valueForColourScale]);
    
    
    if useAllPhiBins
        set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
            0.5:4:length(phiBins), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'YTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'FontName', 'symbol');
    else
        set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
            linspace(0.5, length(phiBins)-0.5, 7), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'YTickLabel', {'-p/20', '-p/30', '-p/60', '0', 'p/60', 'p/30', 'p/20'}, 'FontName', 'symbol');
    end
    
    
    xlabel('\theta (rad)', 'FontName', 'arial'); ylabel('\phi (rad)', 'FontName', 'arial');
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'turning (rad / time step)', 'FontSize', useFontSize);
    axis xy equal tight;
    text(2.5, length(thetaBins) - 2.5, 'B', 'FontName', 'Arial', 'FontSize', 56);
    
    title(sprintf('Neighbour distance r=%0.1f', defaultDist), 'FontName', 'Arial');
    colormap(myColourMap);
    
    figureFileName = 'implemented_turning_angle_vs_phi_theta.eps';
    % print(gcf, '-depsc2', figureFileName);
    
    
    
    
    
    
    
    
    
    
    % figure for acceleration phi theta
    valueForColourScale = (max(max(abs(accelerationPhiTheta))));
    figure,
    set(gcf, 'Position', [1 1 1000 640]);
    imagesc(accelerationPhiTheta(1:end-1, 1:end-1) , [-valueForColourScale, valueForColourScale]);
    
    if useAllPhiBins
        set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
            0.5:4:length(phiBins), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'YTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'FontName', 'symbol');
    else
        set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:4:length(thetaBins), 'YTick', ...
            linspace(0.5, length(phiBins)-0.5, 7), 'TickDir', 'out', 'XTickLabel', {'-p', '-2/3p', '-p/3', '0', 'p/3', '2/3p', 'p'}, ...
            'YTickLabel', {'-p/20', '-p/30', '-p/60', '0', 'p/60', 'p/30', 'p/20'}, 'FontName', 'symbol');
    end
    
    xlabel('\theta (rad)', 'FontName', 'arial'); ylabel('\phi (rad)', 'FontName', 'arial');
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'acceleration (cm/s^2)', 'FontSize', useFontSize);
    axis xy equal tight;
    title(sprintf('Neighbour distance r=%0.1f', defaultDist), 'FontName', 'Arial');
    colormap(myColourMap);
    
    figureFileName = 'implemented_acceleration_vs_phi_theta.eps';
    % print(gcf, '-depsc2', figureFileName);
    
end


    
    return
    
    
    
    
    
    
    
    
    
    
    
    
    
