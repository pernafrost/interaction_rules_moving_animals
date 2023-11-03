

function [] = make_colourful_plots(sf, sn, r, theta, phi, af, turnf, framesPerTimeUnit, timeUnit, spaceUnit, useAllPhiBins, isSavingFigures)
% function make_colourful_plots
% This function plots heatmaps with the responses of a focal individual to
% the position and movement of its neighbours
% Somce care should be taken with the angles, which are assumed to increase
% counterclockwise, but are plotted here in a compass-like style, with zero
% radians on top and angles increasing clockwise.
% Cartesian axes are assumed to be increasing from left to right and from 
% bottom to top (in images the y axis is often considered to start from the
% top, which sometimes produces flipped coordinates when tracking animals 
% from a video)
% 
% INPUT:
% 
% sf: speed of the focal individual (vertical vector)
% sn: speed of the neighbour (vertical vector, or matrix if there are multiple neighbours)
% r: distance to the neighbour (or to the neighbours)
% theta: direction of the neighbour in the frame of reference of the focal
% individual
% phi: relative orientation of the neighobur (or neighbours) with respect
% to the focal individual
% af: acceleration of the focal individual
% turnf: turning of the focal individual
% framesPerTimeUnit: how many frames in a time unit (e.g. frames per 
% second)
% timeUnit: selected unit for time (e.g. 'seconds', 'frames', 'minutes')
% spaceUnit: selected unit for space (e.g. 'cm', 'm', 'pixels')
% useAllPhiBins: if true it should use phi bins from -pi to pi; if false
% it focuses only on the situation when the neighbour is nearly aligned
% isSavingFigures: if 1 than figures are also saved as eps files 
%
% Written by:
% Andrea Perna
% http://www.perna.fr
%
% Date:
% 2014 / 04 / 18


% Parameters of the figures
useFontSize = 18; % font size for the figures

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




%% bins for the figure
deltaXBins = -10:0.5:10;
deltaYBins = -10:0.5:10;
centerDeltaXBins = deltaXBins(1:end-1) + diff(deltaXBins)/2;
centerDeltaYBins = deltaYBins(1:end-1) + diff(deltaYBins)/2;

thetaBins = linspace(-pi, pi, 25); % radians (13 bins make 30 degrees each)
centerThetaBins = thetaBins(1:end-1) + diff(thetaBins)/2;
thetaBins(end) = pi+0.0001; % I add this because histc is of the form EDGES(k) <= X(i) < EDGES(k+1)


if useAllPhiBins
    phiBins = linspace(-pi, pi, 25); % radians (13 bins make 30 degrees each)
    centerPhiBins = phiBins(1:end-1) + diff(phiBins)/2;
    phiBins(end) = pi+0.0001; % I add this because histc is of the form EDGES(k) <= X(i) < EDGES(k+1)
else
    phiBins = linspace(-pi/20, pi/20, 11);
    centerPhiBins = phiBins(1:end-1) + diff(phiBins)/2;
end

rBins = 0:0.5:10; % space units
centerDistanceBins = rBins(1:end-1) + diff(rBins)/2;

largeRBins = [0, 5, inf]; % or [0, 4, 6, inf]; % assuming that the repulsion radius is 5 we bin separately
% events when the neighbour is in the repulsion or in the attraction zone



% Left-right movement depending on the position of the neighbour





[nInDeltaX, valInDeltaX] = histc(r .* cos(theta), deltaXBins);
[nInDeltaY, valInDeltaY] = histc(r .* sin(theta), deltaYBins);

[nInR, valInR] = histc(r, rBins);
[nInLargeR, valInLargeR] = histc(r, largeRBins);
[nInTheta, valInTheta] = histc(theta, thetaBins);
[nInPhi, valInPhi] = histc(phi, phiBins);

clear meanYStepDXDY meanXStepDXDY meanTurnDXDY circSTDDXDY meanAccDXDY nCountsDXDY;
clear meanYStepRTheta meanXStepRTheta meanTurnRTheta circSTDRTheta meanAccRTheta nCountsRTheta;
clear meanYStepPhiTheta meanXStepPhiTheta meanTurnPhiTheta circSTDPhiTheta meanAccPhiTheta nCountsPhiTheta;

% initialize empty arrays
meanTurnDXDY = NaN(length(deltaXBins), length(deltaYBins)); % x is front back, y is left right
meanAccDXDY = NaN(length(deltaXBins), length(deltaYBins));
circSTDDXDY = NaN(length(deltaXBins), length(deltaYBins));
nCountsDXDY = NaN(length(deltaXBins), length(deltaYBins));
meanYStepDXDY = NaN(length(deltaXBins), length(deltaYBins));
meanXStepDXDY = NaN(length(deltaXBins), length(deltaYBins));

meanTurnRTheta = NaN(length(rBins), length(thetaBins));
meanAccRTheta = NaN(length(rBins), length(thetaBins));
circSTDRTheta = NaN(length(rBins), length(thetaBins));
nCountsRTheta = NaN(length(rBins), length(thetaBins));
meanYStepRTheta = NaN(length(rBins), length(thetaBins));
meanXStepRTheta = NaN(length(rBins), length(thetaBins));

meanTurnPhiTheta = NaN(length(phiBins), length(thetaBins));
meanAccPhiTheta = NaN(length(phiBins), length(thetaBins));
circSTDPhiTheta = NaN(length(phiBins), length(thetaBins));
nCountsPhiTheta = NaN(length(phiBins), length(thetaBins));
meanYStepPhiTheta = NaN(length(phiBins), length(thetaBins));
meanXStepPhiTheta = NaN(length(phiBins), length(thetaBins));

% first construct the maps in DeltaX and DeltaY
for cX=1:length(deltaXBins) -1
    cX
    for cY=1:length(deltaYBins) -1
        
        
        allXSteps = [];
        allYSteps = [];
        allTurns = [];
        allNCounts = 0;
        allAcc=[];
        
        for jj = 1:size(r,2)
            valDXDY = intersect(find(valInDeltaX(:,jj)==cX), find(valInDeltaY(:,jj)==cY));
            xStep = (sf(valDXDY) + af(valDXDY)).*cos(turnf(valDXDY)); % the x component of the distance travelled immediately after the interaction is given by the speed between t and t+1
            yStep = (sf(valDXDY) + af(valDXDY)).*sin(turnf(valDXDY)); % the speed between t and t+1 is equal to the speed between t-1 and t + the acceleration
            allYSteps = [allYSteps; yStep];
            allXSteps = [allXSteps; xStep];
            allTurns = [allTurns; turnf(valDXDY)];
            allAcc = [allAcc; af(valDXDY)];
            allNCounts = allNCounts + length(valDXDY);
        end
        
        xTurnComp = nanmean(cos(allTurns));
        yTurnComp = nanmean(sin(allTurns));
        rho = sqrt(xTurnComp^2 + yTurnComp^2);
        meanTurnDXDY(cX,cY) = atan2(yTurnComp, xTurnComp);
        circSTDDXDY(cX,cY) = sqrt(-2*log(rho));
        meanAccDXDY(cX,cY) = nanmean(allAcc);
        nCountsDXDY(cX,cY) = allNCounts;
        meanYStepDXDY(cX,cY) =  nanmean(allYSteps);
        meanXStepDXDY(cX,cY) =  nanmean(allXSteps);
    end
end







for cTheta=1:length(thetaBins)-1
    cTheta
    for cR=1:length(rBins)-1
        
        allXSteps = [];
        allYSteps = [];
        allTurns = [];
        allNCounts = 0;
        allAcc=[];
        
        for jj = 1:size(r,2)
            valRTheta = intersect(find(valInR(:,jj) == cR), find(valInTheta(:,jj) == cTheta));
            xStep = (sf(valRTheta) + af(valRTheta)).*cos(turnf(valRTheta)); % the x component of the distance travelled immediately after the interaction is given by the speed between t and t+1
            yStep = (sf(valRTheta) + af(valRTheta)).*sin(turnf(valRTheta)); % the speed between t and t+1 is equal to the speed between t-1 and t + the acceleration
            allYSteps = [allYSteps; yStep];
            allXSteps = [allXSteps; xStep];
            allTurns = [allTurns; turnf(valRTheta)];
            allAcc = [allAcc; af(valRTheta)];
            allNCounts = allNCounts + length(valRTheta);
        end
        
        xTurnComp = nanmean(cos(allTurns));
        yTurnComp = nanmean(sin(allTurns));
        rho = sqrt(xTurnComp^2 + yTurnComp^2);
        meanTurnRTheta(cR,cTheta) = atan2(yTurnComp, xTurnComp);
        circSTDRTheta(cR,cTheta) = sqrt(-2*log(rho));
        meanAccRTheta(cR,cTheta) = nanmean(allAcc);
        nCountsRTheta(cR,cTheta) = allNCounts;
        meanYStepRTheta(cR, cTheta) =  nanmean(allYSteps);
        meanXStepRTheta(cR, cTheta) =  nanmean(allXSteps);
    end
    
    
    
    
    for cLR = 1:length(largeRBins) - 1
        for cPhi=1:length(phiBins)-1
            
            
            allXSteps = [];
            allYSteps = [];
            allTurns = [];
            allNCounts = 0;
            allAcc=[];
            
            for jj = 1:size(r,2)
                valLRPhiTheta = intersect(find(valInLargeR == cLR), intersect(find(valInPhi(:,jj) == cPhi), find(valInTheta(:,jj) == cTheta)));
                xStep = (sf(valLRPhiTheta) + af(valLRPhiTheta)).*cos(turnf(valLRPhiTheta)); % the x component of the distance travelled immediately after the interaction is given by the speed between t and t+1
                yStep = (sf(valLRPhiTheta) + af(valLRPhiTheta)).*sin(turnf(valLRPhiTheta)); % the speed between t and t+1 is equal to the speed between t-1 and t + the acceleration
                allYSteps = [allYSteps; yStep];
                allXSteps = [allXSteps; xStep];
                allTurns = [allTurns; turnf(valLRPhiTheta)];
                allAcc = [allAcc; af(valLRPhiTheta)];
                allNCounts = allNCounts + length(valLRPhiTheta);
            end
            xTurnComp = nanmean(cos(allTurns));
            yTurnComp = nanmean(sin(allTurns));
            rho = sqrt(xTurnComp^2 + yTurnComp^2);
            meanTurnPhiTheta(cPhi,cTheta,cLR) = atan2(yTurnComp, xTurnComp);
            circSTDPhiTheta(cPhi,cTheta,cLR) = sqrt(-2*log(rho));
            meanAccPhiTheta(cPhi,cTheta,cLR) = nanmean(allAcc);
            nCountsPhiTheta(cPhi,cTheta,cLR) = allNCounts;
            
            meanYStepPhiTheta(cPhi, cTheta, cLR) =  nanmean(allYSteps);
            meanXStepPhiTheta(cPhi, cTheta, cLR) =  nanmean(allXSteps);
        end
    end
    
    
end






%% figure for turning rate DX DY
cAxisLim = (max(max(abs(meanTurnDXDY)))) * framesPerTimeUnit;
figure,
set(gcf, 'Position', [1 1 1000 700]);
% imagesc(meanTurnDXDY* framesPerTimeUnit);
pcolor(deltaXBins,deltaYBins,meanTurnDXDY * framesPerTimeUnit);
caxis([-cAxisLim, cAxisLim]);
hold on;
fill([0, -0.25, 0, 0.25],[0.5, -0.4, 0, -0.4], 'w', 'LineWidth', 1, 'LineStyle', 'none');
axis equal xy tight;
box off;
colormap(myColourMap);

set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
    'FontName', 'Arial');


xlabel(sprintf('left-right distance to neighbour (%s)', spaceUnit)); ylabel(sprintf('front-back distance to neighbour (%s)', spaceUnit));
colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
ylabel(colorbarHandle, sprintf('turning rate (rad / %s)', timeUnit), 'FontSize', useFontSize);
colormap(myColourMap);

if isSavingFigures
    figureFileName = 'turning_rate_vs_deltaX_deltaY.eps';
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end








%% figure for acceleration DX DY
cAxisLim = max(max(abs(meanAccDXDY))) * framesPerTimeUnit^2;
figure,
set(gcf, 'Position', [1 1 1000 700]);
pcolor(deltaXBins,deltaYBins,meanAccDXDY * framesPerTimeUnit^2); caxis([-cAxisLim, cAxisLim]);
hold on;
fill([0, -0.25, 0, 0.25],[0.5, -0.4, 0, -0.4], 'w', 'LineWidth', 1, 'LineStyle', 'none');
axis equal xy tight;
box off;
colormap(myColourMap);

set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
    'FontName', 'Arial');


xlabel(sprintf('left-right distance to neighbour (%s)', spaceUnit)); ylabel(sprintf('front-back distance to neighbour (%s)', spaceUnit));
colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
ylabel(colorbarHandle, sprintf('acceleration (%s/%s^2)', spaceUnit, timeUnit), 'FontSize', useFontSize);
colormap(myColourMap);

if isSavingFigures
    figureFileName = 'acceleration_vs_deltaX_deltaY.eps';
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end





%% figure for number of counts DX DY
figure,
set(gcf, 'Position', [1 1 1000 700]);
pcolor(deltaXBins,deltaYBins,nCountsDXDY);
hold on;
fill([0, -0.25, 0, 0.25],[0.5, -0.4, 0, -0.4], 'w', 'LineWidth', 1, 'LineStyle', 'none');
axis equal xy tight;
box off;
colormap(myColourMap);

set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
    'FontName', 'Arial');


xlabel(sprintf('left-right distance to neighbour (%s)', spaceUnit)); ylabel(sprintf('front-back distance to neighbour (%s)', spaceUnit));
colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
ylabel(colorbarHandle, 'number of counts', 'FontSize', useFontSize);
colormap(myColourMap);

if isSavingFigures
    figureFileName = 'number_of_counts_vs_deltaX_deltaY.eps';
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end







X = rBins'*cos(pi/2 -thetaBins);
Y = rBins'*sin(pi/2 -thetaBins);


%% figure for turning rate R Theta
caxisLim = (max(max(abs(meanTurnRTheta)))) * framesPerTimeUnit;
figure,
set(gcf, 'Position', [1 1 1000 700]);
pcolor(X,Y,meanTurnRTheta * framesPerTimeUnit); caxis([-caxisLim, caxisLim]);
hold on;
fill([0, -0.25, 0, 0.25],[0.5, -0.4, 0, -0.4], 'w', 'LineWidth', 1, 'LineStyle', 'none');
axis equal xy tight;
box off;
colormap(myColourMap);

set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
    'XLim', [-max(rBins)*12/10 max(rBins)*12/10], 'YLim', [-max(rBins)*12/10, max(rBins)*12/10], ... 'XTickLabel', -max(distanceBins):1:max(distanceBins), 'YTickLabel',  -max(distanceBins):1:max(distanceBins), ...
    'FontName', 'Arial');
text(0,max(rBins)*11/10,'\theta = 0', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(max(rBins)*11/10,0,'\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(-max(rBins)*11/10,0,'-\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(0,-max(rBins)*11/10,'\pi', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');

text(-max(rBins), max(rBins), '(c)', 'FontName', 'Arial', 'FontSize', 30);


xlabel(sprintf('distance to neighbour (%s)', spaceUnit)); ylabel(sprintf('distance to neighbour (%s)', spaceUnit));
colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
ylabel(colorbarHandle, sprintf('turning rate (rad / %s)', timeUnit), 'FontSize', useFontSize);
colormap(myColourMap);

if isSavingFigures
    figureFileName = 'turning_rate_vs_R_theta.eps';
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end








%% figure for acceleration R Theta
caxisLim = max(max(abs(meanAccRTheta))) * framesPerTimeUnit^2;
figure,
set(gcf, 'Position', [1 1 1000 700]);
pcolor(X,Y,meanAccRTheta * framesPerTimeUnit^2); caxis([-caxisLim, caxisLim]);
hold on;
fill([0, -0.25, 0, 0.25],[0.5, -0.4, 0, -0.4], 'w', 'LineWidth', 1, 'LineStyle', 'none');
axis equal xy tight;
box off;
colormap(myColourMap);

set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
    'XLim', [-max(rBins)*12/10 max(rBins)*12/10], 'YLim', [-max(rBins)*12/10, max(rBins)*12/10], ... 'XTickLabel', -max(distanceBins):1:max(distanceBins), 'YTickLabel',  -max(distanceBins):1:max(distanceBins), ...
    'FontName', 'Arial');
text(0,max(rBins)*11/10,'\theta = 0', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(max(rBins)*11/10,0,'\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(-max(rBins)*11/10,0,'-\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(0,-max(rBins)*11/10,'\pi', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');

text(-max(rBins), max(rBins), '(b)', 'FontName', 'Arial', 'FontSize', 30);


xlabel(sprintf('distance to neighbour (%s)', spaceUnit)); ylabel(sprintf('distance to neighbour (%s)', spaceUnit));
colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
ylabel(colorbarHandle, sprintf('acceleration (%s/%s^2)', spaceUnit, timeUnit), 'FontSize', useFontSize);
colormap(myColourMap);

if isSavingFigures
    figureFileName = 'acceleration_vs_R_theta.eps';
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end





%% figure for number of counts R Theta
figure,
set(gcf, 'Position', [1 1 1000 700]);
pcolor(X,Y,nCountsRTheta);
hold on;
fill([0, -0.25, 0, 0.25],[0.5, -0.4, 0, -0.4], 'w', 'LineWidth', 1, 'LineStyle', 'none');
axis equal xy tight;
box off;
colormap(myColourMap);

set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
    'XLim', [-max(rBins)*12/10 max(rBins)*12/10], 'YLim', [-max(rBins)*12/10, max(rBins)*12/10], ... 'XTickLabel', -max(distanceBins):1:max(distanceBins), 'YTickLabel',  -max(distanceBins):1:max(distanceBins), ...
    'FontName', 'Arial');
text(0,max(rBins)*11/10,'\theta = 0', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(max(rBins)*11/10,0,'\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(-max(rBins)*11/10,0,'-\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(0,-max(rBins)*11/10,'\pi', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');

text(-max(rBins), max(rBins), '(a)', 'FontName', 'Arial', 'FontSize', 30);

xlabel(sprintf('distance to neighbour (%s)', spaceUnit)); ylabel(sprintf('distance to neighbour (%s)', spaceUnit));
colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
ylabel(colorbarHandle, 'number of counts', 'FontSize', useFontSize);
colormap(myColourMap);

if isSavingFigures
    figureFileName = 'number_of_counts_vs_R_theta.eps';
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end





for cLR = 1:size(meanTurnPhiTheta,3)
    
    %% Figure for turning rate Phi Theta
    
    caxisLim = (max(max(abs(meanTurnPhiTheta(:,:,cLR))))) * framesPerTimeUnit;
    if caxisLim == 0
        caxisLim = 0.000001;
    end
    figure,
    set(gcf, 'Position', [1 1 1000 640]);
    imagesc(meanTurnPhiTheta(1:end-1, 1:end-1, cLR) * framesPerTimeUnit , [-caxisLim, caxisLim]);
    
    
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
    ylabel(colorbarHandle, sprintf('turning (rad / %s^2)', timeUnit), 'FontSize', useFontSize);
    axis xy equal tight;
    title(sprintf('Neighbour distance r=%0.1f - %0.1f', largeRBins(cLR), largeRBins(cLR + 1)), 'FontName', 'Arial');
    colormap(myColourMap);
    
    annotation('textbox', [0.14, 0.9, 0, 0], 'string', '(a)', 'FontSize', 30, 'FontName', 'Arial');
    
    if isSavingFigures
        figureFileName = sprintf('turning_rate_vs_phi_theta_in_range%0.1f-%0.1f.eps', largeRBins(cLR), largeRBins(cLR + 1));
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
    
    
    
    
    
    
    
    
    %% Figure for acceleration Phi Theta
    caxisLim = (max(max(abs(meanAccPhiTheta(:,:,cLR))))) * framesPerTimeUnit^2;
    if caxisLim == 0
        caxisLim = 0.000001;
    end
    figure,
    set(gcf, 'Position', [1 1 1000 640]);
    imagesc(meanAccPhiTheta(1:end-1, 1:end-1, cLR) * framesPerTimeUnit^2 , [-caxisLim, caxisLim]);
    
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
    ylabel(colorbarHandle, sprintf('acceleration (%s/%s^2)', spaceUnit, timeUnit), 'FontSize', useFontSize);
    axis xy equal tight;
    title(sprintf('Neighbour distance r=%0.1f - %0.1f', largeRBins(cLR), largeRBins(cLR + 1)), 'FontName', 'Arial');
    colormap(myColourMap);
    
    annotation('textbox', [0.14, 0.9, 0, 0], 'string', '(b)', 'FontSize', 30, 'FontName', 'Arial');
    
    
    if isSavingFigures
        figureFileName = sprintf('acceleration_vs_phi_theta_in_range%0.1f-%0.1f.eps', largeRBins(cLR), largeRBins(cLR + 1));
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
    
    %% Figure for number of counts Phi Theta
    figure,
    set(gcf, 'Position', [1 1 1000 640]);
    imagesc(nCountsPhiTheta(1:end-1, 1:end-1, cLR));
    
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
    ylabel(colorbarHandle, 'number of counts', 'FontSize', useFontSize);
    axis xy equal tight;
    title(sprintf('Neighbour distance r=%0.1f - %0.1f', largeRBins(cLR), largeRBins(cLR + 1)), 'FontName', 'Arial');
    colormap(myColourMap);
    
    if isSavingFigures
        figureFileName = sprintf('number_of_counts_vs_phi_theta_in_range%0.1f-%0.1f.eps', largeRBins(cLR), largeRBins(cLR + 1));
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
    
end









%% if there are two neighbours, plot the relative influence of each of them
if size(r,2) == 2
    
    
    % bins for the figure
    deltaXBins = -15:0.5:15;
    deltaYBins = -15:0.5:15;
    
    [nInDeltaX, valInDeltaX] = histc(r .* cos(theta), deltaXBins);
    [nInDeltaY, valInDeltaY] = histc(r .* sin(theta), deltaYBins);
    
    
    clear meanYStepDXn1DXn2 meanXStepDXn1DXn2 meanTurnDXn1DXn2 circSTDDXn1DXn2 meanAccDXn1DXn2 nCountsDXn1DXn2;
    clear meanYStepDYn1DYn2 meanXStepDYn1DYn2 meanTurnDYn1DYn2 circSTDDYn1DYn2 meanAccDYn1DYn2 nCountsDYn1DYn2;
    
    % initialize empty arrays
    meanTurnDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins)); % x is front back, y is left right
    meanAccDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins));
    circSTDDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins));
    nCountsDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins));
    meanYStepDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins));
    meanXStepDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins));
    
    % initialize empty arrays
    meanTurnDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins)); % x is front back, y is left right
    meanAccDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins));
    circSTDDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins));
    nCountsDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins));
    meanYStepDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins));
    meanXStepDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins));
    
    
    [~,rankNeighbours] = sort(r,2);
    % some versions of Matlab don't like the previous line; use
    % instead:
    % [idontcareofthis,rankNeighbours] = sort(r,2);
    indexRankNeighbours = NaN(size(r));
    for jj = 1:size(r,2)
        indexRankNeighbours(:,jj) = sub2ind(size(r),(1:size(r, 1))', rankNeighbours(:,jj));
    end
    
    % first construct the maps in DeltaX and DeltaY
    for cXn1=1:length(deltaXBins) -1
        cXn1
        for cXn2=1:length(deltaXBins) -1
            valDXn1DXn2 = intersect(find(valInDeltaX(indexRankNeighbours(:,1))==cXn1), find(valInDeltaX(indexRankNeighbours(:,2))==cXn2));
            xTurnComp = nanmean(cos(turnf(valDXn1DXn2)));
            yTurnComp = nanmean(sin(turnf(valDXn1DXn2)));
            rho = sqrt(xTurnComp^2 + yTurnComp^2);
            meanTurnDXn1DXn2(cXn1,cXn2) = atan2(yTurnComp, xTurnComp);
            circSTDDXn1DXn2(cXn1,cXn2) = sqrt(-2*log(rho));
            meanAccDXn1DXn2(cXn1,cXn2) = nanmean(af(valDXn1DXn2));
            nCountsDXn1DXn2(cXn1,cXn2) = length(valDXn1DXn2);
            xStep = (sf(valDXn1DXn2) + af(valDXn1DXn2)).*cos(turnf(valDXn1DXn2));
            yStep = (sf(valDXn1DXn2) + af(valDXn1DXn2)).*sin(turnf(valDXn1DXn2));
            meanYStepDXn1DXn2(cXn1,cXn2) = nanmean(xStep);
            meanXStepDXn1DXn2(cXn1,cXn2) =  nanmean(yStep);
        end
    end
    
    
    
    % first construct the maps in DeltaX and DeltaY
    for cYn1=1:length(deltaYBins) -1
        cYn1
        for cYn2=1:length(deltaYBins) -1
            valDYn1DYn2 = intersect(find(valInDeltaY(indexRankNeighbours(:,1))==cYn1), find(valInDeltaY(indexRankNeighbours(:,2))==cYn2));
            xTurnComp = nanmean(cos(turnf(valDYn1DYn2)));
            yTurnComp = nanmean(sin(turnf(valDYn1DYn2)));
            rho = sqrt(xTurnComp^2 + yTurnComp^2);
            meanTurnDYn1DYn2(cYn1,cYn2) = atan2(yTurnComp, xTurnComp);
            circSTDDYn1DYn2(cYn1,cYn2) = sqrt(-2*log(rho));
            meanAccDYn1DYn2(cYn1,cYn2) = nanmean(af(valDYn1DYn2));
            nCountsDYn1DYn2(cYn1,cYn2) = length(valDYn1DYn2);
            xStep = (sf(valDYn1DYn2) + af(valDYn1DYn2)).*cos(turnf(valDYn1DYn2));
            yStep = (sf(valDYn1DYn2) + af(valDYn1DYn2)).*sin(turnf(valDYn1DYn2));
            meanYStepDYn1DYn2(cYn1,cYn2) = nanmean(xStep);
            meanXStepDYn1DYn2(cYn1,cYn2) =  nanmean(yStep);
        end
    end
    
    
    
    
    %% figure for turning rate DXn1 DXn2
    cAxisLim = 0.4 * (max(max(abs(meanTurnDXn1DXn2)))) * framesPerTimeUnit; % with 0.4 allow for some saturation
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    % imagesc(meanTurnDXn1DXn2* framesPerTimeUnit);
    pcolor(deltaXBins,deltaXBins,meanTurnDXn1DXn2 * framesPerTimeUnit);
    caxis([-cAxisLim, cAxisLim]);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('front-back distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('front-back distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, sprintf('turning rate (rad / %s)', timeUnit), 'FontSize', useFontSize);
    colormap(myColourMap);
    
    if isSavingFigures
        figureFileName = 'turning_rate_vs_deltaXn1_deltaXn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
    
    
    
    
    
    
    
    %% figure for acceleration DXn1 DXn2
    cAxisLim = 0.4 * max(max(abs(meanAccDXn1DXn2))) * framesPerTimeUnit^2; % with 0.4 allow for some saturation
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(deltaXBins,deltaXBins,meanAccDXn1DXn2 * framesPerTimeUnit^2); caxis([-cAxisLim, cAxisLim]);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('front-back distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('front-back distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, sprintf('acceleration (%s/%s^2)', spaceUnit, timeUnit), 'FontSize', useFontSize);
    colormap(myColourMap);
    hold on;
    [C,h]=contour(deltaXBins, deltaXBins, nCountsDXn1DXn2, [max(max(nCountsDXn1DXn2))/10, max(max(nCountsDXn1DXn2))/10]);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2, 'LineWidth', 1.5, 'LineColor', [0.5 0.5 0.5])
    text_handle = clabel(C,h);
    set(text_handle,'FontSize', 16, 'BackgroundColor',[1 1 1],'Edgecolor',[0 0 0], 'String', '10%')
    
    
    if isSavingFigures
        figureFileName = 'acceleration_vs_deltaXn1_deltaXn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
        save([figureFileName(1:end-4), '.mat'], 'deltaXBins', 'meanAccDXn1DXn2', 'nCountsDXn1DXn2', 'framesPerTimeUnit', 'spaceUnit', 'timeUnit');
    end
    
    
    
    
    
    %% figure for number of counts DXn1 DXn2
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(deltaXBins,deltaXBins,nCountsDXn1DXn2);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('front-back distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('front-back distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'number of counts', 'FontSize', useFontSize);
    colormap(myColourMap);
    
    
    
    if isSavingFigures
        figureFileName = 'number_of_counts_vs_deltaXn1_deltaXn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
    
    
    
    
    %% figure for turning rate DYn1 DYn2
    cAxisLim = 0.4 * (max(max(abs(meanTurnDYn1DYn2)))) * framesPerTimeUnit; % with 0.4 allow for some saturation
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    % imagesc(meanTurnDYn1DYn2* framesPerTimeUnit);
    pcolor(deltaYBins,deltaYBins,meanTurnDYn1DYn2 * framesPerTimeUnit);
    caxis([-cAxisLim, cAxisLim]);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('left-right distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('left-right distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, sprintf('turning rate (rad / %s)', timeUnit), 'FontSize', useFontSize);
    colormap(myColourMap);
    hold on;
    [C,h]=contour(deltaYBins, deltaYBins, nCountsDYn1DYn2, [max(max(nCountsDYn1DYn2))/10, max(max(nCountsDYn1DYn2))/10]);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2, 'LineWidth', 1.5, 'LineColor', [0.5 0.5 0.5])
    text_handle = clabel(C,h);
    set(text_handle,'FontSize', 16, 'BackgroundColor',[1 1 1],'Edgecolor',[0 0 0], 'String', '10%')
    
    
    if isSavingFigures
        figureFileName = 'turning_rate_vs_deltaYn1_deltaYn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
        save([figureFileName(1:end-4), '.mat'], 'deltaYBins', 'meanTurnDYn1DYn2', 'nCountsDYn1DYn2', 'framesPerTimeUnit', 'spaceUnit', 'timeUnit');
    end
    
    
    
    
    
    
    
    
    %% figure for acceleration DYn1 DYn2
    cAxisLim = 0.4 * max(max(abs(meanAccDYn1DYn2))) * framesPerTimeUnit^2; % with 0.4 allow for some saturation
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(deltaYBins,deltaYBins,meanAccDYn1DYn2 * framesPerTimeUnit^2); caxis([-cAxisLim, cAxisLim]);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('left-right distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('left-right distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, sprintf('acceleration (%s/%s^2)', spaceUnit, timeUnit), 'FontSize', useFontSize);
    colormap(myColourMap);
    
    if isSavingFigures
        figureFileName = 'acceleration_vs_deltaYn1_deltaYn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
    
    
    
    
    %% figure for number of counts DYn1 DYn2
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(deltaYBins,deltaYBins,nCountsDYn1DYn2);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('left-right distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('left-right distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'number of counts', 'FontSize', useFontSize);
    colormap(myColourMap);
    
    if isSavingFigures
        figureFileName = 'number_of_counts_vs_deltaYn1_deltaYn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
end














%% if there is only one neighbours, plot the expected influence for two neighbours computed by averaging the individual responses
if size(r,2) == 1
    
    
    % bins for the figure
    deltaXBins = -15:0.5:15;
    deltaYBins = -15:0.5:15;
    
    [nInDeltaX, valInDeltaX] = histc(r .* cos(theta), deltaXBins);
    [nInDeltaY, valInDeltaY] = histc(r .* sin(theta), deltaYBins);
    
    
    clear predictedMeanYStepDXn1DXn2 predictedMeanXStepDXn1DXn2 predictedMeanTurnDXn1DXn2 predictedMeanAccDXn1DXn2 predictedNCountsDXn1DXn2;
    clear predictedMeanYStepDYn1DYn2 predictedMeanXStepDYn1DYn2 predictedMeanTurnDYn1DYn2 predictedMeanAccDYn1DYn2 predictedNCountsDYn1DYn2;
    clear meanYStepDX meanXStepDX meanTurnDX circSTDDX meanAccDX nCountsDX;
    clear meanYStepDY meanXStepDY meanTurnDY circSTDDY meanAccDY nCountsDY;
    
    
    % initialize empty arrays
    meanTurnDX = NaN(length(deltaXBins), 1); % x is front back, y is left right
    meanAccDX = NaN(length(deltaXBins), 1);
    circSTDDX = NaN(length(deltaXBins), 1);
    nCountsDX = NaN(length(deltaXBins), 1);
    meanYStepDX = NaN(length(deltaXBins), 1);
    meanXStepDX = NaN(length(deltaXBins), 1);
    
    % initialize empty arrays
    meanTurnDY = NaN(length(deltaYBins), 1); % x is front back, y is left right
    meanAccDY = NaN(length(deltaYBins), 1);
    circSTDDY = NaN(length(deltaYBins), 1);
    nCountsDY = NaN(length(deltaYBins), 1);
    meanYStepDY = NaN(length(deltaYBins), 1);
    meanXStepDY = NaN(length(deltaYBins), 1);
    
    % initialize empty arrays
    predictedMeanTurnDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins)); % x is front back, y is left right
    predictedMeanAccDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins));
    predictedNCountsDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins));
    predictedMeanYStepDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins));
    predictedMeanXStepDXn1DXn2 = NaN(length(deltaXBins), length(deltaXBins));
    
    % initialize empty arrays
    predictedMeanTurnDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins)); % x is front back, y is left right
    predictedMeanAccDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins));
    predictedNCountsDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins));
    predictedMeanYStepDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins));
    predictedMeanXStepDYn1DYn2 = NaN(length(deltaYBins), length(deltaYBins));
    
    
    
    % first construct the maps in DeltaX
    for cX=1:length(deltaXBins) -1
        valDX = find(valInDeltaX==cX);
        xTurnComp = nanmean(cos(turnf(valDX)));
        yTurnComp = nanmean(sin(turnf(valDX)));
        rho = sqrt(xTurnComp^2 + yTurnComp^2);
        meanTurnDX(cX) = atan2(yTurnComp, xTurnComp);
        circSTDDX(cX) = sqrt(-2*log(rho));
        meanAccDX(cX) = nanmean(af(valDX));
        nCountsDX(cX) = length(valDX);
        xStep = (sf(valDX) + af(valDX)).*cos(turnf(valDX));
        yStep = (sf(valDX) + af(valDX)).*sin(turnf(valDX));
        meanYStepDX(cX) = nanmean(xStep);
        meanXStepDX(cX) =  nanmean(yStep);
    end
    
    
    
    % first construct the maps in DeltaY
    for cY=1:length(deltaYBins) -1
        valDY = find(valInDeltaY==cY);
        xTurnComp = nanmean(cos(turnf(valDY)));
        yTurnComp = nanmean(sin(turnf(valDY)));
        rho = sqrt(xTurnComp^2 + yTurnComp^2);
        meanTurnDY(cY) = atan2(yTurnComp, xTurnComp);
        circSTDDY(cY) = sqrt(-2*log(rho));
        meanAccDY(cY) = nanmean(af(valDY));
        nCountsDY(cY) = length(valDY);
        xStep = (sf(valDY) + af(valDY)).*cos(turnf(valDY));
        yStep = (sf(valDY) + af(valDY)).*sin(turnf(valDY));
        meanYStepDY(cY) = nanmean(xStep);
        meanXStepDY(cY) = nanmean(yStep);
    end
    
    
    % compose predicted effects of two identical neighoburs assuming
    % averaged interactions
    
    [meanYStepDXA, meanYStepDXB] = meshgrid(meanYStepDX);
    predictedMeanYStepDXn1DXn2 = (meanYStepDXA + meanYStepDXB)/2;
    [meanXStepDXA, meanXStepDXB] = meshgrid(meanXStepDX);
    predictedMeanXStepDXn1DXn2 = (meanXStepDXA + meanXStepDXB)/2;
    [meanTurnDXA, meanTurnDXB] = meshgrid(meanTurnDX);
    predictedMeanTurnDXn1DXn2 = atan2(sin(meanTurnDXA) + sin(meanTurnDXB), cos(meanTurnDXA) + cos(meanTurnDXB));
    [meanAccDXA, meanAccDXB] = meshgrid(meanAccDX);
    predictedMeanAccDXn1DXn2 = (meanAccDXA + meanAccDXB)/2;
    [nCountsDXA, nCountsDXB] = meshgrid(nCountsDX);
    predictedNCountsDXn1DXn2 = nCountsDXA + nCountsDXB;
    
    [meanYStepDYA, meanYStepDYB] = meshgrid(meanYStepDY);
    predictedMeanYStepDYn1DYn2 = (meanYStepDYA + meanYStepDYB)/2;
    [meanXStepDYA, meanXStepDYB] = meshgrid(meanXStepDY);
    predictedMeanXStepDYn1DYn2 = (meanXStepDYA + meanXStepDYB)/2;
    [meanTurnDYA, meanTurnDYB] = meshgrid(meanTurnDY);
    predictedMeanTurnDYn1DYn2 = atan2(sin(meanTurnDYA) + sin(meanTurnDYB), cos(meanTurnDYA) + cos(meanTurnDYB));
    [meanAccDYA, meanAccDYB] = meshgrid(meanAccDY);
    predictedMeanAccDYn1DYn2 = (meanAccDYA + meanAccDYB)/2;
    [nCountsDYA, nCountsDYB] = meshgrid(nCountsDY);
    predictedNCountsDYn1DYn2  = nCountsDYA + nCountsDYB;
    
    
    
    %% figure for turning rate DXn1 DXn2
    cAxisLim = 0.8 * (max(max(abs(predictedMeanTurnDXn1DXn2)))) * framesPerTimeUnit; % with 0.8 allow for some saturation
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    % imagesc(meanTurnDXn1DXn2* framesPerTimeUnit);
    pcolor(deltaXBins,deltaXBins,predictedMeanTurnDXn1DXn2 * framesPerTimeUnit);
    caxis([-cAxisLim, cAxisLim]);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('front-back distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('front-back distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, sprintf('predicted turning rate (rad / %s)', timeUnit), 'FontSize', useFontSize);
    colormap(myColourMap);
    
    if isSavingFigures
        figureFileName = 'predicted_turning_rate_vs_deltaXn1_deltaXn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
    
    
    
    
    
    
    
    %% figure for acceleration DXn1 DXn2
    cAxisLim = 0.8 * max(max(abs(predictedMeanAccDXn1DXn2))) * framesPerTimeUnit^2; % with 0.8 allow for some saturation
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(deltaXBins,deltaXBins,predictedMeanAccDXn1DXn2 * framesPerTimeUnit^2); caxis([-cAxisLim, cAxisLim]);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('front-back distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('front-back distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, sprintf('acceleration (%s/%s^2)', spaceUnit, timeUnit), 'FontSize', useFontSize);
    colormap(myColourMap);
    hold on;
    [C,h]=contour(deltaXBins, deltaXBins, predictedNCountsDXn1DXn2, [max(max(predictedNCountsDXn1DXn2))/10, max(max(predictedNCountsDXn1DXn2))/10]);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2, 'LineWidth', 1.5, 'LineColor', [0.5 0.5 0.5])
    text_handle = clabel(C,h);
    set(text_handle,'FontSize', 16, 'BackgroundColor',[1 1 1],'Edgecolor',[0 0 0], 'String', '10%')
    
    
    if isSavingFigures
        figureFileName = 'predicted_acceleration_vs_deltaXn1_deltaXn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
        save([figureFileName(1:end-4), '.mat'], 'deltaXBins', 'predictedMeanAccDXn1DXn2', 'predictedNCountsDXn1DXn2', 'framesPerTimeUnit', 'spaceUnit', 'timeUnit');
    end
    
    
    
    
    
    %% figure for number of counts DXn1 DXn2
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(deltaXBins,deltaXBins,predictedNCountsDXn1DXn2);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('front-back distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('front-back distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'number of counts', 'FontSize', useFontSize);
    colormap(myColourMap);
    
    
    
    if isSavingFigures
        figureFileName = 'predicted_number_of_counts_vs_deltaXn1_deltaXn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
    
    
    
    
    %% figure for turning rate DYn1 DYn2
    cAxisLim = 0.8 * (max(max(abs(predictedMeanTurnDYn1DYn2)))) * framesPerTimeUnit; % with 0.8 allow for some saturation
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    % imagesc(meanTurnDYn1DYn2* framesPerTimeUnit);
    pcolor(deltaYBins,deltaYBins,predictedMeanTurnDYn1DYn2 * framesPerTimeUnit);
    caxis([-cAxisLim, cAxisLim]);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('left-right distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('left-right distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, sprintf('turning rate (rad / %s)', timeUnit), 'FontSize', useFontSize);
    colormap(myColourMap);
    hold on;
    [C,h]=contour(deltaYBins, deltaYBins, predictedNCountsDYn1DYn2, [max(max(predictedNCountsDYn1DYn2))/10, max(max(predictedNCountsDYn1DYn2))/10]);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2, 'LineWidth', 1.5, 'LineColor', [0.5 0.5 0.5])
    text_handle = clabel(C,h);
    set(text_handle,'FontSize', 16, 'BackgroundColor',[1 1 1],'Edgecolor',[0 0 0], 'String', '10%')
    
    
    if isSavingFigures
        figureFileName = 'predicted_turning_rate_vs_deltaYn1_deltaYn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
        save([figureFileName(1:end-4), '.mat'], 'deltaYBins', 'predictedMeanTurnDYn1DYn2', 'predictedNCountsDYn1DYn2', 'framesPerTimeUnit', 'spaceUnit', 'timeUnit');
    end
    
    
    
    
    
    
    
    
    %% figure for acceleration DYn1 DYn2
    cAxisLim = 0.8 * max(max(abs(predictedMeanAccDYn1DYn2))) * framesPerTimeUnit^2; % with 0.8 allow for some saturation
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(deltaYBins,deltaYBins,predictedMeanAccDYn1DYn2 * framesPerTimeUnit^2); caxis([-cAxisLim, cAxisLim]);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('left-right distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('left-right distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, sprintf('acceleration (%s/%s^2)', spaceUnit, timeUnit), 'FontSize', useFontSize);
    colormap(myColourMap);
    
    if isSavingFigures
        figureFileName = 'predicted_acceleration_vs_deltaYn1_deltaYn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
    
    
    
    
    %% figure for number of counts DYn1 DYn2
    figure,
    set(gcf, 'Position', [1 1 1000 700]);
    pcolor(deltaYBins,deltaYBins,predictedNCountsDYn1DYn2);
    axis equal xy tight;
    box off;
    colormap(myColourMap);
    
    set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
        'FontName', 'Arial');
    
    
    xlabel(sprintf('left-right distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('left-right distance to neighbour 1 (%s)', spaceUnit));
    colorbarHandle = colorbar;
    set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
    ylabel(colorbarHandle, 'number of counts', 'FontSize', useFontSize);
    colormap(myColourMap);
    
    if isSavingFigures
        figureFileName = 'predicted_number_of_counts_vs_deltaYn1_deltaYn2.eps';
        print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
    end
    
    
    
    
    
    
    
    % If possible compute the non-pairwise interactions
    
    twoNeighboursTurningMap = 'turning_rate_vs_deltaYn1_deltaYn2.mat';
    if exist(twoNeighboursTurningMap, 'file')
        load(twoNeighboursTurningMap);
        try
            diffFromPredictedMeanTurnDYn1DYn2 = mod(pi + meanTurnDYn1DYn2 - predictedMeanTurnDYn1DYn2, 2*pi) - pi;
            
            %% figure for turning rate DYn1 DYn2
            cAxisLim = 0.7; % 0.8 * (max(max(abs(diffFromPredictedMeanTurnDYn1DYn2)))) * framesPerTimeUnit; % with 0.8 allow for some saturation
            figure,
            set(gcf, 'Position', [1 1 1000 700]);
            % imagesc(meanTurnDYn1DYn2* framesPerTimeUnit);
            pcolor(deltaYBins,deltaYBins,diffFromPredictedMeanTurnDYn1DYn2 * framesPerTimeUnit);
            caxis([-cAxisLim, cAxisLim]);
            axis equal xy tight;
            box off;
            colormap(myColourMap);
            
            set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
                'FontName', 'Arial');
            
            
            xlabel(sprintf('left-right distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('left-right distance to neighbour 1 (%s)', spaceUnit));
            colorbarHandle = colorbar;
            set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%, 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
            ylabel(colorbarHandle, sprintf('turning rate real - predicted(rad / %s)', timeUnit), 'FontSize', useFontSize);
            colormap(myColourMap);
            hold on;
            [C,h]=contour(deltaYBins, deltaYBins, nCountsDYn1DYn2, [max(max(nCountsDYn1DYn2))/10, max(max(nCountsDYn1DYn2))/10]);
            set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2, 'LineWidth', 1.5, 'LineColor', [0.5 0.5 0.5])
            text_handle = clabel(C,h);
            set(text_handle,'FontSize', 16, 'BackgroundColor',[1 1 1],'Edgecolor',[0 0 0], 'String', '10%')
            
            
            if isSavingFigures
                figureFileName = 'diff_between_observed_and_predicted_turning_rate_vs_deltaYn1_deltaYn2.eps';
                print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
            end
        catch
            warndlg('couldn''t compute the observed minus predicted turning for two neighbours');
        end
    end
    
    
    
    twoNeighboursAccMap = 'acceleration_vs_deltaXn1_deltaXn2.mat';
    if exist(twoNeighboursAccMap, 'file')
        load(twoNeighboursAccMap);
        try
            diffFromPredictedMeanAccDXn1DXn2 = meanAccDXn1DXn2 - predictedMeanAccDXn1DXn2;
            
            cAxisLim = 5; % 0.8 * max(max(abs(diffFromPredictedMeanAccDXn1DXn2))) * framesPerTimeUnit^2; % with 0.8 allow for some saturation
            figure,
            set(gcf, 'Position', [1 1 1000 700]);
            pcolor(deltaXBins,deltaXBins,diffFromPredictedMeanAccDXn1DXn2 * framesPerTimeUnit^2); caxis([-cAxisLim, cAxisLim]);
            axis equal xy tight;
            box off;
            colormap(myColourMap);
            
            set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', ... 'XTick', -max(distanceBins):1:max(distanceBins), 'YTick', -max(distanceBins):1:max(distanceBins), ...
                'FontName', 'Arial');
            
            
            xlabel(sprintf('front-back distance to neighbour 2 (%s)', spaceUnit)); ylabel(sprintf('front-back distance to neighbour 1 (%s)', spaceUnit));
            colorbarHandle = colorbar;
            set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');% , 'YTick', linspace(-valueForColourScale, valueForColourScale, 11));
            ylabel(colorbarHandle, sprintf('acceleration real - predicted (%s/%s^2)', spaceUnit, timeUnit), 'FontSize', useFontSize);
            colormap(myColourMap);
            hold on;
            [C,h]=contour(deltaXBins, deltaXBins, nCountsDXn1DXn2, [max(max(nCountsDXn1DXn2))/10, max(max(nCountsDXn1DXn2))/10]);
            set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2, 'LineWidth', 1.5, 'LineColor', [0.5 0.5 0.5])
            text_handle = clabel(C,h);
            set(text_handle,'FontSize', 16, 'BackgroundColor',[1 1 1],'Edgecolor',[0 0 0], 'String', '10%')
            
            
            if isSavingFigures
                figureFileName = 'diff_between_observed_and_predicted_acceleration_vs_deltaXn1_deltaXn2.eps';
                print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
            end
            warndlg('The figures for real - predicted are based on the saved data for three particles');

        catch
            warndlg('couldn''t compute the observed minus predicted acceleration for two neighbours');
        end
    end
end

