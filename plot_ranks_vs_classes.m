

function [] = plot_ranks_vs_classes(consideredVariables, classes, selectedVariable, isSavingFigures, xLabel, yLabel, classesLabel, textLegend)
% function plot_ranks_vs_classes
% given several independent measures of a variable, this plot sorts the
% measures in ascending order and plots them. The markers are coloured
% according to the classes to which each measure belongs.
% For example, if the variable is the directional correlation delay,
% measured for 120 flocks of two particles, the function plots these 120
% values, sorted in ascending order, with dots coloured according to the
% class to which the focal particle belongs (e.g. if the focal was from a
% given population, sex group, or had a preference for a specific position 
% within the flock).
%
% INPUT:
% consideredVariables
% has a column for each observed event and a row for each measure taken of 
% that event. 
% classes has the same number of columns as consideredVariables, and one
% single line, with the "class" of the focal individual, expressed by
% numbers
% selectedVariable
% indicates which row of consideredVariables will be sorted
% isSavingFigures
% specify if the resulting figure should be saved as eps or not
% xLabel
% label for the x axis
% yLabel
% label for the y axis
% classesLabel
% y label for the smaller plot with the colours used for the different
% classes
% textLegend
% text legend to be added to the main plot
%
% Written by:
% Andrea Perna
% http://www.perna.fr
%
% Date:
% 2014 / 04 / 18


% one marker for each measure
markerSelection = 'o^sd';
myMarkers = markerSelection(mod((1:size(consideredVariables,1))-1, length(markerSelection)) +1);

% sort observed events according to selectedVariable
[~, sortingIndices] = sort(consideredVariables(selectedVariable,:));
sortedVariables = consideredVariables(:, sortingIndices);

% if there are missing classes, find how many are still present
[uniqueClasses, iClasses, iUniqueClasses] = unique(classes);
nClasses = length(uniqueClasses);

% a colour for each unique class; classes are also sorted to match the
% sorting of the variables
markerColours = lines(nClasses);
sortedColours = markerColours(iUniqueClasses(sortingIndices),:); % one color for each class
sortedColoursBright = sortedColours * 0.7 + 0.3*ones(size(sortedColours));

% find the observations for selectedVariable that belong to each class
% if consideredVariables is a matrix, it is important to know what is in
% each line: it doesn't always make sense to sort according to
% selectedVariable and plot everything in the same way!
valInClass = cell(0);
for jj = 1:nClasses
    v = find(iUniqueClasses==jj);
    valInClass{jj} = consideredVariables(selectedVariable,v);
end

try
    [p,h] = ranksum(valInClass{1},valInClass{2}) % compare only the first two classes
catch
    warndlg('Problem computing ranksum in plot_ranks_vs_classes; perhaps there is only one class?', 'Problem');
end

clear hplot;
figure,
set(gcf, 'Position', [1 1 1000 500]);
subplot(1,20,1:16);
hold on;
for jj = 1:size(consideredVariables,1) % for each variable
    if jj == selectedVariable % plot slightly different the selectedVariable
        for kk = 1:size(consideredVariables,2) % for each data point
            hplot = plot(kk, sortedVariables(jj, kk), 'Color', sortedColours(kk,:), 'Marker', myMarkers(jj), 'MarkerFaceColor', sortedColours(kk,:), 'MarkerEdgeColor', sortedColours(kk,:)*0.5, 'MarkerSize', 12, 'LineWidth', 1.5);
        end
    else
        for kk = 1:size(consideredVariables,2) % for each data point
            plot(kk, sortedVariables(jj, kk), 'Color', sortedColoursBright(kk,:), 'Marker', myMarkers(jj), 'MarkerFaceColor', sortedColoursBright(kk,:), 'MarkerEdgeColor', sortedColoursBright(kk,:)*0.5, 'MarkerSize', 12, 'LineWidth', 1.5);
        end
    end
end
% set(gca,'XGrid','on', 'XMinorGrid', 'on', 'XMinorTick', 'on', 'MinorGridLineStyle', '-', 'TickDir', 'out')
set(gca, 'FontName', 'Arial', 'FontSize', 18);
% set(gca, 'XLim', [0, 31], 'XTick', 0:5:30, 'XTickLabel', 0:5:30);
set(gca, 'YGrid', 'on')
% set(gca, 'YLim', [0.5, 5]);
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca,'Layer','top')
%  set(0,'DefaultAxesLayer','top')
% box on;
xlabel(gca, xLabel);
ylabel(gca, yLabel);

legend(hplot, textLegend, 'Location','NorthWest');


subplot(1,20,19:20)
hold on;
for kk = 1:nClasses
    plot(0, uniqueClasses(kk), 'Marker', myMarkers(selectedVariable), 'MarkerFaceColor', markerColours(kk,:), 'MarkerEdgeColor', markerColours(kk,:)*0.5, 'MarkerSize', 12);
end
set(gca, 'FontName', 'Arial', 'FontSize', 18);
set(gca, 'YGrid', 'off')
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca, 'YTick', 1:nClasses);
% if nClasses == 2 
%    set(gca, 'YTickLabel', {'back', 'front'}); % back and front are only true when theta == 0
% end
set(gca, 'YLim', [1 - (nClasses-1)*0.1, nClasses + (nClasses-1)*0.1]);
set(gca, 'XTick', []);
set(gca,'Layer','top')
hyl = ylabel(classesLabel);
set(hyl,'Position', get(hyl,'Position') + [0.5 0 0])

if isSavingFigures
    figureFileName = 'figure_directional_correlation_delay.eps';
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end



% figure with the histogram
figure,
set(gcf, 'Position', [1 1 1000 500]);
bph = boxplot(consideredVariables(selectedVariable,:), iUniqueClasses, 'Notch', 'on');
set(gca, 'LineWidth', 1.5);
set(gca, 'TickDir', 'out');
set(gca, 'FontName', 'Arial', 'FontSize', 18);
% if nClasses == 2
%     set(gca, 'XTick', 1:nClasses, 'XTickLabel', {'focal is behind', 'focal is in front'});
% end
ylabel(gca, yLabel);

 
h=findobj(gca,'tag','Outliers'); % Get handles for outlier lines.
set(h,'Marker','o', 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); % Change symbols for all the groups.
set(bph,'LineWidth',2, 'Color', 'k', 'LineStyle', '-');

if isSavingFigures
    figureFileName = 'figure_directional_correlation_delay_boxplot.eps';
    print(gcf, '-depsc2', '-tiff', '-r300', figureFileName);
end



