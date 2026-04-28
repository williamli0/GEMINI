%% Housekeeping
clearvars
close all

%% Plotting MCP Data
% Make sure the folder with the .csv files is on the current Matlab path

files = dir('*.csv'); % struct with a list of all the files in the folder

% Iterate over every file in the folder
for file = files'
    % Rudimentary parsing of the file name
    fileName = convertCharsToStrings(file.name);
    % Read the data from the .csv into an array
    data = readmatrix(fileName);
    a1 = data(:,1); % 1st column = x location
    a2 = data(:,2); % 2nd column = y location
    counts = data(:,3); % 3rd column = counts
    
    % Create scatter plot
    mcpMap = figure(); % Create new figure
    
    % Create scatter plot with counts as the third (color) axis
    scatter(a1, a2,[],counts,'filled')
    a=colorbar;
    a.Label.String = 'Counts';
    colormap turbo
    xlim([0 60])
    ylim([0 60])
    xlabel('a1 (MCP Pixels)')
    ylabel('a2 (MCP Pixels)')
    title(fileName) % Insert file name as figure title (in progress)
    
    % Make plot square
    axesHandles = findobj(get(mcpMap,'Children'), 'flat','Type','axes');
    axis(axesHandles,'square')

    % Save plot
    % INSIDE FOR LOOP = WILL SAVE A GRAPH FOR ALL FILES IN THE FOLDER
    % Basically be careful of file save location before running this script
    % Recommend creating a folder (ex. 'Graphs_Data_05.08.2024" for graphs
    % UNCOMMENT FOLLOWING LINE TO ENABLE SAVING
    %saveas(gcf,'/Users/williamli/Desktop/JP/JP_Testing_FullDay2/Day2Graphs/' + fileName + '.png')

%% Angular Scattering Distribution

    xBin = 25:35; % Range of x values to search over for a bin
    
    % Indices of a1 with x values in the bin of interest
    indices = find(a1>=xBin(1) & a1<=xBin(end));
    
    % Actual values of a1, a2, counts with x values in the bin of interest
    xLocations = a1(indices);
    yLocations = a2(indices);
    countsBin = counts(indices);

    % Combine counts at each y value (if yLocations is the same)
    groupIndex = cumsum(abs(diff([nan;yLocations]))~=0); % Group index
    
    % Total counts summed at each y value in the bin (y axis of bar graph)
    dataCombined = accumarray(groupIndex,countsBin,[]); % y axis
    numCounts = sum(dataCombined);
    
    % Duplicates removed from y locations = ...
    % yCombined corresponds with dataCombined 
    yCombined = unique(yLocations); % x axis
    
    % Converting MCP pixels to scattering angle
    xPsi = atand((yCombined-27)/38.0733); 

    figure()
    bar(xPsi, dataCombined) % Was 28.5 but calibrated manually to 27
    title('Binned Data at x = [' + string(xBin(1)) + ':' + string(xBin(end))' + '], ' + fileName)
    xlabel('Scattering Angle \psi')
    ylabel('Number of Counts')
    hold on
    grid on
    grid minor
    yline(max(dataCombined)/2) % Identifying half-max
    hold off
end


%% Charge State Separation
    yBin = 24:31; % Range of y values to search over for a bin
    
    % Indices of a2 with y values in the bin of interest
    indices2 = find(a2>=yBin(1) & a2<=yBin(end));
    
    % Actual values of a1, a2, counts with y values in the bin of interest
    xLocations2 = a1(indices2);
    yLocations2 = a2(indices2);
    countsBin2 = counts(indices2);
    
    % Duplicates removed from x locations = ...
    % xCombined corresponds with dataCombined 
    xCombined = unique(xLocations2); % y axis
    
    % Empty vector to store counts at each x location
    countsX = zeros(length(xCombined),1);

    for j = 1:length(xCombined)
        for i = 1:length(xLocations2)
            if xLocations2(i) == xCombined(j)
                countsX(j) = countsX(j) + countsBin2(i);
            end
        end
    end

    figure()
    bar(xCombined, countsX)
    title('Binned Data at y = [' + string(yBin(1)) + ':' + string(yBin(end))' + '], ' + fileName)
    xlabel('MCP X (a1) Values')
    ylabel('Number of Counts')


%% Fitting Angular Scattering with  Allegrini et al. 2026
allegrini = figure();

% Data for 0.6 ug/cm^2 carbon foil from Zenodo
% https://zenodo.org/records/18762535
AllegriniX = [5 10 15 20 3 5 10 15 20 5 10 15 20 1 2 3 5 1 7 1.3 5 10 15 20 5 10 15 2 5 10 15 20 5 10 15 20 5 5 10 10 15 20 5 10 15 20].';
AllegriniY = [4.932739 2.522655 1.68352 1.276964 8.46135 5.087358 2.574845 1.698804 1.270938 5.046015 2.568519 1.692684 1.270699 20.043041 11.649077 8.322569 5.036458 20.54845 3.635627 16.102651 5.163912 2.625565 1.753284 1.337823 4.955 2.52476 1.700405 13.269216 4.921402 2.449638 1.653521 1.252816 5.521133 2.576174 1.753599 1.306316 4.888065 4.828199 2.39198 2.410782 1.552558 1.170271 4.766434 2.374497 1.57705 1.200405].';

% Power Law fit to angular scattering distribution
f=fit(AllegriniX,AllegriniY,'power1');
plot(f)

axesHandles = findobj(get(allegrini,'Children'), 'flat','Type','axes');
axis(axesHandles,'square')
xlabel('Incident Energy (keV)')
ylabel('Scattering Angle Half-Width \psi_{1/2}')
title('Scattering Half-Width vs. Incident Energy for 0.6 \mug/cm^2 Carbon Foil')
xscale("log")
yscale("log")
xlim([1 14])
ylim([1 14])
hold on
scatter(AllegriniX, AllegriniY,70,'red','+');
scatter(4, 6, 'filled','o','k') % GEMINI Data
legend("Power Law Fit 20.78x^{-0.8866}",'Allegrini et al. (2026)','This work')