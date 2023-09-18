% Set the file path prefix for the data
prefix = 'Z:/User/pc1aod/'; % Change this to the appropriate path
% subNum = 1; % Uncomment this line if subNum is being used
% chan = 1; % Uncomment this line if chan is being used

% HPC:
% If running on an HPC, use this prefix instead
% prefix = '/shared/dede_group/User/pc1aod/';

% Add paths to necessary code directories
addpath([prefix 'CODE/GEDbounds_clusterImprove'])
addpath([prefix 'CODE/subNetworkDynamics'])
addpath([prefix 'CODE/export_fig_repo'])

% Set the directory for saving figures
savedir = [prefix 'FIGURES/'];

% Set the directory for saving summary data
summaryDatSave = [prefix 'SUMDAT/'];

%% Stitch together results from individual subject runs

% Get a list of filenames in the specified directory that have a .mat extension
filenames = dir(fullfile(summaryDatSave, '*.mat'));

% Load data from the second file to initialize the 'allData' structure
allData = load([filenames(2).folder '/' filenames(2).name]).data;
allFields = fieldnames(allData);

% Filter filenames to only include those that start with 'femaleASD'
filteredstruct = filenames(startsWith({filenames.name}, 'femaleASD'));

% Load example data from a specific file for reference
exampleDat = load('Z:\User\pc1aod\SUMDAT\biomarkCon_1_1.mat').data;
labels = exampleDat.labels;
Th = exampleDat.Th;
Rd = exampleDat.Rd;

% Loop through the filtered filenames and process each file
for ii = 1:length(filteredstruct)
    ii
    
    %% Load the file
    data = load([filteredstruct(ii).folder '/' filteredstruct(ii).name]).data;

    %% Basic demographics
    % Data set
    filteredstruct(ii).dataSet = data.dataSet;
    % Sex
    filteredstruct(ii).sex = data.sex;
    % Age
    filteredstruct(ii).age = data.age;
    % Eyes
    filteredstruct(ii).eyes = data.eyes;
    % Group
    filteredstruct(ii).group = data.group;
    % Sub ID
    filteredstruct(ii).key = data.key;
    % Original file name
    filteredstruct(ii).origfile = data.fn;
    % Original trial count
    filteredstruct(ii).nbTrialOrig = data.nbTrialOrig;
    % Final trial count
    filteredstruct(ii).nbTrialFinal = data.nbTrialFinal;
    % Original channel count
    filteredstruct(ii).nbChanOrig = data.nbChanOrig;
    % Final channel count
    filteredstruct(ii).nbChanFinal = data.nbChanFinal;

    % Switch based on the data set to extract additional information
    switch data.dataSet
        case 'femaleASD'
            filteredstruct(ii).IQ = data.iq1;
            filteredstruct(ii).IQ_measure = data.iq1_measure;
    end

    %% Find matched key
    index{ii} = find(strcmp({filteredstruct.key}, filteredstruct(ii).key) == 1);
end

% index= arrayfun(@(y) find(arrayfun(@(x) strcmp(filteredstruct(x).key, filteredstruct(y).key), ...
%     find(name))), ...
%     find(name), 'uniformoutput', false);

for jj = 1:length(index)
    curIndex = index{jj};
    
    % Check if the current index has exactly two elements
    if length(curIndex) == 2
        % Load data from the two specified files
        temp = open([summaryDatSave filteredstruct(curIndex(1)).name]);
        temp2 = open([summaryDatSave filteredstruct(curIndex(2)).name]);

        % Create new fields in the first file's structure to store the loaded data
        filteredstruct(curIndex(1)).(join([filteredstruct(curIndex(1)).eyes 'filepairpath'])) = temp;
        filteredstruct(curIndex(1)).(join([filteredstruct(curIndex(2)).eyes 'filepairpath'])) = temp2;
    end
end

% Remove rows with empty cells from the 'filteredstruct'
temp = filteredstruct;
out = temp(all(~cellfun(@isempty, struct2cell(temp))));

condition = {'open', 'closed'};

for ii = 1:length(out)
    %% Regional sample entropy
    for jj = 1:length(condition)

        % Extract data and calculate relevant statistics
        data = out(ii).([condition{jj} 'filepairpath']).data;
        sampEnt = data.sampEnt;
        fuzzeyEnt = data.fuzEnt;
        fuzzeyEnt(fuzzeyEnt == 0) = NaN;

%         meanSE = squeeze(mean(sampEnt, 2));
%         out(ii).([condition{jj} 'MeanSE']) = meanSE;
% 
%                 meanFE  = squeeze(mean(fuzzeyEnt, 2, "omitnan"));
%                 out(ii).([condition{jj} 'MeanFE']) = meanFE;
        
        % Calculate relative power and store it in the 'out' structure
        pow = data.power;
        relPow = pow ./ sum(pow, 2);
        out(ii).([condition{jj} 'relPow']) = squeeze(mean(relPow, 3));

    end
end


%% separate data file based on age
out(16) = [];
out(22) = [];
out(23) = [];

% age subgroup
ageBound = [97; 135; 170; 215];
age1out = out([out.age] > ageBound(1) & [out.age] < ageBound(2));
age2out = out([out.age] > ageBound(2) & [out.age] < ageBound(3));
age3out = out([out.age] > ageBound(3) & [out.age] < ageBound(4));

AscDat = struct;
ConDat = struct;
for ii = 1:length(age3out)
    strtemp = strcmpi(age3out(ii).group,"CON");
    if strtemp == 0
        AscDat(ii).dat = age3out(ii);
    else
        ConDat(ii).dat = age3out(ii);
    end

end

temp = AscDat;
AscDatN = temp(~cellfun(@isempty,struct2cell(temp)));

temp = ConDat;
ConDatN = temp(~cellfun(@isempty,struct2cell(temp)));

for aa = 1:length(AscDatN)
    AscFileclosed(aa,:,:) = AscDatN(aa).dat.closedrelPow;
end


for cc = 1:length(ConDatN)
    ConFileclosed(cc,:,:) = ConDatN(cc).dat.closedrelPow;
end

for aa = 1:length(AscDatN)
    AscFileopen(aa,:,:) = AscDatN(aa).dat.openrelPow;
end


for cc = 1:length(ConDatN)
    ConFileopen(cc,:,:) = ConDatN(cc).dat.openrelPow;
end
DiffOC_ASC = AscFileopen - AscFileclosed;
DiffOC_Con = ConFileopen - ConFileclosed;

Va = ConDatN;
for aa = 1:length(Va)
 Temp(aa,:) = Va(aa).dat.IQ;
end
fieldToAverage = [Temp];
mean(fieldToAverage)
std(fieldToAverage)
% out1=rmfield(out,{'openfilepairpath';'closedfilepairpath';'eyes';'name';'bytes';'folder';'date';'isdir';'datenum';'dataSet'});
% writetable(struct2table(out1), 'E:\Research_update\eeg\bigdat_entropy\autismBiomarkersAllDatafemaleASD_withallscale.csv')


%% Reordering electrode data to plot a nice brain result figure

% Define the original and last ordered electrode arrays
original = {'Fp1','AF3','F7','F3','FC1','FC5','T7','C3','CP1','CP5','P7','P3','Pz','PO3','O1','Oz','O2','PO4','P4','P8','CP6','CP2','C4','T8','FC6',...,
    'FC2','F4','F8','AF4','Fp2','Fz','Cz'};
lastordered = {'Fp1','AF3','F3','F7','Fp2','AF4','F4','F8','Fz','FC1','FC5','C3','T7','CP1','CP5',...,
    'FC2','FC6','C4','T8','CP2','CP6','Cz','P3','P7','PO3','O1','Pz','P4','P8','PO4','O2','Oz'};

% Initialize an array to store the index of each electrode in 'original'
for ii = 1:length(lastordered)
    ordIdx(ii) = find(strcmp(lastordered{ii}, original));
end

% Define electrode index mappings
ElecIdx = [1	2	3	4	27	28	29	30	31	5	6	8	7	9	10	21	22	23	24	25	26	32	11	12	13	14	15	16	17	18	19	20];

% Initialize variables
clear Summary Summary1

% Reorder the data based on electrode index mappings
for ii = 1:32
    Summary1(:, ii, :) = Summary(:, ordIdx(ii), :);
end

% Update the 'Summary' variable with the reordered data
clear Summary
Summary = Summary1;

% Calculate the difference between open and closed conditions for ASC and CON
DiffOC_ASC = AscFileOpenRelpow - AscFileClosedRelpow;
DiffOC_Con = ConFileOpenRelpow - ConFileClosedRelpow;


% Define the number of unique colors
n = 2048;

% Define a base color map with color values
base = [0 0 .5; % dark blue
        0 1 0;   % light green
        1 1 0;   % yellow
        1 0.5 0; % orange
        1 0 0];  % red

% Interpolate the base color map to create a custom color map
customMap = interp1(linspace(n, 0, size(base, 1)), base, fliplr(0:n), 'pchip');

% Specify the directory for saving figures
saveDir = 'E:\Research_update\eeg\bigdat_entropy\orderedFiles\lastordered\figure\';

% Create a figure and set its background color to white
figure,
set(gcf, 'color', 'w')

% Display an image with the TFCE_openASC31.TFCE_Obs data using the custom color map
imagesc(TFCE_openASC31.TFCE_Obs)
colormap(customMap)

% Loop through electrode rows and highlight significant points
for elecs = 1:size(TFCE_openASC31.TFCE_Obs, 1)
    hold on;
    sig = find(TFCE_openASC31.P_Values(elecs, :) < 0.05);
    if numel(sig) > 0 % Only plot line if significant indices are found
        line([sig(end), sig(end)], [elecs - 0.5, elecs + 0.5], 'Color', 'red', 'LineWidth', 2);
        sigall(elecs, :) = sig(end);
    else
        sigall(elecs, :) = NaN;
    end
    hold off;
end

% Draw horizontal lines between significant points
for elecs = 1:size(TFCE_openASC31.TFCE_Obs, 1) - 1
    line([sigall(elecs, :), sigall(elecs + 1, :)], [elecs + 0.5, elecs + 0.5], 'Color', 'red', 'LineWidth', 2);
end

% Set axis properties
lim = clim; % Retrieve the color axis limits
axis off;

% Export the figure as a JPG image with a specified resolution
export_fig(join([saveDir 'AAopenclosed_lastordered.jpg'], ''), '-r300')

% Create another figure with a white background
figure,
set(gcf, 'color', 'w')

% Display an image with the CCopenclosed_lastordered.TFCE_Obs data using the custom color map
imagesc(CCopenclosed_lastordered.TFCE_Obs)
clim(lim) % Set color axis limits

% Specify colormap and hide axis
colormap(customMap)
axis off;

% Loop through electrode rows and highlight significant points
for elecs = 1:size(CCopenclosed_lastordered.TFCE_Obs, 1)
    hold on;
    sig = find(CCopenclosed_lastordered.P_Values(elecs, :) < 0.05);
    if numel(sig) > 0 % Only plot line if significant indices are found
        line([sig(end), sig(end)], [elecs - 0.5, elecs + 0.5], 'Color', 'red', 'LineWidth', 2);
        sigall(elecs, :) = sig(end);
    else
        sigall(elecs, :) = NaN;
    end
    hold off;
end

% Draw horizontal lines between significant points
for elecs = 1:size(AAopenclosed_lastordered.TFCE_Obs, 1) - 1
    line([sigall(elecs, :), sigall(elecs + 1, :)], [elecs + 0.5, elecs + 0.5], 'Color', 'red', 'LineWidth', 2);
end

% Export the figure as a JPG image with a specified resolution
export_fig(join([saveDir 'CCopenclosed_lastordered.jpg'], ''), '-r300')


figure,
set(gcf,'color','w')
imagesc(CAclosed_lastordered.TFCE_Obs)
lim1 = clim;
colormap(customMap)
colorbar off
axis off;
hold on;
for elecs = 1:size(CCopenclosed_lastordered.TFCE_Obs,1)
    hold on;
    sig = find(CAclosed_lastordered.P_Values(elecs,:) < 0.05);
    if numel(sig) > 0 % Only plot line if significant indices are found
        line([sig(1), sig(1)], [elecs-0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
        sigall(elecs,:) = sig(1);
    else
        sigall(elecs,:)= NaN;
    end
    hold off;
end
for elecs = 1:size(AAopenclosed_lastordered.TFCE_Obs,1)-1
            line([sigall(elecs,:), sigall(elecs+1,:)], [elecs+0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
end
% hold off;
% cb = colorbar();
% cb.Ruler.Exponent = 3;
export_fig(join([saveDir 'CAclosed_lastordered.jpg'],''), '-r300')

figure,
set(gcf,'color','w')
imagesc(CAopen_lastordered.TFCE_Obs)
clim(lim1)
colormap(customMap)
colorbar off
% cb = colorbar();
% cb.Ruler.Exponent = 3;
axis off;
for elecs = 1:size(AAopenclosed_lastordered.TFCE_Obs,1)
    hold on;
    sig = find(CAopen_lastordered.P_Values(elecs,:) < 0.05);
    if numel(sig) > 0 % Only plot line if significant indices are found
        line([sig(1), sig(1)], [elecs-0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
        sigall(elecs,:) = sig(1);
    else
        sigall(elecs,:)= NaN;
    end
    hold off;
end
for elecs = 1:size(AAopenclosed_lastordered.TFCE_Obs,1)-1
            line([sigall(elecs,:), sigall(elecs+1,:)], [elecs+0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
end
export_fig(join([saveDir 'CAopen_lastordered.jpg'],''), '-r300')


figure,
set(gcf,'color','w')
imagesc(DiffAC_lastordered.TFCE_Obs)
colormap(customMap)
clim(lim1)
axis off;
hold on;
for elecs = 1:size(AAopenclosed_lastordered.TFCE_Obs,1)
    hold on;
    sig = find(DiffAC_lastordered.P_Values(elecs,:) < 0.05);
    if numel(sig) > 0 % Only plot line if significant indices are found
        line([sig(1), sig(1)], [elecs-0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
        sigall(elecs,:) = sig(1);
    else
        sigall(elecs,:)= NaN;
    end
    hold off;
end
for elecs = 1:size(AAopenclosed_lastordered.TFCE_Obs,1)-1
            line([sigall(elecs,:), sigall(elecs+1,:)], [elecs+0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
end
colorbar off
% cb = colorbar();
% cb.Ruler.Exponent = 3;
export_fig(join([saveDir 'DiffAC.jpg'],''), '-r300')


% sample entropy
n = 1026; %number of unique color
base = [0 0 .5;%dark blue
    %     0 0.5 1;% middle blue
    0 1 1;% cyan
    0 1 0;% light green
    1 1 0;% yellow
    1 0.5 0;% orange
    1 0 0]; % red
customMap = interp1(linspace(n,0,size(base,1)), base, fliplr(0:n), 'pchip');

saveDir = 'E:\Research_update\eeg\bigdat_entropy\TFCE_relpow\orderedFile\';
figure,
set(gcf,'color','w')
imagesc(AAopenclosedSE.TFCE_Obs)
colormap(customMap)
colorbar
clim(lim)
export_fig(join([saveDir 'AAopenclosedSE_ordered.jpg'],''), '-r300')

figure,
set(gcf,'color','w')
imagesc(CCopenclosedSE.TFCE_Obs)
clim(lim)
colormap(customMap)
colorbar
export_fig(join([saveDir 'CCopenclosedSE_ordered.jpg'],''), '-r300')
% export_fig(join([saveDir 'colorbar1.jpg'],''), '-r300')


% n = 256; %number of unique color
% base = [0 1 1; 0 1 0]; %blue cyan green
% customMapLow = interp1(linspace(n,0,size(base,1)), base, fliplr(0:n), 'pchip');

figure,
set(gcf,'color','w')
imagesc(CAopenSE.TFCE_Obs)
colormap(customMap)
clim(lim1)
cb = colorbar;
% cb = colorbar('southoutside');
cb.Ruler.Exponent = 3;
export_fig(join([saveDir 'CAopenSE_ordered.jpg'],''), '-r300')
% export_fig(join([saveDir 'colorbarVlow.jpg'],''), '-r300')

figure,
set(gcf,'color','w')
imagesc(CAclosedSE.TFCE_Obs)
clim(lim1)
colormap(customMap)
colorbar
export_fig(join([saveDir 'CAclosedSE_ordered.jpg'],''), '-r300')


%%
% Readin Elecs
fid = fopen('E:\Research_update\eeg\bigdat_entropy\femaleASDplot\elecs1.txt'); %col1 is original elecs
data = textscan(fid,'%s%s%s');
fclose(fid);

for ii = 1:32
    ElecIdx1(ii,:) = find(strcmpi(data{2}{ii},data{1}));
end

for ii = 1:32
    AAopenclosedSE.TFCE_Obs1(ii,:) = AAopenclosedSE.TFCE_Obs(ElecIdx1(ii),:);
end

for ii = 1:32
    ElecIdx2(ii,:) = find(strcmpi(data{3}{ii},data{1}));
end

for ii = 1:32
    AAopenclosedSE.TFCE_Obs2(ii,:) = AAopenclosedSE.TFCE_Obs(ElecIdx2(ii),:);
end

for ii = 1:32
    e_loc1(ii) = e_loc(ElecIdx1(ii),:);
end

%%
saveDir = 'E:\Research_update\eeg\bigdat_entropy\TFCE_relpow\lastordered\figure\';
figure
set(gcf,'color','w')
imagesc(AAopenclosed_lastordered.TFCE_Obs)
colormap(customMap)
colorbar off
lim = clim;
axis off;
hold on;
% Add horizontal lines
yline([9, 22], '-','Color', [0.5 0.5 0.5],  'LineWidth', 2);
% Add vertical lines
xline([19, 38, 53, 73, 87], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);

for elecs = 1:size(AAopenclosed_lastordered.TFCE_Obs,1)
    hold on;
    sig = find(AAopenclosed_lastordered.P_Values(elecs,:) < 0.05);
    sig1 = cellstr(num2str(sig'))';
    all = cellstr(num2str((1:100)'))';
    idx = find(~ismember(all, sig1));

    % Compute the indices where the sequence breaks
    idxi = [1, find(diff(idx)~=1)+1, numel(idx)+1];

    % Extract the start and end values of each group
    result = [idx(idxi(1:end-1))-1', idx(idxi(2:end)-1)+1'];

    if numel(sig) > 0 % Only plot line if significant indices are found
        line([sig(1), sig(1)], [elecs-0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
        line([sig(end), sig(end)], [elecs-0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
    end
        for rr = 1:length(result)
            line([result(rr), result(rr)], [elecs-0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
        end
    hold off;
end


% cb = colorbar();
% cb.Ruler.Exponent = 3;
export_fig(join([saveDir 'AAopenclosed_lastordered.jpg'],''), '-r300')

figure,
set(gcf,'color','w')
imagesc(CCopenclosed_lastordered.TFCE_Obs)
clim(lim)
% lim = clim;
colormap(customMap)
% colorbar
axis off;
hold on;
% Add horizontal lines
yline([9, 22], '-','Color', [0.5 0.5 0.5],  'LineWidth', 2);
% Add vertical lines
% xline([19, 38, 53, 73, 87], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
xline([6, 13], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
for elecs = 1:size(AAopenclosed_lastordered.TFCE_Obs,1)
    hold on;
    sig = find(AAopenclosed_lastordered.P_Values(elecs,:) < 0.05);
    sig1 = cellstr(num2str(sig'))';
    all = cellstr(num2str((1:100)'))';
    idx = find(~ismember(all, sig1));

    % Compute the indices where the sequence breaks
    idxi = [1, find(diff(idx)~=1)+1, numel(idx)+1];

    % Extract the start and end values of each group
    result = [idx(idxi(1:end-1))-1', idx(idxi(2:end)-1)+1'];

    if numel(sig) > 0 % Only plot line if significant indices are found
        line([sig(1), sig(1)], [elecs-0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
        line([sig(end), sig(end)], [elecs-0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
    end
        for rr = 1:length(result)
            line([result(rr), result(rr)], [elecs-0.5, elecs+0.5], 'Color', 'red', 'LineWidth', 2);
        end
    hold off;
end
% cb = colorbar();
% cb.Ruler.Exponent = 3;
export_fig(join([saveDir 'CCopenclosed_lastordered.jpg'],''), '-r300')
%export_fig(join([saveDir 'colorbar1.jpg'],''), '-r300')

figure,
set(gcf,'color','w')
imagesc(ACclosed_lastordered.TFCE_Obs)
lim1 = clim;
colormap(customMap)
colorbar off
axis off;
hold on;
% Add horizontal lines
yline([9, 22], '-','Color', [0.5 0.5 0.5],  'LineWidth', 2);
% Add vertical lines
% xline([19, 38, 53, 73, 87], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
xline([6, 13], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);

% Add shading within the rectangle
x = [0, 18, 20, 0];
y = [0, 0, 32.5, 32.5];
patch(x, y, 'red', 'FaceAlpha', 0.25, 'EdgeColor', 'red')

% Add shading within the rectangle
x = [35, 53, 50 ,41];
y = [0, 0, 32.5, 32.5];
patch(x, y, 'red', 'FaceAlpha', 0.25, 'EdgeColor', 'red')

hold off;
% cb = colorbar();
% cb.Ruler.Exponent = 3;
export_fig(join([saveDir 'ACclosed_lastordered.jpg'],''), '-r300')

figure,
set(gcf,'color','w')
imagesc(ACopen_lastordered.TFCE_Obs)
clim(lim1)
colormap(customMap)
colorbar off
% cb = colorbar();
% cb.Ruler.Exponent = 3;
axis off;
hold on;
% Add horizontal lines
yline([9, 22], '-','Color', [0.5 0.5 0.5],  'LineWidth', 2);
% Add vertical lines
% xline([19, 38, 53, 73, 87], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
xline([6, 13], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
% Add shading within the rectangle
x = [35, 50, 50 ,34];
y = [0, 0, 22, 22];
patch(x, y, 'red', 'FaceAlpha', 0.25, 'EdgeColor', 'red')

hold off;
export_fig(join([saveDir 'ACopen_lastordered.jpg'],''), '-r300')


figure,
set(gcf,'color','w')
imagesc(DiffrelPow.TFCE_Obs)
colormap(customMap)
clim(lim1)
axis off;
hold on;
% Add horizontal lines
yline([9, 22], '-','Color', [0.5 0.5 0.5],  'LineWidth', 2);
% Add vertical lines
% xline([19, 38, 53, 73, 87], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
xline([6, 13], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
% Add shading within the rectangle
x = [0, 15, 12, 0];
y = [0, 0, 32.5, 32.5];
patch(x, y, 'red', 'FaceAlpha', 0.25, 'EdgeColor', 'red')

% Add shading within the rectangle
x = [37, 55, 50 ,41];
y = [0, 0, 32.5, 32.5];
patch(x, y, 'red', 'FaceAlpha', 0.25, 'EdgeColor', 'red')
hold off;
colorbar
cb = colorbar();
cb.Ruler.Exponent = 3;
export_fig(join([saveDir 'colorbarlow.jpg'],''), '-r300')


%% thresholded t value map
% sample entropy
n = 1026; %number of unique color
base = [ 
    0 1 0;% light green
    1 1 0;% yellow
    1 0.5 0;% orange
    1 0 0]; % red
customMap = interp1(linspace(n,0,size(base,1)), base, fliplr(0:n), 'pchip');
saveDir = 'E:\Research_update\eeg\bigdat_entropy\orderedFiles\lastordered\figure\';

cclist = {AAopenclosed_lastordered,CCopenclosed_lastordered};
ccname = {'AAopenclosed_lastordered','CCopenclosed_lastordered'};
for CC = 1:length(cclist)
    MatrixTemp = cclist{CC};
    matrix1 = MatrixTemp.P_Values;
    threshold = 0.05;
    mask = matrix1 < threshold;

    matrix2 = MatrixTemp.TFCE_Obs;
    result = matrix2 .* mask;
    result(result==0)=NaN;

    figure;
    set(gcf, 'color', 'w');

    % Create the heatmap
    imagesc(result);

    % Modify the colormap to show NaN values as white
    customMapWithNaN = customMap;
    nanColor = [1 1 1]; % RGB values for white
    customMapWithNaN = [nanColor; customMapWithNaN];

    colormap(gca, customMapWithNaN);
    caxis([min(result(:)), max(result(:))]);

    % Add colorbar
    colorbar;
    
    % uniform limit 
    lim = [0,9000];
    clim(lim)
    % get off axis
    axis off;
    hold on;
    % Add horizontal lines
    yline([9, 22], '-','Color', [0.5 0.5 0.5],  'LineWidth', 2);
    % Add vertical lines
    xline([19, 38, 53, 73, 87], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
%     xline([6, 13], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
    export_fig(join([saveDir ccname{CC} '.jpg'],''), '-r300')
end


% sample entropy
n = 2056; %number of unique color
base = [ 0 0 .5;%dark blue
    0 0.5 1;% middle blue
    0 1 1;% cyan
    0 1 0;% light green
    1 1 0;% yellow
    1 0.5 0;% orange
    1 0 0]; % red

customMap = interp1(linspace(n,0,size(base,1)), base, fliplr(0:n), 'pchip');

cclist = {TFCE_AAopenclosedOrig,TFCE_CCopenclosedOrig};
ccname = {'TFCE_AAopenclosedOrig','TFCE_CCopenclosedOrig'};
for CC = 1:length(cclist)
    MatrixTemp = cclist{CC};
    matrix1 = MatrixTemp.P_Values;
    threshold = 0.05;
    mask = matrix1 < threshold;

    matrix2 = MatrixTemp.TFCE_Obs;
    result = matrix2 .* mask;
    result(result==0)=NaN;

    figure;
    set(gcf, 'color', 'w');

    % Create the heatmap
    imagesc(result);

    % Modify the colormap to show NaN values as white
    customMapWithNaN = customMap;
    nanColor = [1 1 1]; % RGB values for white
    customMapWithNaN = [nanColor; customMapWithNaN];

    colormap(gca, customMapWithNaN);
%     caxis([min(result(:)), max(result(:))]);

    % Add colorbar
%     colorbar;
    
    % uniform limit 
    lim = [-21000,20000];
    clim(lim)
    % get off axis
    axis off;
    hold on;
    % Add horizontal lines
    yline([9, 22], '-','Color', [0.5 0.5 0.5],  'LineWidth', 2);
    % Add vertical lines
    xline([19, 38, 53, 73, 87], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
%     xline([6, 13], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
    export_fig(join([saveDir ccname{CC} '.jpg'],''), '-r300')
end


% % sample entropy all positive value
% base = [ %0 0 .5;%dark blue
%     %     0 0.5 1;% middle blue
%     0 1 1;% cyan
%     0 1 0;% light green
%     %     1 1 0;% yellow
%     1 0.5 0];% orange
% customMap = interp1(linspace(n,0,size(base,1)), base, fliplr(0:n), 'pchip');

% cclist = {CAopen_lastordered,CAclosed_lastordered,DiffAC_lastordered};
% ccname = {'CAopen_lastordered','CAclosed_lastordered','DiffAC_lastordered'};
n = 1026; %number of unique color
base = [0 1 0;% light green
    1 1 0;% yellow
    1 0.5 0;% orange
    1 0 0]; % red
customMap = interp1(linspace(n,0,size(base,1)), base, fliplr(0:n), 'pchip');
saveDir = 'E:\Research_update\eeg\bigdat_entropy\age\figure\';
% cclist = {TFCE_ASCopen31,TFCE_ASCopen32,TFCE_ASCopen21,TFCE_CONopen31,TFCE_CONopen32,TFCE_CONopen21,...,
%     TFCE_ASCclosed31,TFCE_ASCclosed32,TFCE_ASCclosed21,TFCE_CONclosed31,TFCE_CONclosed32,TFCE_CONclosed21};
% ccname = {'TFCE_ASCopen31','TFCE_ASCopen32','TFCE_ASCopen21','TFCE_CONopen31','TFCE_CONopen32','TFCE_CONopen21',...,
% 'TFCE_ASCclosed31','TFCE_ASCclosed32','TFCE_ASCclosed21','TFCE_CONclosed31','TFCE_CONclosed32','TFCE_CONclosed21'};

cclist = {TFCE_openASC31,TFCE_openASC32,TFCE_openASC21,TFCE_openCON31,TFCE_openCON32,TFCE_openCON21,...,
    TFCE_closedASC31,TFCE_closedASC32,TFCE_closedASC21,TFCE_closedCON31,TFCE_closedCON32,TFCE_closedCON21};
ccname = {'TFCE_openASC31','TFCE_openASC32','TFCE_openASC21','TFCE_openCON31','TFCE_openCON32','TFCE_openCON21',...,
   'TFCE_closedASC31','TFCE_closedASC32','TFCE_closedASC21','TFCE_closedCON31','TFCE_closedCON32','TFCE_closedCON21'};
saveDir = 'E:\Research_update\eeg\bigdat_entropy\age\OrigAge_relpow\Figure\';
for CC = 1:length(cclist)
    MatrixTemp = cclist{CC};
    matrix1 = MatrixTemp.P_Values;
    threshold = 0.05;
    mask = matrix1 < threshold;

    matrix2 = MatrixTemp.TFCE_Obs;
    result = matrix2 .* mask;
    result(result==0)=NaN;

    figure;
    set(gcf, 'color', 'w');

    % Create the heatmap
    imagesc(result);

    % Modify the colormap to show NaN values as white
    customMapWithNaN = customMap;
    nanColor = [1 1 1]; % RGB values for white
    customMapWithNaN = [nanColor; customMapWithNaN];

    colormap(gca, customMapWithNaN);
%     caxis([min(result(:)), max(result(:))]);

    % Add colorbar
    cb = colorbar;

    % Set the font size of the colorbar axis
    cb.FontSize = 12;
%     % Adjust the exponent of the colorbar tick labels
    cb.Ruler.Exponent = 3; % Set the desired exponent value

    % uniform limit 
    lim = [0,3000];
    clim(lim)
    % get off axis
    axis off;
    hold on;
    % Add horizontal lines
    yline([9, 22], '-','Color', [0.5 0.5 0.5],  'LineWidth', 2);
    % Add vertical lines
    xline([19, 38, 53, 73, 87], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
%     xline([6, 13], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
    colorbar off
    export_fig(join([saveDir ccname{CC} '.jpg'],''), '-r300')
%     export_fig(join([saveDir 'colorbarhigh.jpg'],''), '-r300')
end


%%
% cclist = {CAopen_lastordered,CAclosed_lastordered,DiffAC_lastordered};
% ccname = {'CAopen_lastordered','CAclosed_lastordered','DiffAC_lastordered'};
n = 1026; %number of unique color
base = [0 1 0;% light green
    1 1 0;% yellow
    1 0.5 0;% orange
    1 0 0]; % red
customMap = interp1(linspace(n,0,size(base,1)), base, fliplr(0:n), 'pchip');
saveDir = 'E:\Research_update\eeg\bigdat_entropy\age\lastorderedAge\Figure\';

cclist = {Age1CON_ASCclosed,Age1CON_ASCopen,Age2CON_ASCclosed,Age2CON_ASCopen,Age3CON_ASCclosed,Age3CON_ASCopen};
ccname = {'Age1CON_ASCclosed','Age1CON_ASCopen','Age2CON_ASCclosed','Age2CON_ASCopen','Age3CON_ASCclosed','Age3CON_ASCopen'};
for CC = 1:length(cclist)
    MatrixTemp = cclist{CC};
    matrix1 = MatrixTemp.P_Values;
    threshold = 0.05;
    mask = matrix1 < threshold;

    matrix2 = MatrixTemp.TFCE_Obs;
    result = matrix2 .* mask;
    result(result==0)=NaN;

    figure;
    set(gcf, 'color', 'w');

    % Create the heatmap
    imagesc(result);

    % Modify the colormap to show NaN values as white
    customMapWithNaN = customMap;
    nanColor = [1 1 1]; % RGB values for white
    customMapWithNaN = [nanColor; customMapWithNaN];

    colormap(gca, customMapWithNaN);
%     caxis([min(result(:)), max(result(:))]);

    % Add colorbar
    cb = colorbar;

    % Set the font size of the colorbar axis
    cb.FontSize = 12;
%     % Adjust the exponent of the colorbar tick labels
    cb.Ruler.Exponent = 3; % Set the desired exponent value

    % uniform limit 
    lim = [0,3000];
    clim(lim)
    % get off axis
    axis off;
    hold on;
    % Add horizontal lines
    yline([9, 22], '-','Color', [0.5 0.5 0.5],  'LineWidth', 2);
    % Add vertical lines
    xline([6, 13], '--','Color', [0.5 0.5 0.5], 'LineWidth', 2);
    colorbar off
    export_fig(join([saveDir ccname{CC} '.jpg'],''), '-r300')
%     export_fig(join([saveDir 'colorbarhigh.jpg'],''), '-r300')
end
