function [ out ] = berst1( directory, plt, plt_algo, removePlateau )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BeRST 1 Data Analysis Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Sang Min Han, UC Berkeley
% Last modified: 09/13/2022
%
% DESCRIPTION:
% berst1 calculates the fast depolarization events/spiking rate using the
% contrained FOOPSI algorithm inside the CaImAn package
%
% DEPENDENCIES:
% CaImAn-MATLAB (Version 2):
%   https://github.com/flatironinstitute/CaImAn-MATLAB
%
% CVX (Version 2.2, January 2020, Build 1148):
%   http://cvxr.com/cvx/download/
%
% INPUTS:
% 'directory' is the path string of the directory containing the data files
% 'plt' is the binary flag for generating the summary plots
% 'plt_algo' is the binary flag for plotting the algorithm details
% 'removePlateau' is the 0/1 flag for removing extraneous spikes caused by
%   plateaus
%
% OUTPUT:
% 'out' is a cell array of struct arrays. Each struct has the following
% fields:
%   'filePath': path string of the data file
%   'fileName': name string of the data file
%   'times': vector of time stamps in seconds
%   'rawDFF': raw Delta F/F0 fluorescence readings time series vector
%   'dnCalcium': denoised calcium signal vector from constrained_foopsi
%   'dcSpikes': deconvolved spikes vector output from constrained_foopsi
%   'baseline': baseline Delta F/F0 calcium concentration
%   'initC': initial calcium concentraion calculated by constrained_foopsi
%   'timeConstants': discrete AR time constant(s) from constrained_foopsi
%   'noiseSD': noise standard deviation
%   'firingRates': vector of inferred moving average firing rates
%   'spikeRate': overall net firing rate
%   'colNums': data column number (cell number+1)
%   'params': algorithm parameters (pcut_basal, pcut_spikes, last);
%
% There is one struct array (one cell array element) for each subdirectory
% inside 'directory'
% The length of each struct array is the number of files in the
% subdirectory
%
% summary.mat comprised of the out element containing the struct array is
% saved in each respective subdirectory
%
% USAGE:
% out= berst1('/path/to/data/files', 0, 0, 1)
%

%% global initializations

blue= [0, 0.4470, 0.7410]; % default blue color
orange= [0.8500, 0.3250, 0.0980]; % default red
yellow= [0.9290, 0.6940, 0.1250]; % default yellow
green= [0.4660, 0.6740, 0.1880]; % default green

%% parameters

Ts= 21.5e-3; % sampling interval
params= struct('pcut_basal', 0.25, 'pcut_spikes', 1, 'last', []);
% pcut_basal is the top percentile to be used for basal activity
% pcut_spikes is the top percentile to be used for spike activity
% last is the sample end number (leave blank to analyze the entire data)
diffTh= 2; % finite difference threshold
stdTh= 1.2; % standard deviation filter threshold

addpath(genpath(pwd));
dataFiles= dir([directory '/**/*.mat']);
fileNames= extractfield(dataFiles, 'name');
filePaths= extractfield(dataFiles, 'folder');
[diffPaths, ~, ic]= unique(filePaths);

out= cell(1, length(diffPaths));

%% analyze the data files

nBasal= floor(1/Ts)/2; % length of trace for basal activity estimation

% catch and handle data file format exceptions and inconsistencies
ell0= ic(1); pathIdx= 0;
for i= 1:length(fileNames)
    ell= ic(i);

    tables= load([filePaths{i} '/' fileNames{i}]);
    names= fieldnames(tables);
    index= 0;
    fileNamesCmpr= regexprep(fileNames{i}, '-', '');
    fileNamesCmpr= fileNamesCmpr(~isspace(fileNamesCmpr));
    for ii= 1:length(names)
        namesCmpr= regexprep(names{ii}, '-', '');
        namesCmpr= namesCmpr(~isspace(namesCmpr));
        if contains(fileNamesCmpr, namesCmpr)
            index= ii;
        end
    end
    if index == 0
        index= 1;
    end
    tempData= tables.(names{index});
    if istable(tempData)
        data= table2array(tempData);
    else
        data= tempData;
    end
    [M, N]= size(data);

    try
        idx= find(prod(ismissing(data(1:floor(1/Ts),:),'--'),1) == 1);
    catch
        idx= find(prod(ismissing(data(1:floor(1/Ts),:)),1) == 1);
    end

    col= [];
    if isempty(idx)
        try
            col= find(ismissing(data(1,:),'--') == 1);
        catch
            col= find(ismissing(data(1,:)) == 1);
        end
        if isempty(col)
            col= find(ismissing(data(1,:)) == 1);
        end
    end

    if istable(data)
        if isempty(idx) && numel(col) <= 1
            data= table2array(data);
        elseif numel(col) > 1
            row= 1;
            while numel(col) > 1
                row= row+1;
                col= find(ismissing(data(row,:)) == 1);
            end
            data= table2array(data(row:end,:));
        else
            dataL= table2array(data(:,1:idx-1));
            dataR= table2array(data(:,idx+1:N));
            data= [dataL zeros(M,1) dataR];
        end
    end
    data(isnan(data))= 0;

    last= length(data(:,2));
    for m= 1:length(data(:,2))-floor(1/Ts)
        if sum(data(m:m+floor(1/Ts),2)) == 0
            last= m;
            break;
        end
    end

    [val, idx]= min(sum(abs(data)));
    if val ~= 0
        idx= 1;
    end

    % local initializations
    times= data(:,1)./1e3; % in seconds
    deltaT= times(2)-times(1);
    if times(last) == 0
        idx0= find(times == 0, 2);
        last= idx0(end)-floor(1/Ts); % floor(0.99*idx0(end));
    end
    params.last= last;

    if plt
        close all;
        figure('units', 'normalized', 'outerposition', [0.2 0.1 0.6 0.8]);
    end

    % fast depolarization events/spikes inference
    DFF= zeros(last, N-idx);
    c= zeros(last, N-idx);
    s= zeros(last, N-idx);
    g= zeros(10, N-idx);
    firingRates= zeros(last, N-idx);
    netRate= zeros(1, N-idx);
    colNums= zeros(1, N-idx);
    k= 1;
    for j= idx+1:N
        y= data(:,j);

        %%% high-pass filter %%%
        % h= hpf(params.w_c, params.N_taps);
        % y_filtered= highpass(y, params.w_c); % conv(y, h, 'same');
        % DFF(:,k)= y_filtered(1:last);

        DFF(:,k)= y(1:last);

        colNums(k)= j;

        options= struct();
        % struct('noise_method', 'median');
        % struct('p',1);
        % struct('method', 'dual');

        % basal
        stds_basal= movstd(DFF(:,k), [0 nBasal-1], 'omitnan');
        sorted_basal= sort(stds_basal, 'ascend');
        cutoff_basal= floor(params.pcut_basal*length(sorted_basal));
        std_basal= sorted_basal(floor(cutoff_basal/2));
        basalStart= min(find(stds_basal == std_basal), length(DFF(:,k))-...
            nBasal);
        basalEnd= basalStart+nBasal-1;
        basal= DFF(basalStart:basalEnd,k);
        mag_basal= floor(log10(std_basal));
        b= mean(basal);

        % nonbasal
        nonbasal= [DFF(1:basalStart-1,k); DFF(basalEnd+1:end,k)];
        std_nonbasal= std(nonbasal);
        mag_nonbasal= floor(log10(std_nonbasal));

        sn= std_basal; % noise standard deviation

        % run constrained_foopsi
        % [c, b, cin, g, sn, sp]= constrained_foopsi(DFF(:,k), b, c1, g, ...
        %     sn, options);
        [c(:,k), b, c1, g_k, sn, s(:,k)]= constrained_foopsi(DFF(:,k), ...
            b, [], [], sn, options);

        g(1:length(g_k),k)= g_k; % discrete AR time constants

        spikes= s(:,k);
        spikes= spikes(spikes ~= 0);
        sorted_spikes= sort(spikes, 'descend');
        cutoff_spikes= floor(params.pcut_spikes*length(spikes));
        mean_spikes= mean(sorted_spikes(1:cutoff_spikes));
        std_spikes= std(sorted_spikes(1:cutoff_spikes));

        potentialPlateau= 0;
        if abs(std_nonbasal-std_basal)/std_basal >= 0.75
            thresh= max(0, b) + 1.5*std_basal;
            fprintf('threshold is mean_spikes+1.5*std_spikes\n');
            potentialPlateau= 1;
        else
            thresh= max(0, b) + 4*std_basal;
            fprintf('threshold is b+4*std_basal\n');
        end
        sp= s(:,k).*(s(:,k) > thresh);

        % plateau
        if removePlateau && potentialPlateau
            finiteDiff= [diff(DFF(:,k)); 0];
            stdDFF= stdfilt(DFF(:,k));
            plateau= (abs(finiteDiff) < diffTh) & (stdDFF < stdTh) & ...
                (DFF(:,k) >= mean(DFF(:,k)) + std(DFF(:,k)));
            sp(find(plateau == 1))= 0;
        end

        T= times(last); % in seconds

        % rate calculation
        netRate(k)= sum(sp > 0)/T;

        firingRates(:,k)= movsum(sp > 0, [floor(1/Ts)-1 0], 'omitnan')/...
            (deltaT*(floor(1/Ts)-1));
        fprintf('Neuron %i/%i:\t sn = %f\t Net Rate = %f\n', j-1, N-1, ...
            sn, netRate(k));

        if plt % plot
            clf;

            subplot(2,1,1) %%
            hold on;
            plot(times(1:last), DFF(:,k), 'Color', green, 'LineWidth', 1.2)
            plot(times(basalStart:basalEnd), ...
                DFF(basalStart:basalEnd,k), 'k', 'LineWidth', 1.2)
            plot(times(1:last), s(:,k), 'Color', orange, 'LineWidth', 1.2)
            plot(times(1:last), sn*ones(1,last), 'k--', 'LineWidth', 1)
            plot(times(1:last), std_nonbasal*ones(1,last), 'Color', ...
                green, 'LineStyle', '--', 'LineWidth', 1)
            plot(times(1:last), thresh*ones(1,last), 'Color', yellow, ...
                'LineStyle', '--', 'LineWidth', 1)
            hold off;
            xlim([times(1) times(last)]);
            legend({'$F$', 'Basal $F$', 'Spikes', 'Basal $\sigma$', ...
                'Nonbasal $\sigma$', 'Threshold'}, 'Interpreter', 'LaTeX');
            ttl= ['File ' fileNames{i}];
            ttl= strrep(ttl, '_', '\_');
            title(ttl);

            subplot(2,1,2) %%
            plot(times(1:last), sp, 'Color', blue, 'LineWidth', 1.2);
            xlim([times(1) times(last)]);
            title(sprintf('Neuron %i/%i: Net Rate = %.4f', j-1, N-1, ...
                netRate(k)), 'Interpreter', 'LaTeX');

            waitforbuttonpress;
        end

        k= k+1;
    end

    if ell == ell0
        pathIdx= pathIdx+1;
    else
        pathIdx= 1;
    end

    ell0= ell;

    % out cell array of struct arrays
    out{1,ic(i)}(pathIdx)= struct('filePath', filePaths{i}, 'fileName', ...
        fileNames{i}, 'times', times, 'rawDFF', DFF, 'dnCalcium', c, ...
        'dcSpikes', s, 'baseline', b, 'initC', c1, 'timeConstants', g, ...
        'noiseSD', sn, 'firingRates', firingRates, 'spikeRate', ...
        netRate, 'colNums', colNums, 'params', params);

    if plt_algo % plot algorithm detail
        if N-idx <= 3
            R= 1;
            C= N-idx;
        else
            R= floor((N-idx)/2);
            C= ceil((N-idx)/2);
        end

        f1= figure; clf; %%
        for k= 1:N-idx
            subplot(R,C,k) %
            plot(times(1:last), DFF(:,k), 'LineWidth', 1);
            xlim([times(1) times(last)]);
            xlabel('Time (s)', 'Interpreter', 'LaTeX');
            ylabel('$\frac{\Delta F}{F_0}$', 'Interpreter', 'LaTeX');
            title(sprintf('Column %i', colNums(k)), 'Interpreter', ...
                'LaTeX');
        end
        set(gcf, 'Position', get(0, 'Screensize'));

        f2= figure; clf; %%
        for k= 1:N-idx
            subplot(R,C,k) %
            plot(times(1:last), c(:,k), 'LineWidth', 1);
            xlim([times(1) times(last)]);
            xlabel('Time (s)', 'Interpreter', 'LaTeX');
            ylabel('$\frac{\Delta F}{F_0}$', 'Interpreter', 'LaTeX');
            title(sprintf('Column %i', colNums(k)), 'Interpreter', ...
                'LaTeX');
        end
        set(gcf, 'Position', get(0, 'Screensize'));

        f3= figure; clf; %%
        for k= 1:N-idx
            subplot(R,C,k) %
            plot(times(1:last), s(:,k), 'LineWidth', 1);
            xlim([times(1) times(last)]);
            xlabel('Time (s)', 'Interpreter', 'LaTeX');
            ylabel('$\frac{\Delta F}{F_0}$', 'Interpreter', 'LaTeX');
            title(sprintf('Column %i', colNums(k)), 'Interpreter', ...
                'LaTeX');
        end
        set(gcf, 'Position', get(0, 'Screensize'));

        % save figure
        exportgraphics(f1, filePaths{i}, [fileNames{i} '_raw'], '.eps', ...
            'BackgroundColor', 'none', 'ContentType', 'vector');
        exportgraphics(f2, filePaths{i}, [fileNames{i} '_dn'], '.eps', ...
            'BackgroundColor', 'none', 'ContentType', 'vector');
        exportgraphics(f3, filePaths{i}, [fileNames{i} '_dc'], '.eps', ...
            'BackgroundColor', 'none', 'ContentType', 'vector');
        close all;
    end
    fprintf('file %s (%i/%i) complete...\n', fileNames{i}, i, ...
        length(fileNames));
end
fprintf('finished.\n');

%% save workspace

for i= 1:length(diffPaths)
    match= 0;
    for j= 1:length(out)
        path= extractfield(out{j}, 'filePath');
        path= path{1};
        if strcmp(diffPaths{i}, path)
            match= 1;
            break;
        end
    end
    if match
        summary= out{j};
        save([diffPaths{i} '/summary.mat'], 'summary', '-v7.3');
    end
end
fprintf('saved.\n');

end

