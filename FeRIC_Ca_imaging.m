% %Ca2+ imaging analysis
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cytosolic Ca2+ levels were monitored in cells expressing GCaMP6 
%     and TRP ion channels
% GCaMP6 fluorescence was acquired at 1 Hz.
% 
% Read NIfTI-1 experiment time series data from working directory
% Data Format: NIfTI-1. 
% 
% Basal Session (F0): 1-61 s;
% RF Session/ Stimulus session: 62-301 s; 
% Agonist Session: 302-421 s.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Working Directory Example (contains 2 or more experiments under 
%                             same experimental condition): 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                          
%       pwd/
%           experiment_1.hdr
%           experiment_1.img
%           experiment_2.hdr
%           experiment_2.img
%           ...
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
% OUTPUTs:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	experiment_1.mat --
%           'cell_data_corrected': 
%                 Photobleaching-corrected time-series data;
% 
%           'ratio': 
%                 Responsiveness = (number of RF- or agonist- responsive cells) 
%                                 / (number of the total analyzed cells);
% 
%           'rfmask': 
%                 The binary mask to identify RF- or agonist-responsive cells;
% 
%           'ncell' : 
%                 The number of detected cells;
% 
%           'elim_cell': 
%                 The number of cells excluded from the analysis (when fitted
%                 baseline has growth factor abs(b)> 0.002);
% 
%           'elim_cell_data:': 
%                 The time-series profile of GCaMP6 deltaF/F0 of cells 
%                 excluded from the analysis;
% 
%           'L': 
%                 imshow(label2rgb(reshape(L,[image_height, image_width]))) 
%                 will show the cell segmentation image;
% 
%           'rftrsh': 
%                 The threshold value to stablish if a cell is RF- or 
%                 agonist- responsive (default 10 std of basal session, F0);
% 
%           'std_basal': 
%                 The standard deviation of the basal session(F0) for 
%                 each cell;
% 
%           'sig_ave_rf': 
%                 The averaged time-series profile of GCaMP6 deltaF/F0 of 
%                 RF- or agonist-responsive cells;
% 
%           'sig_ave_noresp': 
%                 The averaged time-series profile of GCaMP6 deltaF/F0 of 
%                 non-responsive cells;
% 
%           'sig_ave_all': 
%                 The averaged time-series profile of GCaMP6 deltaF/F0 of 
%                 all analyzed cells;
% 
%           'sig_sem_rf':
%                 The standard error of the mean(SEM)of the time-series 
%                 GCaMP6 deltaF/F0 of the RF- or agonist 
%                 responsive cells;
% 
%           'sig_sem_noresp': 
%                 The SEM of the time-series of GCaMP6 deltaF/F0 of 
%                 non-responsive cells;
% 
%           'sig_sem_all': 
%                 The SEM of the time-series of GCaMP6 deltaF/F0 of 
%                 all analyzed cells;
% 
%           'rfauc': 
%                 The averaged area under the curve over time of the GCaMP6 
%                 deltaF/F0 during RF or agonist session (62-301 s) for all 
%                 analyzed cells;
% 
% 	experiment_1_all.fig -- 
%                 The figure of the averaged time-series profile of GCaMP6 
%                 deltaF/F0(±SEM) of all analyzed cells;
% 
% 	experiment_1_rf.fig -- 
%                 The figure of the averaged time-series profile of GCaMP6 
%                 deltaF/F0(±SEM) of RF- or agonist-responsive cells;
%     
% 	experiment_1_noresp.fig -- 
%                 The figure of the averaged time-series profile of GCaMP6 
%                 deltaF/F0(±SEM) of non-responsive cells;
%     
% 	experiment_2.mat
% 	experiment_2_all.fig
% 	experiment_2_rf.fig
% 	experiment_2_noresp.fig
%     ...
% 
%     experiment.mat --  
%           'AUC_rf': 
%                 The averaged area under the curve (AUC) for all responsive 
%                 cells in the working directory;
% 
%           'AUC_all': 
%                 The averaged area under the curve (AUC) for all analyzed 
%                 cells in the working directory;
% 
%           'exp_num': 
%                 The number of independent experiments analyzed in the 
%                 working directory;
% 
%           'actNum': 
%                 The number of RF- or agonist-responsive cells from all 
%                 independent experiments in the working directory;
% 
%           'allNum': 
%                 The total number of cells analyzed from all experiments 
%                 in the working directory;
% 
%           'ratio_tol': 
%                 The percent of RF- or agonist-responsiveness from all 
%                 experiments in the working directory;
% 
%           'ave_coll_rf': 
%                 The averaged time-series profile of the GCaMP6 
%                 deltaF/F0(±SEM) of RF- or agonist-responsive cells in the 
%                 working directory;
% 
%           'ave_coll_noresp': 
%                 The averaged time-series profile of the GCaMP6 
%                 deltaF/F0(±SEM) of non-responsive cells in the working 
%                 directory;
% 
%           'ave_coll_all': 
%                 The averaged time-series profile of the GCaMP6 
%                 deltaF/F0(±SEM) of analyzed cells in the working directory;
% 
%           'sem_coll_rf': 
%                 The SEM of the time-series profile of the GCaMP6 deltaF/F0 
%                 of RF- or agonist-responsive cells in the working directory;
% 
%           'sem_coll_noresp': 
%                 The SEM of the time-series profile of the GCaMP6 deltaF/F0 
%                 of non-responsive cells in the working directory;
% 
%           'sem_coll_all': 
%                 The SEM of the time-series profile of the GCaMP6 deltaF/F0 
%                 of all analyzed cells in the working directory;
% 
%           'max_value': 
%                 The averaged maximum value of the GCaMP6 deltaF/F0 for the 
%                 RF or agonist session of all analyzed cells in the working 
%                 directory;
% 
%           'tail_value': 
%                 The averaged value at the end (average over 290-295 s) for 
%                 the RF or agonist session of all analyzed cells in the 
%                 working directory;
% 
% 
%     experiment.fig --  
%             The figure of the averaged time-series profile of the GCaMP6 
%             deltaF/F0(±SEM)of analyzed cells from all experiments conducted 
%             under the same experimental conditions in the working directory.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
% Area under the curve (AUC) = sum(time_series_value);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Last Updated: Dec 15 2019, J Chen @ UC Berkeley.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lst = dir(pwd);
for indf = 1:length(lst)
    fname = lst(indf).name;
    if contains(fname,'.hdr')
        %% load data
        % ~27 sec
        filename = fname(1:length(fname)-4);
        fprintf('%s ', filename);
        data = double(squeeze(niftiread([filename,'.hdr'])));
        sz = size(data);
        %% cell parcellation
        [L,mask] = cell_parc(data);
        ncell = length(unique(L))-2;
        elim_cell = [];
        data = reshape(data,[],sz(3));
        L = reshape(L,[],1);
        for indc = 2:max(L(:))
            cell_data(indc-1,:) = mean(data(L==indc,:),1);
        end
        indn = 0;
        for indc = 1:size(cell_data,1)
            sig = squeeze(cell_data(indc,:));
            [sig,err] = j_baseline(sig,6:60);
            if err > 0
                elim_cell = [elim_cell,indc];
            else
                indn = indn + 1;
                cell_data_corrected(indn,:) = sig;
            end
        end
        h = figure;
        subplot(2,2,1)
        for indc = 1:size(cell_data_corrected,1)
            plot(cell_data_corrected(indc,:))
            hold on
        end
        ylabel('fluorescence Intensity');
        xlabel('Time (sec)')
        title([filename, ' Time Profile of Each Cell'])
        hold off
        %% area under the curve (normalized by time)
        rfauc = sum(cell_data_corrected(:,62:301),2)/length(62:301); % change if experiment rf session changes
        %rftrsh = 0.1;
        rftrsh = 10*std(squeeze(cell_data_corrected(:,6:61)),0,2); % can change
        resp_data_cell = cell_data_corrected(rfauc >= rftrsh,:);
        nonresp_data_cell = cell_data_corrected(rfauc < rftrsh,:);
        elim_data_cell = cell_data(elim_cell,:);
        
        ratio = size(resp_data_cell,1)/(ncell-length(elim_cell));
        fprintf('ratio = %f \n',ratio)
        %%  averaged curve
        sig_ave_noresp = mean(nonresp_data_cell,1);
        sig_ave_all =  mean([resp_data_cell;nonresp_data_cell],1);
        sig_ave_rf = mean(resp_data_cell,1);
        
        sig_sem_noresp = std(nonresp_data_cell,0,1)/sqrt(size(nonresp_data_cell,1));
        sig_sem_all =  std(cell_data_corrected,0,1)/sqrt(size(cell_data_corrected,1));
        sig_sem_rf = std(resp_data_cell,0,1)/sqrt(size(resp_data_cell,1));
        
        subplot(2,2,3)
        errorbar(sig_ave_rf,sig_sem_rf)
        hold on
        title(['rf responsive', ' (ratio = ', num2str(ratio), ')'])
        ylabel('fluorescence Intensity');
        xlabel('Time (sec)')
        subplot(2,2,4)
        errorbar(sig_ave_noresp,sig_sem_noresp)
        title(['no responsive',' (ratio = ', num2str(ratio), ')'])
        ylabel('fluorescence Intensity');
        xlabel('Time (sec)')
        subplot(2,2,2)
        errorbar(sig_ave_all,sig_sem_all)
        title(['all', ' (ratio = ', num2str(ratio), ')'])
        ylabel('fluorescence Intensity');
        xlabel('Time (sec)')
        
        saveas(h, [filename '.fig']);
        
        save([filename '.mat'], 'filename', 'ratio', 'sig_ave_rf', 'sig_ave_noresp', 'sig_ave_all', ...
            'nonresp_data_cell','resp_data_cell','elim_cell','elim_data_cell'...
            ,'sig_sem_noresp','sig_sem_all','sig_sem_rf'...
            ,'mask', 'L' ,'rftrsh','cell_data_corrected', 'cell_data','ncell','rfauc','-v7.3');
        clear('filename', 'ratio', 'sig_ave_rf', 'sig_ave_noresp', 'sig_ave_all', ...
            'nonresp_data_cell','resp_data_cell','elim_cell','elim_data_cell'...
            ,'sig_sem_noresp','sig_sem_all','sig_sem_rf'...
            ,'mask', 'L' ,'rftrsh','cell_data_corrected', 'cell_data','ncell','rfauc')
    end
end

%Collectiong all the experiment data
% searching thru the folder of data of single experiment condition
% Last updated: Dec 15 2019, Jingjia Chen

clear
pathnow = pwd;
loc = find(pathnow == '\',1,'last');
exp_name = pathnow(loc + 1 : end);
lst = dir(pathnow);

sig_coll_rf  = [];
sig_coll_noresp  = [];
actNum = 0;
allNum = 0;
norespNum = 0;
expNum = 0;
ratio_coll = [];
fprintf('%s \n', exp_name)
for ind = 1:length(lst)
    fname = lst(ind).name;
    if ~contains(fname,'.mat')
        continue;
    end
    if isequal(fname,[exp_name,'.mat'])
        continue;
    end
    load(fname,'nonresp_data_cell','resp_data_cell','ratio');
    norespcellNum = size(nonresp_data_cell,1);
    respcellNum = size(resp_data_cell,1);
    allcellNum = norespcellNum + respcellNum;
    ratio_coll = [ratio_coll, ratio];
    sig_coll_rf = [sig_coll_rf;resp_data_cell];
    sig_coll_noresp = [sig_coll_noresp;nonresp_data_cell];
    allNum = allNum + allcellNum;
    actNum = actNum + respcellNum;
    expNum = expNum + 1;
end

sem_coll_rf = std(sig_coll_rf,0,1) ./ sqrt(size(sig_coll_rf,1));
ave_coll_rf = mean(sig_coll_rf,1);

sem_coll_noresp = std(sig_coll_noresp,0,1) ./ sqrt(size(sig_coll_noresp,1));
ave_coll_noresp = mean(sig_coll_noresp,1);

sig_coll_all = [sig_coll_rf;sig_coll_noresp];

sem_coll_all = std(sig_coll_all,0,1) ./ sqrt(size(sig_coll_all,1));
ave_coll_all = mean(sig_coll_all,1);

ratio_tol = actNum/allNum;
sem_ratio = std(ratio_coll)./sqrt(expNum);

h = figure;
subplot(1,2,1)
errorbar(ave_coll_rf,sem_coll_rf)
hold on
errorbar(ave_coll_noresp,sem_coll_noresp)
title([exp_name, ' (ratio = ', num2str(ratio_tol),')'])
legend('rf responsive','no responsive')
subplot(1,2,2)
errorbar(ave_coll_all,sem_coll_all)
title([exp_name, ' all (ratio = ', num2str(ratio_tol),')'])

saveas(h, [exp_name, '.fig']);

fprintf('ratio in total = %f \n',ratio_tol)
fprintf('SEM of ratio = %f \n',sem_ratio)

AUC_rf = sum(ave_coll_rf(:,62:301)); % change if experiment rf session changes
AUC_all = sum(ave_coll_all(:,62:301)); % change if experiment rf session changes
max_value = max(ave_coll_rf(:,62:301)); % change if experiment rf session changes
tail_value = mean(ave_coll_rf(:,290:295));

save(exp_name,'ratio_tol', 'sem_ratio', 'expNum', ...
    'actNum','ave_coll_rf','sem_coll_rf',...
    'allNum','ave_coll_all','sem_coll_all',...
    'ave_coll_noresp','sem_coll_noresp',...
    'AUC_rf','AUC_all','tail_value','max_value')


function [sigcorr,err_flag] = j_baseline(sig, tspan)
sig = reshape(sig,[],1);
tspan = reshape(tspan,[],1);
err_flag = -1;
pe = fit(tspan,sig(tspan),'exp1');
incd = [1: length(sig)]';

bslne = pe.a * exp(pe.b * incd);
if (abs(pe.b) > 0.002)
    sigcorr = (sig - mean(sig(tspan)))./mean(sig(tspan));
    err_flag = 1;
else
    sigcorr = (sig - bslne)./pe.a;
end
end

function [L_out,mask_out] = cell_parc(data_in)
img4mask = max(data_in,[],3);
img4mask = img4mask - min(img4mask(:));
img4mask = img4mask/max(img4mask(:));
masktrs = 0.1;
mask_out = img4mask > masktrs;
mask_out = imopen(mask_out,ones(5,5));
mask_out = imerode(mask_out,strel('disk',5));
mask_perim = bwperim(mask_out,8);

img4mask_e = imerode(img4mask,strel('disk',5));
img4mask_obr = imreconstruct(img4mask_e,img4mask);
img4mask_obrd = imdilate(img4mask_obr,strel('disk',5));
img4mask_obrcbr = imreconstruct(imcomplement(img4mask_obrd),imcomplement(img4mask_obr));
img4mask_obrcbr = imcomplement(img4mask_obrcbr);
mask_em = imregionalmax(img4mask_obrcbr,8);
overlay2 = imoverlay(img4mask,mask_perim | mask_em,'green');

img4mask_c = imcomplement(img4mask);
mask_mod = imimposemin(img4mask_c, ~mask_out | mask_em);
L_out = watershed(mask_mod);

figure
subplot(1,2,1)
imshow(overlay2)
subplot(1,2,2)
imshow(label2rgb(L_out))
end