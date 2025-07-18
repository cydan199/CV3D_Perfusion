function cov_visualization(app)
app.CancelFlag = false;

outPut(app, ["Final CoV projection images are saved in: ", app.save_comp]);

% initialize variables
all_cov_mean = cell(length(app.timepoints),1);
all_cov_median = cell(length(app.timepoints),1);

all_cov_fromV = cell(length(app.timepoints),1);
all_cov_fromZ = cell(length(app.timepoints),1);

for m = 1:length(app.timepoints)
    app.tp = app.timepoints{m};
    cm = getcolormap_sup();
    % create all saving folders
    app.savefns = fullfile(app.save_folder, app.tp);%,'Extraction');
    app.save_cov = fullfile(app.savefns,'CoV');
    app.save_mat = fullfile(app.savefns,'mat');
    app.save_mat_proj = fullfile(app.save_mat,app.proj_layer);
    app.save_cov_proj = fullfile(app.save_cov,app.proj_layer);

    % save_csv = fullfile(savefns,'csv');
    % if ~exist(save_csv)
    %     mkdir(save_csv)
    % end 

    save_nifti = fullfile(app.savefns,'nifti');

    noiseCorr_offset = importdata(fullfile(app.save_mat,'offset.mat'));

    % save_projOCTA_fn = fullfile(save_cov_proj,strcat('Binarized_avgOCTA_proj_offset_',strrep(num2str(noiseCorr_offset),'.','-'),'.mat'));
    save_covMedian_fn = fullfile(app.save_mat_proj,strcat('CoV_proj_fromBscans_median_offset_',num2str(noiseCorr_offset),'.mat'));
    
    save_covMean_fn = fullfile(app.save_mat_proj,strcat('CoV_proj_fromBscans_mean_offset_',num2str(noiseCorr_offset),'.mat'));
    
    save_mask_fn = fullfile(save_nifti,'avg_mask_all_corr.nii');
    if ~exist(save_mask_fn,'file')
        save_mask_fn = fullfile(save_nifti,'avg_mask_all.nii');
    end

    savefn_mat_V = fullfile(app.save_mat_proj,strcat(app.tp,'_OCTA_fromVolume_CoV','.mat'));
    savefn_mat_Z = fullfile(app.save_mat_proj,strcat(app.tp,'_OCTA_fromZeiss_CoV','.mat'));


    %%% read CoV map and draw historgam
    avg_mask_all = double(niftiread(save_mask_fn));
    covBscan_mean = importdata(save_covMean_fn);
    covBscan_median = importdata(save_covMedian_fn);

    cov_fromV = importdata(savefn_mat_V);
    cov_fromZ = importdata(savefn_mat_Z);

    covBscan_mean(isnan(covBscan_mean)) = 0;
    covBscan_median(isnan(covBscan_median)) = 0;

    cov_fromV(isnan(cov_fromV)) = 0;
    cov_fromZ(isnan(cov_fromZ)) = 0;

    all_cov_mean{m} = mat2gray(covBscan_mean.*avg_mask_all);
    all_cov_median{m} = mat2gray(covBscan_median.*avg_mask_all);

    all_cov_fromV{m} =  mat2gray(cov_fromV.*avg_mask_all);
    all_cov_fromZ{m} =  mat2gray(cov_fromZ.*avg_mask_all);

    %%% choose the threshold value based on histogram
    covBscan_MED_mask = covBscan_median.*avg_mask_all;

    cov_fromV_mask = cov_fromV.*avg_mask_all;
    cov_fromZ_mask = cov_fromZ.*avg_mask_all;

    %%% compute mean and std
    measurement = struct;
    measurement.mu_customV = mean(nonzeros(covBscan_MED_mask));
    measurement.std_customV = std(nonzeros(covBscan_MED_mask));

    measurement.mu_customE = mean(nonzeros(cov_fromV_mask));
    measurement.std_customE = std(nonzeros(cov_fromV_mask));

    measurement.mu_zeissE = mean(nonzeros(cov_fromZ_mask));
    measurement.std_zeissE = std(nonzeros(cov_fromZ_mask));

    th = struct;

    %%% Plots & threshold
    %%% CoV from Zeiss
    close all
    figure('Visible','off');
    histogram(nonzeros(cov_fromZ_mask),'NumBins',150,'Normalization','probability');
    title('CoV values from instrument-processed data','FontName','Times New Roman','FontSize',12); 
    xlabel('Pixel Values [a.u.]','FontName','Times New Roman')
    ylabel('Normalized # of Pixels','FontName','Times New Roman')
    xlim([0,1])
    fZ_hist = gca;
    xt = get(fZ_hist, 'YTick');
    set(fZ_hist, 'YTick', xt, 'YTickLabel', round(xt/max(xt),1))
    savefm_hist =  fullfile(app.save_comp, strcat(app.tp,'_CoV_fromZeiss_histogram.png'));
    exportgraphics(fZ_hist,savefm_hist)
    close all

    ptitle = 'CoV values from PlexElite-generated data';
    q = strcat(app.tp," - what is the threshold value: ");
    threshold_Z = inputDialogWithPlot(cov_fromZ_mask, ptitle, q); %input("What is the threshold value: ");
    outPut(app,['Threshold value (Zeiss_E) = ', num2str(threshold_Z)])

    cov_fromZ_mask(cov_fromZ_mask > threshold_Z) = threshold_Z;

    th.threshold_Z = threshold_Z;

    fZ = figure('Visible','off');
    imshow(ind2rgb(round(mat2gray(cov_fromZ_mask).*255),cm),cm);
    savefn_img = fullfile(app.save_comp,strcat(app.tp,'_CoV_fromZeiss_threshold_',strrep(num2str(threshold_Z),'.','-'),'.png'));
    exportgraphics(fZ,savefn_img)
    close all

    % CoV from volume
    % close all
    figure('Visible','off');
    histogram(nonzeros(cov_fromV_mask),'NumBins',150,'Normalization','probability');% 'Normalization','count'
    title('CoV values from volume-extracted data','FontName','Times New Roman','FontSize',12);% 
    ylabel('Normalized # of Pixels','FontName','Times New Roman');xlabel('CoV Value [a.u.]','FontName','Times New Roman')
    xlim([0,1])
    fV_hist = gca;
    xt = get(fV_hist, 'YTick');
    set(fV_hist, 'YTick', xt, 'YTickLabel', round(xt/max(xt),1))
    savefm_hist = fullfile(app.save_comp, strcat(app.tp,'_CoV_fromVolume_histogram.png'));
    exportgraphics(fV_hist,savefm_hist)
    close all

    ptitle = 'CoV values from volume-extracted data';
  
    threshold_V = inputDialogWithPlot(cov_fromV_mask, ptitle, q); %input("What is the threshold value: ");

    outPut(app,['Threshold value (Custom_E) = ', num2str(threshold_V)])
    cov_fromV_mask(cov_fromV_mask > threshold_V) = threshold_V;

    th.threshold_V = threshold_V;

    fV = figure('Visible','off');
    imshow(ind2rgb(round(mat2gray(cov_fromV_mask).*255),cm),cm);
    savefn_img = fullfile(app.save_comp,strcat(app.tp,'_CoV_fromVolume_threshold_',strrep(num2str(threshold_V),'.','-'),'.png'));
    exportgraphics(fV,savefn_img)

    %%% Median CoV from Bscans
    figure('Visible','off');
    histogram(nonzeros(covBscan_MED_mask),'NumBins',150,'Normalization','probability');
    title('CoV values from B-scan basis CoV computation','FontName','Times New Roman','FontSize',12);%'
    xlabel('Pixel Values [a.u.]','FontName','Times New Roman')
    ylabel('Normalized # of Pixels','FontName','Times New Roman')
    xlim([0,1])
    fmedian_hist = gca;
    xt = get(fmedian_hist, 'YTick');
    set(fmedian_hist, 'YTick', xt, 'YTickLabel', round(xt/max(xt),1))
    savefm_hist =  fullfile(app.save_comp, strcat(app.tp,'_medianCoV_histogram.png'));
    exportgraphics(fmedian_hist,savefm_hist)
    close all

    ptitle = 'CoV values from B-scan basis CoV computation';
    threshold_Median = inputDialogWithPlot(covBscan_MED_mask, ptitle, q);% input("What is the threshold value: ");
    outPut(app,['Threshold value (Custom_V) = ', num2str(threshold_Median)])

    
    covBscan_MED_mask(covBscan_MED_mask > threshold_Median) = threshold_Median; % covBscan_MED_mask_norm

    th.threshold_Median = threshold_Median;

    fmedian = figure('Visible','off');
    savefn_img = fullfile(app.save_comp,strcat(app.tp,'_CoV_proj_fromVolume_median_threshold_',strrep(num2str(threshold_Median),'.','-'),'.png'));
    imshow(ind2rgb(round(mat2gray(covBscan_MED_mask).*255),cm),cm);  
    exportgraphics(fmedian,savefn_img)
    close all

    %%% differenciate vessels at different layers
    covBscan_MED_sup = mat2gray(covBscan_MED_mask);

    q = strcat(app.tp," - what is the threshold value to separate the deep vessels: ");
    deep_MED_threshold = findTurningPointLognorm(app, covBscan_MED_mask, ptitle, q);% inputdlg("What is the threshold value for deep vessels: ");
    outPut(app,['Threshold value to separate the deep vessel = ', num2str(deep_MED_threshold)]);
    th.deep_MED_threshold = deep_MED_threshold;

    covBscan_MED_sup(covBscan_MED_sup > deep_MED_threshold) = 0; % (deep_MED_threshold/threshold_Median) if normalized

    th.deep_MED_threshold = deep_MED_threshold;
    measurement.mu_customeV_sup = mean(nonzeros(covBscan_MED_sup));
    measurement.std_customeV_sup = std(nonzeros(covBscan_MED_sup));


    fmedian_sup = figure('Visible','off');
    imshow(ind2rgb(round((covBscan_MED_sup).*255),cm),cm);
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.8;
    savefn_img = fullfile(app.save_comp,strcat(app.tp,'_supVessel_CoV_proj_fromVolume_median_threshold_',strrep(num2str(deep_MED_threshold),'.','-'),'.png'));
    exportgraphics(fmedian_sup,savefn_img)
    close all

    covBscan_MED_deep = covBscan_MED_mask;
    covBscan_MED_deep(covBscan_MED_deep <= deep_MED_threshold) = 0; % (deep_MED_threshold/threshold_Median) if mat2gray

    measurement.mu_customeV_deep = mean(nonzeros(covBscan_MED_deep));
    measurement.std_customeV_deep = std(nonzeros(covBscan_MED_deep));

    fmedian_deep = figure('Visible','off');
    imshow(ind2rgb(round((covBscan_MED_deep).*255),cm),cm); % mat2gray
    ax = gca;
    ax.TitleFontSizeMultiplier = 0.8;
    savefn_img = fullfile(app.save_comp,strcat(app.tp,'_deepVessel_CoV_proj_fromVolume_median_threshold_',strrep(num2str(deep_MED_threshold),'.','-'),'.png'));
    exportgraphics(fmedian_deep,savefn_img)
    close all

    save(fullfile(app.save_comp, strcat(app.tp,'_userinput_threshold.mat')),'th','-v7.3');
    save(fullfile(app.save_comp, strcat(app.tp,'_measurement.mat')),"measurement",'-v7.3');

end

end



