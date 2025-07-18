function imageProjection(app)

saveloc = fullfile(app.loadloc,'cross_comparison');
if ~exist(saveloc,'dir')
    mkdir(saveloc)
end

cleanOutPut(app,'***Start of Noise Removal & Image Projection***');
ftemp = 156;

for m = 1:length(app.timepoints)

    app.tp = app.timepoints{m};

    app.savefns = fullfile(app.save_folder, app.tp);

    %%% define saving folders
    app.loadloc_enface = fullfile(app.loadloc, app.disease, app.subject, app.tp);
    app.loadloc_cube = fullfile(app.loadloc, app.disease, app.subject, strcat(app.subject,"_CUBE"), app.tp);

    app.save_registered = fullfile(app.loadloc_cube,'CropRegistered');
    app.save_segmentation =fullfile(app.loadloc_cube,'CropSegmentation');
    app.save_volume = fullfile(app.savefns,'Cropped_V');
    app.save_mat = fullfile(app.savefns,'mat');
    app.save_cov = fullfile(app.savefns,'CoV');
    app.save_mat_proj = fullfile(app.save_mat,app.proj_layer);
    app.save_cov_proj = fullfile(app.save_cov,app.proj_layer);

    if exist(fullfile(app.save_mat,'offset.mat'),'file')
        offset = changeOffsetValue(importdata(fullfile(app.save_mat,'offset.mat')));
    else

        urinput = inputdlg(strcat("what is the offset value for noise estimation at: ",app.tp),'Enter Offset Value');
        if ~isnan(str2double(urinput))
            offset = str2double(urinput);
        else
            warndlg('Invalid input. Please enter a valid number.','Input Error');
        end
    end
    save(fullfile(app.save_mat,'offset.mat'),'offset','-v7.3');
    app.showOffset.Value = offset;
    outPut(app, [strcat(app.tp, " offset value is:"), num2str(offset)]);

    if exist(fullfile(app.save_mat,'cropped_amount.mat'),'file')
        app.cropped_amount = importdata(fullfile(app.save_mat,'cropped_amount.mat'));
    end

    app.depthROI_CoV = importdata(fullfile(app.save_mat,'depthROI_CoV.mat'));

    % read segmentation file
    [segmentation_fns, segmentation_loc] = uigetfile(fullfile(app.save_segmentation,...
        strcat(app.subject,'_',app.tp,'_Angiography_3x3_cube_z_all'),...
        'Boundary.nii.gz'),...
        'Choose layer segmentation (generate by OCTSurfer)');
    outPut(app, ['Segmentation directory: ', segmentation_fns]);
    app.LayersegmentationFn.Value = fullfile(segmentation_loc,segmentation_fns);

    fig = app.UIFigure;
    d = uiprogressdlg(fig,'Title','Generating projection images...',...
        'Indeterminate','on','Cancelable','on');
    drawnow
    d.Message = strcat("Processing data at timepoint ", app.tp);

    segmentation = double(read_DNN_Segmentation(fullfile(segmentation_loc,segmentation_fns)));
    segmentation = segmentation(app.depthROI_CoV(1):app.depthROI_CoV(2),app.cropped_amount+1:end-app.cropped_amount,app.cropped_amount+1:end-app.cropped_amount);

    segILM = segmentation;
    segILM(~(segILM == 1))= 0;

    segNFL = segmentation;
    segNFL(~(segNFL == 2))= 0;
    segNFL((segNFL == 2))= 1;

    segOPL = segmentation;
    segOPL(~(segOPL == 3))= 0;
    segOPL(segOPL == 3)= 1;

    segONL = segmentation;
    segONL(~(segONL == 4))= 0;
    segONL(segONL == 4)= 1;

    if d.CancelRequested == 1
        % app.CancelFlag = true;
        outPut(app, 'Operation cancelled by user.');
        return
    end

    vol_check = [];

    between_segmentation = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);
    between_segmentation_sup = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);
    between_segmentation_deep = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);

    idx_between_segmentation = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);
    idx_between_segmentation_sup = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);
    idx_between_segmentation_deep = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);

    OCTA_fns = dir(fullfile(app.save_registered, '*FlowCube*.img'));

    retina_reg_V = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount,length(OCTA_fns));
    sup_reg_V = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount,length(OCTA_fns));
    deep_reg_V = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount,length(OCTA_fns));

    maxVal_idx = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount,length(OCTA_fns));
    maxVal_sup_idx = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount,length(OCTA_fns));
    maxVal_deep_idx = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount,length(OCTA_fns));


    OCTA_averaged = zeros(app.depthROI_CoV(2)-app.depthROI_CoV(1)+1,app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);
    for i = 1:length(OCTA_fns)
        if d.CancelRequested == 1
            outPut(app, 'Operation cancelled by user.');
            return
        end
        close all
        if app.CancelFlag
            outPut(app, 'Operation cancelled by user.');
            return
        end

        [~,name,~] = fileparts(fullfile(OCTA_fns(i).folder, OCTA_fns(i).name));

        OCTA_vol = double(read_OCT_vol(app, fullfile(OCTA_fns(i).folder, OCTA_fns(i).name)));
        OCTA_vol = OCTA_vol(app.depthROI_CoV(1):app.depthROI_CoV(2),app.cropped_amount+1:end-app.cropped_amount,app.cropped_amount+1:end-app.cropped_amount);

        for ii=1:size(segNFL,3)
            if d.CancelRequested == 1
                % app.CancelFlag = true;
                outPut(app, 'Operation cancelled by user.');
                return
            end

            if app.CancelFlag
                outPut(app, 'Operation cancelled by user.');
                return
            end

            [row_NFL, ~] = layer_correction(segNFL, ii); 
            [row_OPL, ~] = layer_correction(segOPL, ii); 
            [row_ONL, ~] = layer_correction(segONL, ii); 

            row_NFL_smooth= row_NFL;
            row_OPL_smooth = row_OPL;
            row_ONL_smooth = row_ONL;

            for jj=1:size(row_NFL_smooth,1)
                if app.CancelFlag
                    outPut(app, 'Operation cancelled by user.');
                    return
                end
                %%% only for the purpose of wrong segmentation
                if (row_NFL_smooth(jj)>row_ONL_smooth(jj))
                    Values_between_Segmentation = OCTA_vol(row_ONL_smooth(jj):row_NFL_smooth(jj), jj, ii);
                    disp('segmentation error!')
                else
                    Values_between_Segmentation =  OCTA_vol(row_NFL_smooth(jj):row_ONL_smooth(jj), jj, ii);
                end

                if (row_NFL_smooth(jj)>row_OPL_smooth(jj))
                    val_between_sup = OCTA_vol(row_OPL_smooth(jj):row_NFL_smooth(jj), jj, ii);
                    disp('segmentation error at OPL!')
                else
                    val_between_sup = OCTA_vol(row_NFL_smooth(jj):row_OPL_smooth(jj), jj, ii);
                end

                if (row_OPL_smooth(jj)>row_ONL_smooth(jj))
                    val_between_deep = OCTA_vol(row_ONL_smooth(jj):row_OPL_smooth(jj), jj, ii);
                    disp('segmentation error at ONL!')
                else
                    val_between_deep = OCTA_vol(row_OPL_smooth(jj):row_ONL_smooth(jj), jj, ii);
                end


                [Max_OCTA_Column_fast, max_OCTA_idx] = max(Values_between_Segmentation);
                [Max_OCTA_sup_fast, max_OCTA_sup_idx] = max(val_between_sup);
                [Max_OCTA_deep_fast, max_OCTA_deep_idx] = max(val_between_deep);


                between_segmentation(ii,jj) = Max_OCTA_Column_fast;
                between_segmentation_sup(ii,jj)=Max_OCTA_sup_fast;
                between_segmentation_deep(ii,jj)=Max_OCTA_deep_fast;


                idx_between_segmentation(ii,jj) = max_OCTA_idx;
                idx_between_segmentation_sup(ii,jj) = max_OCTA_sup_idx;
                idx_between_segmentation_deep(ii,jj) = max_OCTA_deep_idx;
            end

            % save  of three volumes for check
            if ii == ftemp
                if d.CancelRequested == 1
                    app.CancelFlag = true;
                    return
                end
                vol_check(:,:,i) = OCTA_vol(:,:,ii);
            end
        end

        saveAs3DTiffStack(vol_check,fullfile(app.save_folder,app.tp,'Cropped_V',strcat('frame_',num2str(ftemp),'_fromVolume.tif')));
        save(fullfile(app.save_mat,strcat('frame_',num2str(ftemp),'_fromVolume.mat')),"vol_check","-v7.3");

        retina_reg_V(:,:,i) = imrotate(fliplr(between_segmentation),90);
        sup_reg_V(:,:,i) = imrotate(fliplr(between_segmentation_sup),90);
        deep_reg_V(:,:,i) = imrotate(fliplr(between_segmentation_deep),90);

        maxVal_idx(:,:,i) = imrotate(fliplr(idx_between_segmentation),90);% ALTERNATIVE: fliplr(above_segmentation_temp);
        maxVal_sup_idx(:,:,i) = imrotate(fliplr(idx_between_segmentation_sup),90);% ALTERNATIVE: fliplr(above_segmentation_temp);
        maxVal_deep_idx(:,:,i) = imrotate(fliplr(idx_between_segmentation_deep),90);% ALTERNATIVE: fliplr(above_segmentation_temp);

        OCTA_averaged = OCTA_averaged + OCTA_vol; clear OCTA_vol
    end

    % compute CoV
    %%% Indices %%%
    savefn_mat_maxValIdx = fullfile(app.save_mat_proj,strcat(app.tp,'_maxValIdx_fromVolume_CoV','.mat'));
    compute_CoV_2D(maxVal_idx,savefn_mat_maxValIdx);

    savefn_mat_maxValIdx_sup = fullfile(app.save_mat_proj,strcat(app.tp,'_maxValIdx_sup_fromVolume_CoV','.mat'));
    compute_CoV_2D(maxVal_sup_idx,savefn_mat_maxValIdx_sup);

    savefn_mat_maxValIdx_deep = fullfile(app.save_mat_proj,strcat(app.tp,'_maxValIdx_deep_fromVolume_CoV','.mat'));
    compute_CoV_2D(maxVal_deep_idx,savefn_mat_maxValIdx_deep);

    %%% Intensity %%%
    savefn_mat_V = fullfile(app.save_mat_proj,strcat(app.tp,'_OCTA_fromVolume_CoV','.mat'));
    compute_CoV_2D(retina_reg_V,savefn_mat_V);

    savefn_mat_V_sup = fullfile(app.save_mat_proj,strcat(app.tp,'_OCTA_sup_fromVolume_CoV','.mat'));
    compute_CoV_2D(sup_reg_V,savefn_mat_V_sup);

    savefn_mat_V_deep = fullfile(app.save_mat_proj,strcat(app.tp,'_OCTA_deep_fromVolume_CoV','.mat'));
    compute_CoV_2D(deep_reg_V,savefn_mat_V_deep);

    retina_reg_Z = importdata(fullfile(app.save_mat,strcat('registered_',app.tp,'_enfaceOCTA_fromZeiss.mat')));
    retina_reg_Z = retina_reg_Z(app.cropped_amount+1:end-app.cropped_amount,app.cropped_amount+1:end-app.cropped_amount,:);
    savefn_mat_Z = fullfile(app.save_mat_proj,strcat(app.tp,'_OCTA_fromZeiss_CoV','.mat'));
    compute_CoV_2D(retina_reg_Z,savefn_mat_Z);

    save(fullfile(app.save_mat_proj,strcat('registered_',app.tp,'_maxValIDX_fromVolume.mat')),"maxVal_idx",'-v7.3');
    save(fullfile(app.save_mat_proj,strcat('registered_',app.tp,'_maxValIDX_sup_fromVolume.mat')),"maxVal_sup_idx",'-v7.3');
    save(fullfile(app.save_mat_proj,strcat('registered_',app.tp,'_maxValIDX_deep_fromVolume.mat')),"maxVal_deep_idx",'-v7.3');

    save(fullfile(app.save_mat_proj,strcat('registered_',app.tp,'_enfaceOCTA_fromVolume.mat')),"retina_reg_V",'-v7.3');
    saveAs3DTiffStack(retina_reg_V,fullfile(app.save_cov_proj,strcat('registered_',app.tp,'_enfaceOCTA_fromVolume.tif')));

    save(fullfile(app.save_mat_proj,strcat('registered_',app.tp,'_enfaceOCTA_sup_fromVolume.mat')),"sup_reg_V",'-v7.3');
    save(fullfile(app.save_mat_proj,strcat('registered_',app.tp,'_enfaceOCTA_deep_fromVolume.mat')),"deep_reg_V",'-v7.3');

    outPut(app, ['Custome en face OCTA images are saved as: ', fullfile(app.save_cov_proj,strcat('registered_',app.tp,'_enfaceOCTA_fromVolume.tif'))]);

    OCTA_averaged = OCTA_averaged./length(OCTA_fns);
    save(fullfile(app.save_mat,'averaged_regOCTA_vol.mat'),'OCTA_averaged','-v7.3');

    clearvars savefn_mat_V savefn_mat_Z retina_reg_V retina_reg_Z

    %%% Compute the CoV of OCTA Bscans
    % define/read the depthROI for CoV computation
    if exist(fullfile(app.save_mat,'depthROI_CoV.mat'),'file')
        app.depthROI_CoV = importdata(fullfile(app.save_mat,'depthROI_CoV.mat'));
    else
        errordlg('The ROI for CoV computation is missing.','File Not Found');
    end

    %%% Make binarized map from averaged OCTA B-scans
    threshold_BGD = zeros(size(OCTA_averaged,3),1);
    binarized_OCTA = zeros(size(OCTA_averaged)); % (depthROI(1):depthROI(2),app.cropped_amount+1:end-app.cropped_amount,app.cropped_amount+1:end-app.cropped_amount);
    avg_OCTA_masked = OCTA_averaged;

    % load Volumetric CoV
    if ~exist(fullfile(app.save_mat,'CoV_fromBscans.mat'),'file')
        errordlg('No Volumetric CoV is found.','File Not Found');
    else
        Bscans_CoV = importdata(fullfile(app.save_mat,'CoV_fromBscans.mat'));
        bscans_cov_crop = Bscans_CoV(:,app.cropped_amount+1:end-app.cropped_amount,app.cropped_amount+1:end-app.cropped_amount);
        bscans_cov_crop(isnan(bscans_cov_crop))=0;
        if d.CancelRequested == 1
            outPut(app, 'Operation cancelled by user.');
            return
        end
        CoV_Bscans_masked = bscans_cov_crop;
    end

    VRI_shifting = 10;
    VRI_thickness = 20;

    if exist(fullfile(app.save_volume,strcat('Binarized_avgOCTA_mask_offset_',strrep(num2str(offset),'.','-'),'.tif')),'file')
        delete(fullfile(app.save_volume,strcat('Binarized_avgOCTA_mask_offset_',strrep(num2str(offset),'.','-'),'.tif')));
    end
    if exist(fullfile(app.save_volume,strcat('Binarized_avgOCTA_vol_offset_',strrep(num2str(offset),'.','-'),'.tif')),'file')
        delete(fullfile(app.save_volume,strcat('Binarized_avgOCTA_vol_offset_',strrep(num2str(offset),'.','-'),'.tif')));
    end

    if exist(fullfile(app.save_volume,strcat('Segmentation_',app.proj_layer,'.tif')),'file')
        delete(fullfile(app.save_volume,strcat('Segmentation_',app.proj_layer,'.tif')));
    end

    if exist(fullfile(app.save_volume,strcat('Binarized_CoV_fromOCTA_Bscans_offset_',strrep(num2str(offset),'.','-'),'.tif')),'file')
        delete(fullfile(app.save_volume,strcat('Binarized_CoV_fromOCTA_Bscans_offset_',strrep(num2str(offset),'.','-'),'.tif')));
    end

    if exist(fullfile(app.save_volume,'CoV_fromOCTA_Bscans.tif'),'file')
        delete(fullfile(app.save_volume,'CoV_fromOCTA_Bscans.tif'));
    end

    %%% define the variable for CoV projection image
    cm = getcolormap_sup();

    between_segmentation_median = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);
    between_segmentation_mean = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);
    covBscan_median = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);
    covBscan_mean = zeros(app.numAscans-2*app.cropped_amount,app.numBscans-2*app.cropped_amount);

    for num_frame = 1:size(OCTA_averaged,3)
        if d.CancelRequested == 1
            outPut(app, 'Operation cancelled by user.');
            return
        end

        % binarize averaged OCTA per B-scan location
        BGD_intensity = [];
        [row_ILM, col_ILM] = layer_correction(segILM, num_frame); 
        row_ILM_smooth= row_ILM;

        for jj=1:size(row_ILM_smooth,1)
            if d.CancelRequested == 1
                outPut(app, 'Operation cancelled by user.');
                return
            end

            upper_boundary_row = row_ILM_smooth(jj)-VRI_shifting-VRI_thickness+1;
            lower_bloundary_row = row_ILM_smooth(jj)-VRI_shifting;

            if upper_boundary_row <= 0 || lower_bloundary_row <= 0
               upper_boundary_row = 1;
               lower_bloundary_row = 20;
            end

            Values_between_Segmentation =  OCTA_averaged(upper_boundary_row:lower_bloundary_row, jj, num_frame);
            BGD_intensity = [BGD_intensity,Values_between_Segmentation(:)];
        end

        threshold_BGD(num_frame) = mean(BGD_intensity(:),'all') + std(BGD_intensity(:),0,'all') + offset;

        avg_OCTA_mask_temp = OCTA_averaged(:,:,num_frame);
        avg_OCTA_mask_temp(avg_OCTA_mask_temp <= threshold_BGD(num_frame)) = 0;
        avg_OCTA_mask_temp(avg_OCTA_mask_temp > threshold_BGD(num_frame)) = 1;

        binarized_OCTA(:,:,num_frame) = avg_OCTA_mask_temp;

        % multiply mask with OCTA
        avg_OCTA_masked(:,:,num_frame) = OCTA_averaged(:,:,num_frame).*binarized_OCTA(:,:,num_frame);
        CoV_Bscans_masked(:,:,num_frame) = bscans_cov_crop(:,:,num_frame).*binarized_OCTA(:,:,num_frame);

        [row_NFL, col_NFL] = layer_correction(segNFL, num_frame); 
        [row_OPL, col_OPL] = layer_correction(segOPL, num_frame); 
        [row_ONL, col_ONL] = layer_correction(segONL, num_frame); 

        row_NFL_smooth= row_NFL;
        row_OPL_smooth= row_OPL;
        row_ONL_smooth = row_ONL;

        for jj=1:size(row_NFL_smooth,1)
            if d.CancelRequested == 1
                outPut(app, 'Operation cancelled by user.');
                return
            end
            %%% only for the purpose of wrong segmentation
            if (row_NFL_smooth(jj) > row_ONL_smooth(jj))
                between_segmentation(num_frame,jj) = 0;
                between_segmentation_median(num_frame,jj) = 0;
                between_segmentation_mean(num_frame,jj) = 0;

                between_segmentation_sup(num_frame,jj) = 0;
                between_segmentation_deep(num_frame,jj) = 0;

                disp('segmentation error!')
            else
                Values_between_Segmentation_OCTA =  avg_OCTA_masked(row_NFL_smooth(jj):row_ONL_smooth(jj), jj, num_frame);

                val_between_CoV = nonzeros(CoV_Bscans_masked(row_NFL_smooth(jj):row_ONL_smooth(jj), jj, num_frame));
                val_between_sup_CoV = nonzeros(CoV_Bscans_masked(row_NFL_smooth(jj):row_OPL_smooth(jj), jj, num_frame));
                val_between_deep_CoV = nonzeros(CoV_Bscans_masked(row_OPL_smooth(jj):row_ONL_smooth(jj), jj, num_frame));



                Max_OCTA_Column_fast = max(Values_between_Segmentation_OCTA);
                between_segmentation(num_frame,jj) = Max_OCTA_Column_fast;
                %%% median projection %%%
                median_CoV_colume_fast = median(val_between_CoV);
                between_segmentation_median(num_frame,jj) = median_CoV_colume_fast;

                median_CoV_colume_sup_fast = median(val_between_sup_CoV);
                between_segmentation_sup_median(num_frame,jj) = median_CoV_colume_sup_fast;

                median_CoV_colume_deep_fast = median(val_between_deep_CoV);
                between_segmentation_deep_median(num_frame,jj) = median_CoV_colume_deep_fast;
                %%% mean projection %%%
                Mean_CoV_Column_fast = mean(val_between_CoV);
                between_segmentation_mean(num_frame,jj) = Mean_CoV_Column_fast;
            end
        end

        imwrite(mat2gray(CoV_Bscans_masked(:,:,num_frame)).*255,cm,fullfile(app.save_volume,strcat('Binarized_CoV_fromOCTA_Bscans_offset_',strrep(num2str(offset),'.','-'),'.tif')),'WriteMode','append');
        imwrite(mat2gray(bscans_cov_crop(:,:,num_frame)).*255,cm,fullfile(app.save_volume,'CoV_fromOCTA_Bscans.tif'),'WriteMode','append');

        %%% save segmentation results to check the performance
        with_seg_zeros = mat2gray(avg_OCTA_masked(:,:,num_frame));
        ind = sub2ind(size(with_seg_zeros),row_NFL_smooth,col_NFL);
        with_seg_zeros(ind)=1;

        with_seg_ones = mat2gray(avg_OCTA_masked(:,:,num_frame));
        ind = sub2ind(size(with_seg_ones),row_ONL_smooth,col_ONL);
        with_seg_ones(ind)=1;
        fused=imfuse(with_seg_zeros,with_seg_ones);
        imwrite(fused,fullfile(app.save_volume,strcat('Segmentation_',app.proj_layer,'.tif')),'WriteMode','append');

        imwrite(mat2gray(avg_OCTA_masked(:,:,num_frame)),fullfile(app.save_volume,strcat('Binarized_avgOCTA_vol_offset_',strrep(num2str(offset),'.','-'),'.tif')),'WriteMode','append');
        imwrite(mat2gray(avg_OCTA_mask_temp),fullfile(app.save_volume,strcat('Binarized_avgOCTA_mask_offset_',strrep(num2str(offset),'.','-'),'.tif')),'WriteMode','append');
    end

    d.Message = strcat("Saving results at timepint ", app.tp);

    CoV_Bscans_masked_ref = CoV_Bscans_masked(:,:,ftemp);

    % save binarized CoV as mat
    % save(fullfile(app.save_mat,strcat('Binarized_CoV_fromOCTA_Bscans_offset_',strrep(num2str(offset),'.','-'),'.mat')),"CoV_Bscans_masked","-v7.3");
    save(fullfile(app.save_mat,strcat('Ref_Binarized_CoV_fromOCTA_Bscans_offset_',strrep(num2str(offset),'.','-'),'.mat')),"CoV_Bscans_masked_ref","-v7.3");


    save(fullfile(app.save_mat,strcat('Binarization_avgOCTA_threshold_offset_',strrep(num2str(offset),'.','-'),'.mat')),'threshold_BGD','-v7.3')
    save_maskOCTA_fn = fullfile(app.save_mat,strcat('Binarized_avgOCTA_mask_offset_',strrep(num2str(offset),'.','-'),'.mat'));
    save(save_maskOCTA_fn,'binarized_OCTA','-v7.3');

    h = figure('Visible','off');
    plot(threshold_BGD);
    xlabel('Frame Number')
    ylabel('Signal Intensity')
    save_img_fn = fullfile(app.save_mat, strcat('Binarization_avgOCTA_threshold_offset_',strrep(num2str(offset),'.','-'),'.png'));
    exportgraphics(h, save_img_fn)

    avg_OCTA_proj = imrotate(fliplr(between_segmentation),90);
    imwrite(mat2gray(avg_OCTA_proj),fullfile(app.save_cov_proj,strcat('Binarized_avgOCTA_proj_offset_',strrep(num2str(offset),'.','-'),'.tif')));

    save_projOCTA_fn = fullfile(app.save_cov_proj,strcat('Binarized_avgOCTA_proj_offset_',strrep(num2str(offset),'.','-'),'.mat'));
    save(save_projOCTA_fn,'avg_OCTA_proj','-v7.3');

    save_binarizedOCTA_fn = fullfile(app.save_mat,strcat('Binarized_avgOCTA_offset_',strrep(num2str(offset),'.','-'),'.mat'));
    save(save_binarizedOCTA_fn,'avg_OCTA_masked','-v7.3');

    covBscan_median = imrotate(fliplr(between_segmentation_median),90);
    covBscan_median_sup = imrotate(fliplr(between_segmentation_sup_median),90);
    covBscan_median_deep = imrotate(fliplr(between_segmentation_deep_median),90);

    covBscan_mean = imrotate(fliplr(between_segmentation_mean),90);

    save_covMedian_fn = fullfile(app.save_mat_proj,strcat('CoV_proj_fromBscans_median_offset_',num2str(offset),'.mat'));
    save_covMedian_sup_fn = fullfile(app.save_mat_proj,strcat('CoV_proj_sup_fromBscans_median_offset_',num2str(offset),'.mat'));
    save_covMedian_deep_fn = fullfile(app.save_mat_proj,strcat('CoV_proj_deep_fromBscans_median_offset_',num2str(offset),'.mat'));

    save_covMean_fn = fullfile(app.save_mat_proj,strcat('CoV_proj_fromBscans_mean_offset_',num2str(offset),'.mat'));

    save(save_covMedian_fn,'covBscan_median','-v7.3');
    save(save_covMedian_sup_fn,'covBscan_median_sup','-v7.3');
    save(save_covMedian_deep_fn,'covBscan_median_deep','-v7.3');

    save(save_covMean_fn,'covBscan_mean','-v7.3');
end
d.close;
outPut(app, "Image projection completed");
end
