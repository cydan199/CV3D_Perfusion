function [row_shift_all,col_shift_all] = matchAxial(refOCT, tarOCT, numBatch, scale)

refOCT_log  = refOCT;
tarOCT_log  = tarOCT;

numAscans   = size(refOCT_log,2);
wdtBatch    = numAscans/numBatch;

for k = 1:size(refOCT_log,3)
    
    for j = 1:numBatch
        refImg = imadjust(mat2gray(squeeze(refOCT_log(:,wdtBatch*(j-1)+1:wdtBatch*j,k))));
        tarImg = imadjust(mat2gray(squeeze(tarOCT_log(:,wdtBatch*(j-1)+1:wdtBatch*j,k))));
        
        if isnan(tarImg)
            continue;
        end

        idxBatch(j) = (wdtBatch*(j-1)+1) + floor(wdtBatch/2); % index of the centre frame in each batch
        
        [output]=dftregistration(fft2(refImg),fft2(tarImg));

        if length(output) == 2
            row_shift(j) = 0;
            col_shift(j) = 0;
        else
            row_shift(j) = output(3);
            col_shift(j) = output(4);
        end

    end
    row_shift_smooth = round(smooth(row_shift,0.5,'rlowess'));
    col_shift_smooth = round(smooth(col_shift, 0.5,'rlowess'));

    % apply shift to the centre frame in the batch only, and extrapolate
    % the other frames
    row_shift_interp = round(interp1(idxBatch,row_shift_smooth,1:numAscans,'linear','extrap'))';
    col_shift_interp = round(interp1(idxBatch,col_shift_smooth,1:numAscans,'linear','extrap'))';
    
    plot(1:numAscans,row_shift_interp,'r');hold on
    plot(1:numAscans, col_shift_interp, 'b');hold off
    legend('row shift','col shift')

    row_shift_all(:,k) = row_shift_interp;
    col_shift_all(:,k) = col_shift_interp;
end

