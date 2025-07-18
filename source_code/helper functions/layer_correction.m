function [row, col] = layer_correction(raw_segmentation_vol, frame)
    % define the initial volume
    row = zeros(size(raw_segmentation_vol,2),1);
    col = zeros(size(raw_segmentation_vol,3),1);

    Current_Segmentation_layer_fast = raw_segmentation_vol(:,:,frame); % Selecting the First BM Segmentation
    [row_layer, col_layer] = find(Current_Segmentation_layer_fast==1); % Finding the segmentation 1 value locations
    
    if length(col_layer) ~= length(col)
        padding_size = length(col) - length(col_layer);
        col_layer = cat(1, col_layer, zeros(padding_size,1));
    end

    count = 0;

    for jj = 1:length(row)
        
        if col_layer(jj-count) ~= jj
            
            if (jj-count) == 1
                row(jj)=row_layer(jj-count);
                col(jj) = jj;
            elseif (jj-count) > length(row_layer)
                row(jj)= row(jj-1);
                col(jj)=jj;
            else
                row(jj) = mean(row(jj-1),row_layer(jj-count));
                col(jj) = jj;
            end
            count = count+1; % location of this
        else
            row(jj) = row_layer(jj-count);
            col(jj) = col_layer(jj-count);
        end
    end

end