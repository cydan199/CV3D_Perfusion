function data_permute_flip = read_OCT_vol(app, OCT_vol_fn)
datasize = app.datasize;

% Loading Volume
data = uint8(zeros(datasize));

% fid output is PCZMI1742760577_Angio (3mmx3mm)_2-8-2022_14-26-15_OD_sn16842_cube_z_template.img
% Open file
fid = fopen(OCT_vol_fn);

for m = 1:datasize(3)
    % Read binary data from file
    data(:,:,m) = uint8(fread(fid,[datasize(1),datasize(2)],'uint8'));
end

data_permute = permute(data,[2,1,3]);
data_permute_flip = flipud(data_permute);
fclose(fid);

end