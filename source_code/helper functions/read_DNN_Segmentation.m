function data_permute = read_DNN_Segmentation(loadloc)

% Loading Volume
data = niftiread(loadloc);

data_permute = permute(data,[2,1,3]);
data_permute_flip = flipdim(data_permute,1);

end