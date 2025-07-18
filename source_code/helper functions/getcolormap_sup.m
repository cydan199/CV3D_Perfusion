function [cm] = getcolormap_sup()

mat = zeros(256,3); % first row [0 0 0] for black (no variance)
mat(2:64,3)= linspace(0.5, 1, 63)'; %dark blue to blue
mat(65:96,2)= linspace(0, 1, 32)'; %blue to cyan
mat(65:96,3)= ones(32,1)'; %blue to cyan

mat(97:160,2) = ones(64,1)'; %cyan to green
mat(97:160,3) = linspace(1,0,64)';

mat(161:224,2) = ones(64,1)'; %green to yellow
mat(161:224,1) = linspace(0,1,64)';

mat(225:256,1) = ones(32,1)'; %yellow to red
mat(225:256,2) = linspace(1,0,32)';

cm=mat;
end

