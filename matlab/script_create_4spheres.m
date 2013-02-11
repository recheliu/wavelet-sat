clear all;
clear java;

n_dims = 3;
V_len = 64;
V_size = V_len * ones(1, n_dims);
V = zeros(V_size);
[X, Y, Z] = meshgrid(1:V_size(2), 1:V_size(1), 1:V_size(3));

radius = [8 4 2 1 1];
centers = [48 24 12 6 2];
n_spheres = length(centers);
for si = 1:n_spheres
    sigma = 1.1 * radius(si);
    G = exp(-((X - centers(si)).^2 + (Y - centers(si)).^2 + (Z - centers(si)).^2) / (2 * sigma^2)); % /(sigma * sqrt(2 * pi));
    G = G / max(G(:));
    V = V + G;
%     M = zeros(V_size);
%     M(G>0.05) = 1;
%     V = V + M;
end
    
filename = '4spheres';
fid = fopen([filename '.raw'], 'wb');
fwrite(fid, V, 'single');
fclose(fid);

fid = fopen([filename '.nhdr'], 'wt');
fprintf(fid, 'NRRD0001\n');
fprintf(fid, 'content: %s.raw\n', filename);
fprintf(fid, 'dimension: %d\n', n_dims);
fprintf(fid, 'type: float\n');
fprintf(fid, 'sizes: %d %d %d\n', V_size(2), V_size(1), V_size(3));
fprintf(fid, 'spacings: 1 1 1\n');
fprintf(fid, 'endian: little\n');
fprintf(fid, 'data file: ./%s.raw\n', filename);
fprintf(fid, 'encoding: raw\n');
fclose(fid);










