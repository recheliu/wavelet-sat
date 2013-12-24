clear all;
clear java;

% nc_filepath = 'D:\projects\WaveletSAT\build\vs2010\x64\output\VEL.t.20c.20010201.b_64.wnc4';
% nc_filepath = 'D:\projects\WaveletSAT\build\vs2010\x64\output\nucleon.b_32.z_1.wnc4';
% nc_filepath = 'D:\data\volumes\rawfiles\simulated\nucleon.nhdr.wnc4';
nc_filepath = 'D:\data\volumes\rawfiles\test\4spheres.nhdr.wnc4';
min_bin = 3;

%% load the coefficients
coef_counts     = ncread(nc_filepath, 'COEF_COUNT');
coef_offsets    = ncread(nc_filepath, 'COEF_OFFSET');
coef_bins       = ncread(nc_filepath, 'COEF_BIN');
coef_values     = ncread(nc_filepath, 'COEF_VALUE');

%% load the dimension information
coef_size = size(coef_counts);
n_dims = length(coef_size);
ncid = netcdf.open(nc_filepath,'NC_NOWRITE');
data_size = [];
for dim = 1:n_dims
    data_dim_name = sprintf('DATA_DIM_%d', dim - 1);
    data_dim_id = netcdf.inqDimID(ncid, data_dim_name);
    [data_dim_name, data_dim_len] = netcdf.inqDim(ncid, data_dim_id);
    data_size(end+1) = data_dim_len;
end
[bin_dim_name, n_bins] = netcdf.inqDim(ncid, netcdf.inqDimID(ncid, 'BIN'));
netcdf.close(ncid);

%% now plot the histogram
n_levels = 6;
figure;
for li = 1:n_levels
    if( li == 1 ) 
        coef_base = 0;
        coef_local_len = 1;
    else
        coef_base = 2^(li - 2);
        coef_local_len = 2^(li - 2);
    end
        
    n_local_hists = coef_local_len^n_dims;
    coef_hists = zeros(n_local_hists, n_bins);
    for lhi = 1:n_local_hists
        lh_subs = zeros(1, n_dims);
        switch n_dims
            case 2
                [lh_subs(1), lh_subs(2)] = ind2sub(coef_local_len * ones(1, n_dims), lhi);
                lh_global_subs = lh_subs + coef_base;
                lh_global_id = sub2ind(coef_size, lh_global_subs(1), lh_global_subs(2));
            case 3
                [lh_subs(1), lh_subs(2), lh_subs(3)] = ind2sub(coef_local_len * ones(1, n_dims), lhi);
                lh_global_subs = lh_subs + coef_base;
                lh_global_id = sub2ind(coef_size, lh_global_subs(1), lh_global_subs(2), lh_global_subs(3));
        end
        coef_ids = cast(coef_offsets(lh_global_id), 'double') + (0:coef_counts(lh_global_id)-1) + 1;
        coef_hists(lhi, coef_bins(coef_ids) + 1) = coef_values(coef_ids);
    end
    coef_hists = coef_hists ./ (sum(coef_hists, 2) * ones(1, n_bins));
    subplot(n_levels, 1, li);
    set(gca, 'FontSize', 25.0);
    plot(1:n_bins, coef_hists);
    xlim([min_bin, n_bins]);
    if( 1 == li  )
        title(nc_filepath);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2013 Teng-Yok Lee
%
% See the file LICENSE.txt for copying permission.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
