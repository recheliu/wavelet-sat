import numpy as np;    
import math;            
import sat_dwt_decoder;

waveletsat_decoder = sat_dwt_decoder.decoder();
# waveletsat_decoder.load("D:/data/WaveletSAT/video/vis08.wnc4");
waveletsat_decoder.load("D:/data/WaveletSAT/video/AVSS_AB_Easy.wnc4");

data_size = [];
waveletsat_decoder.get_size(data_size);
print data_size;

data_level = [];
waveletsat_decoder.get_level(data_level);
print data_level;

n_bins = waveletsat_decoder.get_n_bins();
print n_bins;

query_time = 50;
query_coord = [27, 81];
query_size = [16, 32];

query_level = [];
query_local_coord = [];
for d in range(0, len(query_size)):
    qsize = query_size[d];
    query_local_coord.append(int(math.floor(query_coord[d]/qsize)));
    query_level.append(data_level[d] - int(math.floor(math.log(qsize, 2))));
    
# append the 3rd dim for time
query_local_coord.append(0);  
query_level.append(0);  

# import matplotlib.pyplot as plt;        
# sparse_coef = {};
# waveletsat_decoder.get_coef_sparse([4, 4, 4], [3, 3, 3], sparse_coef);
#     
# value_sum = sum(sparse_coef.values());
# coef_hist = [];
# for b in range(0, n_bins):
#     value = 0;
#     if( sparse_coef.has_key(b) ):
#         value = sparse_coef[b]/value_sum;
#     coef_hist.append(value);
# 
# plt.subplot(data_level[2], 1, tlevel);
# plt.bar(range(0, n_bins), coef_hist);
# plt.xlim(0, n_bins)
# plt.show(); 
    
import matplotlib.pyplot as plt;        
time_window = 1 << (data_level[2]);
prev_cdf = [];
level_dist = [];

plt.figure();
for tlevel in range(0, data_level[2]):
 
    if( tlevel > 0 ):
        time_window >>= 1;
     
    query_level[2] = tlevel;
    query_local_coord[2] = int(math.floor(query_time / time_window));
 
    print time_window;
    print query_level;
    print query_local_coord;
     
    sparse_coef = {};
    waveletsat_decoder.get_coef_sparse(query_level, query_local_coord, sparse_coef);
     
    value_sum = sum(sparse_coef.values());
    coef_hist = np.zeros( (n_bins, 1) );
    for b in range(0, n_bins):
        value = 0;
        if( sparse_coef.has_key(b) ):
            value = sparse_coef[b]/value_sum;
        coef_hist[b] = value;
        
    coef_cdf = np.cumsum(coef_hist);
    if( tlevel > 0 ):
        coef_diff = np.subtract(coef_cdf, prev_cdf);
        level_dist.append(math.sqrt(np.sum(np.multiply(coef_diff, coef_diff)))); 
        
    prev_cdf = coef_cdf;
     
    plt.subplot(data_level[2], 1, tlevel + 1);
    plt.bar(range(0, n_bins), coef_hist);
    plt.xlim(0, n_bins);
    plt.title("Level %d" % (tlevel))
 
plt.show(); 
     
plt.figure();
plt.plot(level_dist);
plt.show();
