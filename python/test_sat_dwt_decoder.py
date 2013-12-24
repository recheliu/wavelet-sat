# ADD-BY-LEETEN 2013/09/24-BEGIN
import numpy as np;
import time;
# ADD-BY-LEETEN 2013/09/24-END
import sat_dwt_decoder;

# ADD-BY-LEETEN 2013/09/24-BEGIN
time_queue = [];

time_queue.append(time.time());
# ADD-BY-LEETEN 2013/09/24-END
waveletsat_decoder = sat_dwt_decoder.decoder();
# waveletsat_decoder.load("D:/data/WaveletSAT/video/vis08.wnc4");
waveletsat_decoder.load("D:/data/WaveletSAT/video/AVSS_PV_Hard.wnc4");
time_queue[-1] = time.time() - time_queue[-1];	# ADD-BY-LEETEN 2013/09/24

# # MOD-BY-LEETEN 2013/09/14-FROM:
# data_size = waveletsat_decoder.get_size();
# # MOD-BY-LEETEN 2013/09/14-TO:
data_size = [];
waveletsat_decoder.get_size(data_size);
# # MOD-BY-LEETEN 2013/09/14-END
print data_size;

# ADD-BY-LEETEN 2013/09/14-BEGIN
n_bins = waveletsat_decoder.get_n_bins();
print n_bins;

time_hist = np.ndarray(shape=(n_bins, data_size[2]-1), dtype=float, order='C');
# ADD-BY-LEETEN 2013/09/14-END

# # MOD-BY-LEETEN 2013/09/14-FROM:
# for t in range(0, data_size[2]-1):
# 	region_hist = waveletsat_decoder.get_region_histogram([0, 0, t], [data_size[0]-1, data_size[1]-1, t+1]);
# 	print region_hist;
# # MOD-BY-LEETEN 2013/09/14-TO:
spatial_left = [110, 102];
spatial_size = [44, 33];
spatial_right = [];
for d in range(0, len(spatial_left)):
	spatial_right.append(spatial_left[d] + spatial_size[d]);
	
time_queue.append(time.time());
for t in range(0, data_size[2]-1):
	region_hist = []; # ADD-BY-LEETEN 2013/10/20
	waveletsat_decoder.get_region_histogram([spatial_left[0], spatial_left[1], t], [spatial_right[0], spatial_right[1], t+1], region_hist);
	for b in range(0, n_bins):
		time_hist[b, t] = region_hist[b];
time_queue[-1] = time.time() - time_queue[-1];
print time_queue;

import matplotlib.pyplot as plt;		
plt.imshow(time_hist, origin='lower');
plt.show(); 
# # MOD-BY-LEETEN 2013/09/14-END

# import matplotlib.pyplot as plt;
# plt.bar(range(0, len(region_hist)), region_hist);
# plt.show();


############################################################
# Copyright (c) 2013 Teng-Yok Lee
#
# See the file LICENSE.txt for copying permission.
############################################################
