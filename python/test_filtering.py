import numpy as np;	# ADD-BY-LEETEN 2013/09/14
import time;
import sat_dwt_decoder;

time_queue = [];
time_queue.append(time.time());

waveletsat_decoder = sat_dwt_decoder.decoder();
# waveletsat_decoder.load("D:/data/image/WaveletSAT/saltandpepper.b_256.wnc4");
waveletsat_decoder.load("D:/data/image/WaveletSAT/flower_L.b_256.wnc4");
time_queue[-1] = time.time() - time_queue[-1];

time_queue.append(time.time());
data_size = [];
waveletsat_decoder.get_size(data_size);
print data_size;

n_bins = waveletsat_decoder.get_n_bins();
print n_bins;
time_queue[-1] = time.time() - time_queue[-1];

def median_filter(hist):
	weight = sum(hist);
	median = -1;
	if( weight == 0 ):
		return median;
	c = 0.0;
	n_bins = len(hist);
	prev_c = 0.0; 
	for b in range(0, n_bins): 
		c = c + hist[b] / weight;
		if( c >= 0.5 ):
			median = b + (0.5 - prev_c)/(c - prev_c);
			break;
		prev_c = c;
	return median;

def non_zero_filter(hist):
	weight = sum(hist);
	if( weight != 0.0 ):
		for b in range(0, n_bins): 
			prob = hist[b] / weight; 
			if( prob >= 0.9 ):
				return b;
	return 0;
	
n_dims = len(data_size);
print time_queue;

import matplotlib.pyplot as plt;
for w in (1, 3, 5, 10, 20, 50):
	time_queue = [];
	
	left_offset = [];
	right_offset = [];
	for d in range(0, n_dims):
		left_offset.append(	-w);
		right_offset.append(+w);
	
	
	time_queue.append(time.time());
	result = [];
	waveletsat_decoder.apply_filter(median_filter, left_offset, right_offset, result, True);
	# waveletsat_decoder.apply_filter(non_zero_filter, [-1 -1], [0, 0], result);
	time_queue[-1] = time.time() - time_queue[-1];
	
	print time_queue;

 	if( n_dims == 2 ):
 		filtered_result = np.array(result);
		filtered_array = filtered_result.reshape(data_size[1], data_size[0], order='C').copy();
		
	 	plt.figure();
	 	plt.imshow(filtered_array);
	 	plt.set_cmap("gray");
	 	plt.colorbar();
	 	plt.title(("[%d %d] - [%d %d]" % (left_offset[0], left_offset[1], right_offset[0], right_offset[1])));

plt.show(); 
    

# if( n_dims == 2 ):
# 	filtered_result = np.array(result);
# 	filtered_array = filtered_result.reshape(data_size[1], data_size[0], order='C').copy();
        
# 	plt.figure();
# 	plt.imshow(filtered_array);
# 	plt.set_cmap("gray");
# 	plt.colorbar();
# 	plt.show(); 
#    
