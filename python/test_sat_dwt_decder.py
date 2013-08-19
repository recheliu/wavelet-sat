import sat_dwt_decoder;

waveletsat_decoder = sat_dwt_decoder.decoder();
waveletsat_decoder.load("D:/data/WaveletSAT/video/vis08.wnc4");

data_size = waveletsat_decoder.get_size();
print data_size;

for t in range(0, data_size[2]-1):
	region_hist = waveletsat_decoder.get_region_histogram([0, 0, t], [data_size[0]-1, data_size[1]-1, t+1]);
	print region_hist;

# import matplotlib.pyplot as plt;
# plt.bar(range(0, len(region_hist)), region_hist);
# plt.show();
