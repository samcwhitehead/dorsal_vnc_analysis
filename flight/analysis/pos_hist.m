function counts = pos_hist(x_pos)
	  %divide the xpos timeseries in to 16 bins
	counts = zeros(1,16);
	for count_idx = 1:16
	    counts(count_idx) = sum(x_pos == count_idx-1);
	end
end
