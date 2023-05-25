function width = calc_hwm(pos_bins)
	% given a position time series, return the minimum histogram width
	% from center to contain 50% of the values
	% returns half of the width
	total_pos = sum(pos_bins);
	center_sum = 0;
	width = 0;
	while center_sum < 0.5 * total_pos
	    width = width + 1;
	    center_sum = sum(pos_bins(8-width:7+width));
	end
end
