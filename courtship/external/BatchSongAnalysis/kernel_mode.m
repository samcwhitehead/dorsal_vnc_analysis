function [kmode] = kernel_mode(vector,bins)
%%Function to get the mode from a kernel distribution.
%kernel_mode(vector, 21:1:199) - vector is any vector whose mode is going
%to be obtained.
%bins - a sequence to indicate which is the size of the intervals used to
%get the mode.

x = bins;
f=ksdensity(vector,x);
max_f_col=find(f==max(f));
max_x=x(max_f_col);
kmode = max_x;
        
        