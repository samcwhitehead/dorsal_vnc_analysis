function index = argmax(varargin)
%ARGMAX  Index of maximum element.
% ARGMAX(X) returns an index I such that X(I) == MAX(X(:)).
%

[~,index] = max(varargin{:});