function [sortA,i2,j] = myunique(A)
%% MYUNIQUE the same input and output as unique 
%
% Solve the change of unique command in different version of MATLAB
%
% See also: unique
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

matlabversion = version;
if 1
    [sortA, i2, j] = unique(A,'rows','legacy'); %#ok<*ASGLU>
else
    [sortA, i2, j] = unique(A,'rows');
end