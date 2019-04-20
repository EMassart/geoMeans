function [ M, tTot, varargout] = meanArithmetic( A)
%MEANARITHMETIC(A) computes the arithmetic mean of a set of matrices.
%A is a 3D array containing the matrices along the third dimension
%M is the arithmetic mean and tTot is the CPU time used by the function

tStart = cputime;
[~,~,n] = size(A);
M = sum(A,3)./n;
tTot = cputime-tStart;

end