function [ M, tTot, varargout] = meanHarmonic( A)
%MEANHARMONIC(A) computes the harmonic mean of a set of matrices.
%A is a 3D array containing the matrices along the third dimension
%M is the harmonic mean and tTot is the CPU time used by the function

tStart = cputime;
[m,~,n] = size(A);
MH = zeros(m,m);
for i = 1:n
    MH = MH + inv(A(:,:,i))./n;
end
M = inv(MH);
tTot = cputime-tStart;

end