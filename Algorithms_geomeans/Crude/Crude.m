function [ M, tTot] = Crude( A)
%CRUDE(A) computes the Crude mean of a set of matrices.
%A is a 3D array containing the matrices along the third dimension
%M is the Crude mean and tTot is the CPU time used by the function

%References : 
%Ben Jeuris and Raf Vandebril. Geometric mean algorithms based on harmonic and
%arithmetic iterations. In Geometric Science of Information, pages 785–793. Springer,
%2013.

tStart = cputime;
[m,~,n] = size(A);
IA = zeros(m,m,n);                      

for i = 1:n
    IA(:,:,i) = inv(A(:,:,i));
end

MA = sum(A,3)./n;                       %arithmetic mean
MH = inv(sum(IA,3)./n);                 %harmonic mean
M = MGeom(MA,MH,0.5);                   

tTot = cputime-tStart;

end

