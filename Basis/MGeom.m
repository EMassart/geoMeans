function [ m ] = MGeom(A,B,t)
%Geometric mean of two matrices

s = sqrtm(A);
m = s\B/s;
m = m^t;
m = s*m*s;
end