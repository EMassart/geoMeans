function [ A ] = defPos4( n )
%DEFPOS4(n) : generates a symmetric positive definite matrix of size n and
%norm equal to one, according to the Wishart distribution.
%A is the matrix obtained

A = randn(n, n);
A = A'*A;

nor = norm(A,'fro');
A = A./nor;

lambdas = eig(A);
e = find(lambdas<=0);

if(size(e)~= 0)
    disp('The matrix is not positive definite');
end

end




