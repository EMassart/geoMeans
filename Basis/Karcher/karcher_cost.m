function y=karcher_cost(X,matrices,varargin)
%Authors : B. Jeuris, R. Vandebril and B. Vandereycken
%Implementation coming from : http://people.cs.kuleuven.be/~raf.vandebril/

% karcher_cost: Evaluation of the Karcher cost function
%   karcher_cost(X,mat) returns the value of the Karcher cost function, 
%   determined by the matrices given in mat, evaluated in the point X.

if nargin >= 3
    w = varargin{1};
else 
    w = ones(1,length(matrices));
end

y=0;
sqrtX=sqrtm(X);
for i=1:length(matrices)
    arg=sqrtX\matrices{i}/sqrtX;
    if (norm(imag(eig(arg)),'fro')>1e-15)
        y=Inf;
        break;
    elseif (any(real(eig(arg))<0))
        y=Inf;
        break;
    end
    y=y+w(i)*norm(logm(arg),'fro')^2;
end