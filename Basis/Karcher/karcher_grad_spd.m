function [grad,logsum]=karcher_grad_spd(X,matrices,varargin)
%Authors : B. Jeuris, R. Vandebril and B. Vandereycken
%Implementation coming from : http://people.cs.kuleuven.be/~raf.vandebril/
% karcher_grad_spd: The gradient of the Karcher cost function
%   [grad,logsum]=karcher_grad_spd(X,mat)
%   grad is the gradient of the Karcher costfunction, determined by the 
%   matrices in mat, at the point X when the manifold is endowed with inner 
%   product inpro_spd. The output logsum contains the sum of all terms 
%   logm(matrices{i}\X), which can be used in further constructions.

%les poids sont normalisés

if nargin >= 3
    w = varargin{1};
else 
    w = ones(1,length(matrices));
end

logsum=zeros(size(X,1));

for i=1:length(matrices)
    logsum=logsum + w(i)*logm(matrices{i}\X);
end

grad=2*X*logsum;

grad=(grad+grad')/2;