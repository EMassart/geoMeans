function [X1,m]=armijo_ls_spd(X0,grad,dir,alpha,beta,sigma,mmax,retraction,matrices, varargin)
%Authors : B. Jeuris, R. Vandebril and B. Vandereycken
%Implementation coming from : http://people.cs.kuleuven.be/~raf.vandebril/

% armijo_ls_spd: Armijo line search algorithm
%   [X1,m]=armijo_ls_spd(X0,grad,dir,alpha,beta,sigma,mmax,retraction,mat)
%   X1 is the point at greatest distance in the direction of tangent vector 
%   dir from the initial point X0, with gradient grad, that satisfies the 
%   Armijo condition determined by alpha, beta, and sigma (as in [1]). The 
%   output m is the power to which beta was taken and thus gives the 
%   decrease of the stepsize, of which the maximum is set to mmax. The 
%   argument retraction contains the desired type of retraction (as in the 
%   function retr) and the argument mat contains the given matrices in the 
%   mean.
%
%    [1] Optimization Algorithms on Matrix Manifolds, P.-A. Absil, R.
%        Mahony, R. Sepulchre, 2008


if nargin >= 10
    w = varargin{1};
else 
    w = ones(1,length(matrices));
end

cost0=karcher_cost(X0,matrices,w);

m=0;
X1=retr(X0,alpha*dir,retraction);
if (any(any(isnan(X1))))
    dist_decrease=-Inf;
    X1=X0;
else
    dist_decrease=cost0-karcher_cost(X1,matrices,w);
end

inpro_grad_dir=-sigma*alpha*inpro_spd(grad,dir,X0);

while ((dist_decrease<beta^m*inpro_grad_dir || any(eig(X1)<=0)) && m<=mmax)
    m=m+1;
    X1=retr(X0,beta^m*alpha*dir,retraction);
    if (any(any(isnan(X1))))
        dist_decrease=-Inf;
        X1=X0;
    else
        dist_decrease=cost0-karcher_cost(X1,matrices,w);
    end
end