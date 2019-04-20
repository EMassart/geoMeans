function [ A ] = defPosCond( n, cond, varargin )
%DEFPOSCOND(n,cond,varargin) : generates a symmetric positive definite matrix 
%of size n, norm equal to one and condition number about equal to cond, 
%in the same way as in [1].
%A is the matrix obtained
%varargin = 1 -> an eigenvalue considerably smaller than the others
%varargin = 2 -> an eigenvalue considerably larger than the others
%varargin = 3 -> the first half of the eigenvalues are small while the 
%                others are considerably larger

%[1] : B. Jeuris, R. Vandebril, B. Vandereycken : A survey and comparison of
%contemporary algorithms for computing the matrix geometric mean.

if nargin <= 2
    version = 1;
else 
    version = varargin{1};
end

if version == 1
%     disp('Version1');
    [Q,~] = qr(rand(n));
    D = diag([rand(1,n-1)+1,10^(-cond)]);
    A = Q*D*Q';
elseif version == 2
%     disp('Version2');
    [Q,~] = qr(rand(n));
    D = diag([1,10^(-cond)+10^(-cond)*rand(1,n-1)]);
    A = Q*D*Q';
elseif version == 3  
%     disp('Version3');
    [Q,~] = qr(rand(n));
    if mod(n,2) == 0
         D = diag([1+rand(1,n/2),10^(-cond)+10^(-cond)*rand(1,n/2)]);
    else
         D = diag([1+rand(1,floor(n/2)+1),10^(-cond)+10^(-cond)*rand(1,floor(n/2))]);
    end
    A = Q*D*Q';
end

nor = norm(A,'fro');
A = A./nor;

end




