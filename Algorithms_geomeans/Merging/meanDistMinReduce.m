function [ A, tTot, varargout] = meanDistMinReduce( A )
%DISTMINREDUCE(A) estimates the Karcher mean using the DistMax merging algorithm.
%This algorithm consists in merging successively the matrices, until there
%is only one matrix left. The pair of matrices merged are made of the two 
%closest matrices.

%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%M is the estimation of the Karcher mean and tTot is the CPU time used by
%the function

tStart = cputime;
[~,~,n] = size(A);                          % number of matrices
w = ones(1,n);                              % weights to give to each node
ncur = n;                                   % current number of nodes
D = zeros(n,n);
iter = 0;
% tRecord = zeros(1,n-1);

for i = 1:n-1
    for j = i+1:n
        [D(i,j),~] = dist(A(:,:,i),A(:,:,j));          
    end
end

while (ncur > 1)
    iter  = iter + 1;
%     t = cputime;

    %choice of the pair of matrices (the two closest ones)
    [v,i1vect] = min(D+tril(ones(ncur)+NaN));
    [~,i2] = min(v);
    i1 = i1vect(i2);
    ivect = sort([i1 i2]);
    i1 = ivect(1);
    i2 = ivect(2);

    A(:,:,i1) = MGeom(A(:,:,i1),A(:,:,i2),w(i2)/(w(i1)+w(i2)));              %weighted mean of the two matrices
    A(:,:,i2) = [];
    
    %update of the weights, of the distance matrix and of the number of nodes
    w(i1) = w(i1)+w(i2);
    w(i2) = [];
    
    D(i2,:) = [];
    D(:,i2) = [];
    
    for i = i1+1:ncur-1
        [D(i1,i),~] = dist(A(:,:,i1),A(:,:,i));
    end
    for i = 1:i1-1
        [D(i,i1),~] = dist(A(:,:,i),A(:,:,i1));
    end
%     tRecord(iter) = cputime-t;
    ncur = ncur - 1;
    
end

tTot = cputime-tStart;
end