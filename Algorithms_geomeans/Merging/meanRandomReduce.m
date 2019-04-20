function [ A, tTot, varargout ] = meanRandomReduce( A )
%MEANRANDOMREDUCE(A) estimates the Karcher mean using the random merging algorithm.
%This algorithm consists in merging successively the matrices, until there
%is only one matrix left. The pair of matrices merged are chosen randomly
%among the remaining matrices.

%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%M is the estimation of the Karcher mean and tTot is the CPU time used by
%the function

tStart = cputime;
[~,~,n] = size(A);                                      %number of matrices
w = ones(1,n);                                          % weights to give to each node
ncur = n;                                               % current number of nodes
iter = 0;
% tRecord = zeros(1,n-1);


while (ncur > 1)
    
    iter  = iter + 1;
%     t = cputime;
    per = randperm(ncur);
    ivect = sort([per(1) per(2)]);
    i1 = ivect(1);
    i2 = ivect(2);

    A(:,:,i1) = MGeom(A(:,:,i1),A(:,:,i2),w(i2)/(w(i1)+w(i2)));         %weighted mean of the two matrices
    A(:,:,i2) = [];
    w(i1) = w(i1)+w(i2);
    w(i2) = [];
    
%     tRecord(iter) = cputime-t;   
    ncur = ncur - 1;
    
end

% varargout{1} = tRecord;
tTot = cputime-tStart;

end