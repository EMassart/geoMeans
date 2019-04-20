function [ B ] = meanProgressExpansion_sub( A, G )
%MEANPROGRESSEXPANSION_SUB(A) estimates the Karcher mean using the ProgressiveExpansion
%merging algorithm. 
%This algorithm consists in merging successively the matrices, until there
%is only one matrix left. The pair of matrices merged are chosen so that
%the mean computed at one iteration is automatically merged during the next
%iteration.

%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%B is the estimation of the Karcher mean

[~,~,n] = size(A);
B = A(:,:,G(1));
alpha = 1/2;                                    %weight to give to the new node (updated during the execution of the algorithm)
% count = 0;                                      %records the total number of means of two matrices computed

for i = 2:n
    B = MGeom(B,A(:,:,G(i)),alpha);
    alpha = 1/(i+1);
%     count = count + 1;
end
    
end