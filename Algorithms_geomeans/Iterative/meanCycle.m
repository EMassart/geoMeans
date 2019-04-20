function [ M, tTot, varargout] = meanCycle( A )
%MEANCYCLE(A) estimates the Karcher mean of a set of positive definite matrices with a
%Cycle algorithm. The Cycle algorithm is an iterative algorithm in which
%the matrices are represented as the nodes of a graph. At each iteration, a
%cycle is randomly chosen among these nodes, and the two nodes incident to each
%edge are replaced by their mean. The replacement of the nodes is immediately 
%performed, meaning that the replacement of a node in the beginning of the
%cycle will impact all the subsequent nodes.
%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%M is the estimation of the Karcher mean and tTot is the CPU time used by the function


tStart = cputime;
[~,~,n] = size(A);
D = zeros(n,n);
%tRecord = zeros(1,nIterMax);
nIterMax = 1000;
iter = 0;
count = 0;

%the computation of the distance separating the matrices is not formally
%required by the algorithm, therefore the CPU time required will be removed
%from the total 

tCumulStart = cputime;
for i = 1:n-1
    for j = i+1:n
        D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
    end
end
tCumulStop = cputime - tCumulStart;

while (max(max(D))>= 1e-8 && iter<nIterMax)
    
    iter = iter+1;
%     t = cputime;
    nodes = randperm(n);
 
    for i = 1:n-1
        i1 = nodes(i);
        i2 = nodes(i+1);
        A(:,:,i1) = MGeom(A(:,:,i1),A(:,:,i2),0.5);
        A(:,:,i2) = A(:,:,i1);
    end
    
    A(:,:,nodes(n)) = MGeom(A(:,:,nodes(n)),A(:,:,nodes(1)),0.5);
    count = count + n;
    A(:,:,nodes(1)) = A(:,:,nodes(n));
    
%     tRecord(iter) = (cputime-t)/n;
    
    tCumulStart = cputime;
    D = zeros(n,n);
    for i = 1:n-1
        for j = i+1:n
            D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
        end
    end
    tCumulStop = tCumulStop + cputime - tCumulStart;
end
    
M = A(:,:,1);
varargout{1} = count;
% tRecord = tRecord(1:iter);
% varargout{2} = tRecord;
if(iter == nIterMax),  fprintf('Maximum number of iterations reached \n'); end

tTot = cputime - tStart - tCumulStop;

end