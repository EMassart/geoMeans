function [ M, tTot, varargout ] = meanProgressExpansion( A, varargin )
%MEANPROGRESSEXPANSION(A) estimates the Karcher mean using the ProgressiveExpansion algorithm.
%This algorithm is a kind of Repetitive merging algorithm in which the 
%underlying merging algorithm is chosen to follow a progressive expansion.

%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%nIterMax is the number of iterations allowed to the algorithm
%M is the estimation of the Karcher mean and tTot is the CPU time used by
%the function

if nargin >= 2
    nIterMax = varargin{1};
else
    nIterMax = 1;                               %by default, only one iteration of the algorithm is performed, and then 
                                                %the arithmetic mean of the matrices obtained is computed
end

tStart = cputime;
[m,~,n] = size(A);
D = zeros(n,n); 
D_geod = zeros(n,n);
B = zeros(m,m,n);
count = 0;                                      %records the total number of means of two matrices computed 
iter = 0;
dMaxRecord = zeros(1,nIterMax);                 %records the maximum distance separating the new iterates

%the computation of the distance separating the matrices is not formally
%required by the algorithm, therefore the CPU time required will be removed
%from the total 
tRemStart = cputime;
for i = 1:n-1
    for j = i+1:n
        D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
        D_geod(i,j) = dist(A(:,:,i),A(:,:,j));
    end
end
dMaxRecord(iter+1) = max(max(D_geod));
tRem = cputime-tRemStart;


while (max(max(D))>= 1e-8 && iter<nIterMax)
    
    iter = iter+1;
    
    for i = 1:n
        G = mod(i:1:i+n-1,n);                                 %order in which the nodes will be merged
        G(G==0) = n;
        [B(:,:,i)] = meanProgressExpansion_sub(A,G);
        count = count + n-1;
    end
    
    A = B;
    
    tRemStart = cputime;
    for i = 1:n-1
        for j = i+1:n
            D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
            D_geod(i,j) = dist(A(:,:,i),A(:,:,j));
        end
    end  
    dMaxRecord(iter+1) = max(max(D_geod));
    tRem = tRem + cputime-tRemStart;
end
   

%Final step
M = sum(A,3)./n;
dMaxRecord = dMaxRecord(1:iter+1);
varargout{1} = count;
varargout{2} = dMaxRecord;

% if(iter == nIterMax),  fprintf('MeanProgressExpansion : Maximum number of iterations reached \n'); end

tTot = cputime-tStart-tRem;
end