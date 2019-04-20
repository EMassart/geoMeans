function [ M, tTot, varargout] = meanIterativeIdealMapCheap( A )
%MEANITERATIVEIDEALMAPCHEAP(A) estimates the Karcher mean using an algorithm based on
%the Iterative algorithm proposed in [1].
%The cycle is chosen at each iteration in order to maximize the
%communication speed between the nodes.
%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%M is the estimation of the Karcher mean and tTot is the CPU time used by
%the function

%References : 
%[1] Miklós Pálfia. A multivariable extension of two-variable matrix means. SIAM Journal
%on Matrix Analysis and Applications, 32(2):385–393, 2011.

tStart = cputime;
[~,~,n] = size(A);
C = zeros(n,n);                                     %will contain the weights of the influence of each input matrix in the current nodes
                                                    %(used to compute the new cycle at each iteration)
nIterMax = 1000;
iter = 0;
count = 0;                                          %records the number of means of two matrices evaluated


%the computation of the distance separating the matrices is not formally
%required by the algorithm, therefore the CPU time required will be removed
%from the total 
tCumulStart = cputime;
D = zeros(n,n);
for i = 1:n-1
    for j = i+1:n
        D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
    end
end
D = D+D';
tCumulStop = cputime - tCumulStart;

while( max(max(D))>=1e-8 && iter<nIterMax)
    
    iter = iter+1;
    %computation of the cycle (the assumption is made here that the
    %sequence of cycles is precomputed, therefore the required time is not recorded)
    if (iter==1)
        
        tCumulStart = cputime;  
        G = 1:n;                   
        for i = 1:n-1                       %C is firstly defined as the transposed of the incidence matrix
            C(i,G(i)) = 1;
            C(i,G(i+1)) = 1;
        end
        C(n,G(n)) = 1;
        C(n,G(1)) = 1;
        tCumulStop = tCumulStop + cputime - tCumulStart;
        
    else
        
        tCumulStart = cputime; 
        D_cycle = triu(squareform(pdist(C)));
        G = idealMap(D_cycle);                      %computes the new cycle
        Cbis = C;                                   %updates the matrix C 
        for i = 1:n-1
            C(i,:) = (Cbis(G(i),:) + Cbis(G(i+1),:))./2;  
        end
        C(n,:) = (Cbis(G(n),:) + Cbis(G(1),:))./2;
        tCumulStop = tCumulStop + cputime - tCumulStart;
    end
    
    B = A;
    for i = 1:n-1
        A(:,:,i) = MGeom(B(:,:,G(i)), B(:,:,G(i+1)), 0.5);
    end
    A(:,:,n) = MGeom(B(:,:,G(n)), B(:,:,G(1)), 0.5);
    
    count = count + n;
    
    tCumulStart = cputime;
    D = zeros(n,n);
    for i = 1:n-1
        for j = i+1:n
            D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
        end
    end
    D = D+D';
    tCumulStop = tCumulStop + cputime - tCumulStart;

end

M = A(:,:,1);
varargout{1} = iter;
if(iter == nIterMax),  fprintf('Maximum number of iterations reached \n'); end

tTot = cputime-tStart-tCumulStop;
end





