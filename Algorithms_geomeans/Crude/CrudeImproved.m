function [M_res,tTot] = CrudeImproved(A,varargin)
%CRUDEIMPROVED(A,n_cl) : Adaptation of the Crude mean that was proposed in [1]
%to get a more accurate estimation of the Karcher mean by splitting
%the input matrices into sets regarding their condition number.
%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%M_res is the estimation of the Karcher mean and tTot is the CPU time used by the function

%References : 
%[1] Ben Jeuris and Raf Vandebril. Geometric mean algorithms based on harmonic and
%arithmetic iterations. In Geometric Science of Information, pages 785–793. Springer,
%2013.

if nargin >= 2
    n_cl = varargin{1};                                 %number of sets among which the input matrices will be split
else
    n_cl = 2;
end

tStart = cputime;
[m,~,n] = size(A);
lco = zeros(n,1);                                   
for i = 1:n                                             
    lco(i) = log10(cond(A(:,:,i)));
end

lcoRange = linspace(0,8,n_cl+1);                        %defines the range of condition numbers associated to each set (making the assumption 
                                                        %that the logarithms of the condition numbers are mostly located between 10^0 and 10^8).
                                                        %The matrices having a condition number whose logarithm is lower than 10^0 will be added 
                                                        %to the first set and the ones having a condition number bigger than 10^8 will be added to 
                                                        %the last one.                                                       

ind = 0;                                                %will record the effective number of sets (the number of non-empty sets)      
sum_indices = 0;
M = zeros(m,m,n);                                       %will contain the Crude mean, computed for each non-empty set 
                                                        %(there are at most n sets because the total number of matrices is n).
                                                        
w = zeros(1,n);                                         %number of matrices belonging to each non-empty set

% X0 = zeros(m,m);                                      %only required if the final step is performed with a steepest descent algorithm

for i = 2:length(lcoRange)
    if i<length(lcoRange)                               
        indices = find(lco >= lcoRange(i-1) & lco < lcoRange(i));
    else 
        indices = find(lco >= lcoRange(i-1));
    end
        
    if ~isempty(indices)
        ind = ind + 1;
        M(:,:,ind) =  Crude( A(:,:,indices));
        w(ind) = length(indices);
%         X0 = X0 + w(ind)*M(:,:,ind);
    end
    sum_indices = sum_indices + length(indices);
end

% X0 = X0./sum(w);
M = M(:,:,1:ind);                                        
w = w(1:ind);

%final step : combines the Crude means of each set
if ind==2                                                
    M_res = MGeom(M(:,:,1),M(:,:,2),w(2)/(w(1)+w(2)));
else
    [ M_res,~] = meanDistMaxReduce( M );
    
% or with the Karcher mean to perform the last step (considerably more costly!)
%     MCell = cell(1,length(w));
%     for j = 1:length(w);
%         MCell{j} = M(:,:,j);
%     end
%     [ M_res,~ ] = karcher_sd_spd(X0,'approx2',MCell,w);
end
    
tTot = cputime-tStart;

end
