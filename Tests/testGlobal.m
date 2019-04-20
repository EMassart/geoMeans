function [] = testGlobal( str )
%Gives an overview of the global performances of the different algorithms.
%n = number of matrices per set
%m = size of the matrices
% WARNING : the execution of this code may take a few hours.

n = 5;     %number of matrices in each set
k = 1000;  %number of sets used in the test
m = 10;     %size of the matrices

list_of_methods = {'cheap',
    'meanHarmonic',
    'meanArithmetic',
    'meanRandomReduce',
    'meanDistMinReduce',
    'meanDistMaxReduce',
    'meanTree',
    'meanProgressExpansion',
    'meanIterativeRandom',
    'meanIterativeIdealMap',
    'meanIterativeIdealMapCheap',
    'meanStar',
    'meanCycle',
    'Crude',
    'CrudeImproved',
    'CrudeIterative'};

line_types = {'*b','*m','*k','ob','om','og','dk','dm','sr','sg','sk','*r','*g','^b','^k','^r'};

l = length(list_of_methods);
tTot_rec = zeros(l,k);
tTotM = zeros(l,1);
dist_rec = zeros(l,k);
distM = zeros(l,1);
d_ref = zeros(1,k);

for indx = 1:k
    fprintf('Set number %d\n',indx);
    A = zeros(m,m,n);
    Acell = cell(1,n);
    for jLoc = 1:n
        A(:,:,jLoc) = defPos4(m);
        Acell{jLoc} = A(:,:,jLoc);
    end
    [meanKarcher,~ ] = karcher_sd_spd(sum(A,3)./n,'approx2',Acell);
    
   %runs the methods and computes the distances between the solutions
   %returned and the Karcher mean computed with a steepest descent
   %algorithm
    for i = 1:l
     [M, tTot_rec(i,indx)] = eval([list_of_methods{i},'( A );']);
     dist_rec(i,indx) = dist(M,meanKarcher);
    end
    d_ref(indx) = dist(A(:,:,1),meanKarcher);    
end

%averages the values obtained over the sets
for i = 1:l
    tTotM(i) = mean(tTot_rec(i,:));
    distM(i) = mean(dist_rec(i,:)./d_ref);
end
save(str);

%plot
figure;
for idx_method  = 1:l
  plot(tTotM(idx_method),distM(idx_method),line_types{idx_method}); hold on;
end
xlabel('CPU time');
ylabel('Error E_{rel}');
axis([-0.005 0.06 -0.1 1.3]);
legend(list_of_methods)
% additive axes
vert = graph2d.constantline(0, 'LineStyle',':', 'Color',[.5 .5 .5]);
changedependvar(vert,'x');
hori = graph2d.constantline(0, 'LineStyle',':', 'Color',[.5 .5 .5]);
changedependvar(hori,'y');


end