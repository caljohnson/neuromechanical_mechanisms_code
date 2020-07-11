%understanding the mechanical coupling matrix D4^-1 (D4inv)

addpath('./src');

dims = [1:12 15 20 25]';
gridsz = 1;

%fixed chain length
for jj = 1:size(dims,1)
    dim = dims(jj);
    %second-difference matrix D2
    delX = 1/(gridsz*dim);
    e = ones(gridsz*dim,1);
    D2 = (1/delX^2).*spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
    D4 = D2*D2';
    D4inv = inv(full(D4));
    
    %max row sum (max coupling input to a given module)
    maxrowsum(jj) = max(sum(D4inv,2));
    minrowsum(jj) = min(sum(D4inv,2));
    
end

%fixed module length
for jj = 1:size(dims,1)
    dim = dims(jj);
    %second-difference matrix D2
    delX = 1/6;
    e = ones(gridsz*dim,1);
    D2 = (1/delX^2).*spdiags([e -2*e e], [0 1 2], gridsz*dim, gridsz*dim+2);
    D4 = D2*D2';
    D4inv = inv(full(D4));
   
    
    %max row sum (max coupling input to a given module)
    maxrowsum2(jj) = max(sum(D4inv,2));
    minrowsum2(jj) = min(sum(D4inv,2));

    
end

figure(1); clf;
semilogy(dims,maxrowsum,'bo','MarkerSize',15,'MarkerFaceColor', 'b'); hold on
semilogy(dims,maxrowsum2,'ro','MarkerSize',15,'MarkerFaceColor', 'r'); 
h = legend('fixed body length $L$', 'fixed module length $\ell$');
h.Interpreter = 'latex';
semilogy(dims,minrowsum,'bo','MarkerSize',15);%,'MarkerFaceColor', 'b'); hold on
semilogy(dims,minrowsum2,'ro','MarkerSize',15);%,'MarkerFaceColor', 'r'); 
shade(dims,minrowsum,dims,maxrowsum,'FillType',[1 2; 2 1],'Color','b')
shade(dims,minrowsum2,dims,maxrowsum2,'FillType',[1 2; 2 1],'Color','r')
xlabel('# of modules (N)'); 
ylabel({'range of mechanical', 'coupling strengths'});
set(gca, 'FontSize',30)