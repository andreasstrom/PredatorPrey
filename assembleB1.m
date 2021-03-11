function [A,b] = assembleB1(p,t,f)
N = size(p,2);
A = sparse(N,N);
b = zeros(N,1);
for K = 1:size(t,2)
    nodes = t(1:3,K); 
    x = p(1,nodes);
    y = p(2,nodes);
    coeff = [[1 1 1]' x' y']\eye(3);
    area_K = polyarea(x,y);
    AK = (coeff(2,:)'*coeff(2,:) + coeff(3,:)'*coeff(3,:)).*area_K;
    bK = f(mean(x),mean(y)).*area_K./3.*ones(3,1);
    A(nodes,nodes) = A(nodes,nodes) + AK;
    b(nodes) = b(nodes) + bK;
end