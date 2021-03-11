function [A,M] = assembleB2(p,t)
N = size(p,2);
A = sparse(N,N);
M = sparse(N,N); % allocate mass matrix
for K = 1:size(t,2)
    nodes = t(1:3,K);
    x = p(1,nodes);
    y = p(2,nodes);
    coeff = [[1 1 1]' x' y']\eye(3);
    area_K = polyarea(x,y);
    AK = (coeff(2,:)'*coeff(2,:) + coeff(3,:)'*coeff(3,:)).*area_K;
    A(nodes,nodes) = A(nodes,nodes) + AK;
    MK = [2 1 1; 1 2 1; 1 1 2]./12*area_K; % element mass matrix
    M(nodes,nodes) = M(nodes,nodes) + MK; % add element masses to M
end