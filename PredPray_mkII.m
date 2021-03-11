close all
T=2; a=4; d=.01; tn=.01; ti=0:tn:T; h=1./[5 20 40]; h=h(1);
[p, e, t]=initmesh(@circleg,'hmax',h);
[A,M] = assembleB2(p,t);
u = zeros(length(p),length(ti));
u(:,1) = 1+20*rand(length(p),1);
pRate = zeros(1,length(ti));
for K = 1:length(t)
    nodes = t(1:3,K);
    pRate(1) = pRate(1)+polyarea(p(1,nodes),p(2,nodes))./3.*...
        (sum(u(nodes,1),'all'));
end
for n=2:length(ti)
    u(:,n) = (M/tn-M/2+d*A/2)\((M/tn+M/2-d*A/2)*u(:,n-1)...
        -M*(u(:,n-1).^2+u(:,n-1)./(u(:,n-1)+a)));
    for K = 1:length(t)
        nodes = t(1:3,K);
        pRate(n) = pRate(n)+polyarea(p(1,nodes),p(2,nodes))./3.*...
            sum(u(nodes,n),'all');
    end
end
figure(1)
plot(ti,pRate)
xlabel('time')
ylabel('population rate')
title(strcat('h_{max}=',num2str(h)))
figure(2)
pdesurf(p,t,u(:,end))
title(strcat('Solution for t=2 and h_{max}=',num2str(h)))