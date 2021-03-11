close all
f = @(x,y)(8*pi^2.*sin(2*pi*x)*sin(2*pi*y));
uex = @(x,y)(sin(2*pi*x).*sin(2*pi*y));
h = 1./[2 4 8 16 32];
EnE = zeros(length(h),1);
pvec = cell(1,5);
uvec = cell(1,5);
for i = 1:length(h)
    [p, e, t]=initmesh(@circleg,'hmax',h(i));
    [A,b] = assembleB1(p,t,f);
    I = eye(length(p));
    A(e(1,:),:) = I(e(1,:),:);
    b(e(1,:)) = uex(p(1,e(1,:)),p(2,e(1,:)));
    u = A\b;
    err = uex(p(1,:),p(2,:))-u';
    EnE(i) = sqrt(err*A*err');  %Energy norm Error
    uvec{i} = u;                %for plotting
    pvec{i} = p;
end

%plotting
figure(1)
convRate = polyfit(log(h),log(EnE),1);
loglog(h,EnE,'r',h,h.^convRate(1),'b')
legend('EnE',strcat('h_{max}^p, p=',num2str(convRate(1))))
title('EnE & h_{max}')
xlabel('h_{max}')
figure(2)
%plot3(pvec{1}(1,:),pvec{1}(2,:),uvec{1},'*')
[X,Y] = meshgrid(linspace(max(pvec{1}(1,:)),min(pvec{1}(1,:)),length(pvec{1}(1,:))),...
    linspace(max(pvec{1}(2,:)),min(pvec{1}(2,:)),length(pvec{1}(2,:))));
mesh(griddata(pvec{1}(1,:),pvec{1}(2,:),uvec{1},X,Y))
figure(3)
%plot3(pvec{5}(1,:),pvec{5}(2,:),uvec{5},'*')
%too computationally heavy to plot the full matrix from griddata() using 
%mesh() for my laptop, hence the number of elements in linspace() has been
%decreased, since it only affects the plotting and not the calculations
[X,Y] = meshgrid(linspace(max(pvec{5}(1,:)),min(pvec{5}(1,:)),length(pvec{5}(1,:))./5),...
    linspace(max(pvec{5}(2,:)),min(pvec{5}(2,:)),length(pvec{5}(2,:))./5));
mesh(griddata(pvec{5}(1,:),pvec{5}(2,:),uvec{5},X,Y))
