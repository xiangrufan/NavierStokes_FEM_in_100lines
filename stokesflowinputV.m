addpath distmesh
clear all
% [p,t,b] from distmesh tool
% make sure your matlab path includes the directory where distmesh is installed.
clear all
fd=@(p)  drectangle(p,-1,1,-1,1);
[p,t]=distmesh2d(fd,@huniform,0.1,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);
b=unique(boundedges(p,t));
MS=mesh(p,t);
N=size(p,1);T=size(t,1);
K=sparse(N,N);
BGLS=K;Z=K;Kpx=K;Kpy=K;Kxp=K;Kyp=K;
BGLS2=K; BGLS2x=K;BGLS2y=K;BGLSx=K;
F=zeros(N,1);  T1=F;T2=F; EmptyCol=F;
sv=0.01;
for e=1:T
    nodes=t(e,:);
    Pe=[ones(3,1),p(nodes,:)]; % 3 by 3 matrix with rows=[1 xcorner ycorner]
    Area=abs(det(Pe))/2;
    C=inv(Pe);
    grad=C(2:3,:);Ke=sv*Area*grad'*grad;
    a1=C(2,1);a2=C(2,2);a3=C(2,3);
    b1=C(3,1);b2=C(3,2);b3=C(3,3);
    tau=1/3*Area/(4*sv);
    KpxL= [1 1 1]'*[a1,a2,a3]*Area/3;
    KpyL =  [1,1,1]'*[b1,b2,b3]*Area/3;
    KxpL=  [a1,a2,a3]*[1 1 1]'*Area/3;
    KypL =  [b1,b2,b3]*[1,1,1]'*Area/3;
    pos=p(nodes,:);
    length1=(norm(pos(2,:)-pos(3,:)));
    length2=(norm(pos(1,:)-pos(3,:)));
    length3=(norm(pos(1,:)-pos(2,:)));
    lapu=[(a1-(a2+a3)/2)/length1,(a2-(a1+a3)/2)/length2,(a3-(a2+a1)/2)/length3];
    K(nodes,nodes)=K(nodes,nodes)+Ke;
    BGLS(nodes,nodes)=BGLS(nodes,nodes)-tau*grad'*grad;
    BGLSx(nodes,nodes)=BGLSx(nodes,nodes)-tau*lapu'*lapu;
    Kpx(nodes,nodes)=Kpx(nodes,nodes)+KpxL;
    Kpy(nodes,nodes)=Kpy(nodes,nodes)+KpyL;
    Kxp(nodes,nodes)=Kxp(nodes,nodes)+KxpL;
    Kyp(nodes,nodes)=Kyp(nodes,nodes)+KypL;
    F(nodes)=F(nodes);
    T1(nodes)=T1(nodes)+grad(1,:)';
    T2(nodes)=T2(nodes)+grad(2,:)';
    AA=1 ;
    
    BGLS2x(nodes,nodes)=BGLS2x(nodes,nodes)+AA*tau*[a1,a2,a3]'* [a1,a2,a3];
    BGLS2y(nodes,nodes)=BGLS2y(nodes,nodes)+AA*tau*[b1,b2,b3]'* [b1,b2,b3];
    
end

%Finding BC nodes
bc1=MS.rec_selector([-1.05,-1.05],[-0.95,1.05]);
bc2=MS.rec_selector([0.95,-1.05],[1.05,1.05]);
WallBC=  setdiff(b,bc2);

K(WallBC,:)=0; % Eliminate freedom in Essential BC
Kpx(WallBC,:) = 0;
Kpy(WallBC,:) =0;
K(WallBC,WallBC)=speye(length(WallBC),length(WallBC));
BGLS2x(WallBC,:) =0;
BGLS2y(WallBC,:) =0;

Kpp= BGLS;
xx=bc2(2) ;
Kpp(xx,:) = 0;
Kpp(xx,xx)=1;
Kxp(xx,:)=0; Kyp(xx,:)=0;
F1= F;   F2= F; F3= F;
F1(bc1)= 0.1*(1-abs(p(bc1,2))).^2;
BigK= [K+BGLS2x,Z,Kpx; Z,K+BGLS2y,Kpy;Kxp,Kyp,BGLS ];
BigF=[F1;F2;F3];
U=BigK\BigF;

figure
trisurf(t,p(:,1),p(:,2),0*p(:,1),U(1:N),'edgecolor','k','facecolor','interp');
view(2),axis([-1 1 -1 1]),axis equal,colorbar

figure
trisurf(t,p(:,1),p(:,2),0*p(:,1),U(2*N+1:3*N),'edgecolor','k','facecolor','interp');
view(2),axis([-1 1 -1 1]),axis equal,colorbar

figure
quiver(p(:,1),p(:,2),U(1:N),U(N+1:2*N));