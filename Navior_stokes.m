addpath distmesh
clear all

fd=@(p) ddiff(drectangle(p,-1,1,-1,1),dcircle(p,0,0,0.5));
fh=@(p) 0.05+0.3*dcircle(p,0,0,0.5);
[p,t]=distmesh2d(fd,fh,0.05,[-1,-1;1,1],[-1,-1;-1,1;1,-1;1,1]);


b=unique(boundedges(p,t));

MS=mesh(p,t);

N=size(p,1);T=size(t,1);

K=sparse(N,N);
BGLS=K;
Z=K;
Kpx=K;
Kpy=K;
Kxp=K;
Kyp=K;


F=zeros(N,1);
sv=0.01;
for e=1:T
    nodes=t(e,:);
    Pe=[ones(3,1),p(nodes,:)]; % 3 by 3 matrix with rows=[1 xcorner ycorner]
    Area=abs(det(Pe))/2;
    C=inv(Pe);
    
    grad=C(2:3,:);Ke=sv*Area*grad'*grad;
    
    Fe=0;
    
    a1=C(2,1);a2=C(2,2);a3=C(2,3);
    b1=C(3,1);b2=C(3,2);b3=C(3,3);
    
    tau=1*(1/3)*Area/(4*sv);
    
    KpxL= [1 1 1]'*[a1,a2,a3]*Area/3;
    KpyL =  [1,1,1]'*[ b1,b2,b3]*Area/3;
    KxpL=  [1,1,1]'*[a1,a2,a3]*Area/3;
    KypL =   [1,1,1]'*[b1,b2,b3]*Area/3;
    
    K(nodes,nodes)=K(nodes,nodes)+Ke;
    BGLS(nodes,nodes)=BGLS(nodes,nodes)+tau*0.001*grad'*grad;
    Kpx(nodes,nodes)=Kpx(nodes,nodes)+KpxL;
    Kpy(nodes,nodes)=Kpy(nodes,nodes)+KpyL;
    Kxp(nodes,nodes)=Kxp(nodes,nodes)+KxpL;
    Kyp(nodes,nodes)=Kyp(nodes,nodes)+KypL;
    
    F(nodes)=F(nodes)+Fe;
    
end


%Finding BC nodes
bc1=MS.rec_selector([-1.05,-1.05],[-0.95,1.05]);
bc2=MS.rec_selector([0.95,-1.05],[1.05,1.05]);
WallBC=  setdiff(b,bc2);

K(WallBC,:)=0; % Eliminate freedom in Essential BC
Kpx(WallBC,:) = 0;
Kpy(WallBC,:) =0;
K(WallBC,WallBC)=speye(length(WallBC),length(WallBC));
% I do not know why they have K(:,b)=0. This would invalidate the laplacian


Kpp= BGLS;
xx=bc2 ;
Kpp(xx,:) = 0;
Kpp(xx,xx)=1;



% select the left most reactangle BC
F1= F;
F3= F;
F2=F;
%similar to this but below code is not possible

F1(bc1)= 100*(1-abs(p(bc1,2))).^2;
%  F1(bc1)=1;
%  BigK= [K,Z,Kpx; Z,K,Kpy;-Kxp,-Kyp,BGLS]; % Wrong, Because in row 3, P virtual is not restrained in place other then Bc1 and Bc2
BigK= [K,Z,Kpx ; Z,K,Kpy ;Kxp,Kyp,Kpp];
BigF=[F1;F2;F3];
%%
U=BigK\BigF;  % stokes solution
rho=1;


it=0;
err=100;
while err>1e-8 & it<50
    
    ududx=zeros(N,N);
    vdudy=zeros(N,N);
    udvdx=zeros(N,N);
    vdvdy=zeros(N,N);
    dudx=zeros(N,N);
    dudy=zeros(N,N);
    dvdx=zeros(N,N);
    dvdy=zeros(N,N);
    
    for  e=1:T
        nodes=t(e,:);
        Pe=[ones(3,1),p(nodes,:)]; % 3 by 3 matrix with rows=[1 xcorner ycorner]
        Area=0.02*abs(det(Pe))/2;
        C=inv(Pe);
        a1=C(2,1);a2=C(2,2);a3=C(2,3);
        b1=C(3,1);b2=C(3,2);b3=C(3,3);
        u1=U(nodes(1));u2=U(nodes(2));u3=U(nodes(3));
        v1=U(N+nodes(1));v2=U(N+nodes(2));v3=U(N+nodes(3));
        
        ududx(nodes,nodes)=ududx(nodes,nodes)+Area*[u1,u2,u3]'*C(2,:);
        vdudy(nodes,nodes)=vdudy(nodes,nodes)+Area*[u1,u2,u3]'*C(3,:);
        udvdx(nodes,nodes)=udvdx(nodes,nodes)+Area*[v1,v2,v3]'*C(2,:);
        vdvdy(nodes,nodes)=vdvdy(nodes,nodes)+Area*[v1,v2,v3]'*C(3,:);
        dudx(nodes,nodes)=dudx(nodes,nodes)+Area*[1,1,1]'*C(2,:);
        dudy(nodes,nodes)=dudy(nodes,nodes)+Area*[1,1,1]'*C(3,:);
        dvdx(nodes,nodes)=dvdx(nodes,nodes)+Area*[1,1,1]'*C(2,:);
        dvdy(nodes,nodes)=dvdy(nodes,nodes)+Area*[1,1,1]'*C(3,:);
    end
    ududx(WallBC,:) = 0;
    vdudy(WallBC,:) = 0;
    udvdx(WallBC,:) = 0;
    vdvdy(WallBC,:) = 0;
    NonlinearRes= BigK-[ududx,vdudy,Z;udvdx,Z+vdvdy,Z;Z,Z,Z];
    Jacobian = BigK+ [diag(dudx*U(1:N)),diag(dudy*U(1:N)),Z;diag(dvdx*U(N+1:2*N)),diag(dvdy*U(N+1:2*N)),Z;Z,Z,Z];
    % U can add an damper like damping=0.2; U = Uold*(damping) + (U-Uold)*(1-damping)
    Uold=U;
    U=NonlinearRes\BigF;
    it=it+1;
    err= norm(Uold-U);
end




%   figure
% trisurf(t,p(:,1),p(:,2),0*p(:,1),U(1:N),'edgecolor','k','facecolor','interp');
% view(2),axis([-1 1 -1 1]),axis equal,colorbar
%
figure
trisurf(t,p(:,1),p(:,2),0*p(:,1),(U(1:N).^2+U(N+1:2*N).^2).^0.5,'edgecolor','k','facecolor','interp');
view(2),axis([-1 1 -1 1]),axis equal,colorbar
figure
quiver(p(:,1),p(:,2),U(1:N),U(N+1:2*N));
%
%
% figure
% trisurf(t,p(:,1),p(:,2),0*p(:,1),U(N+1:2*N),'edgecolor','k','facecolor','interp');
% view(2),axis([-1 1 -1 1]),axis equal,colorbar
% sum( BGLS*U(2*N+1:3*N))
% figure
% trisurf(t,p(:,1),p(:,2),0*p(:,1),U(2*N+1:3*N),'edgecolor','k','facecolor','interp');
% view(2),axis([-1 1 -1 1]),axis equal,colorbar