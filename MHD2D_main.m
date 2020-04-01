%% MHD2D_main
% 二维可压缩磁流体模拟（用MATLAB语言复原工作）
% 陈宗华(Zong-hua Chen, 814484233@qq.com,15307756701) 
% 玉林师范学院（东校区）物理与电信工程学院
% 广西玉林市教育东路1303号，537000
% 2019-08
%%
close;clear;clc;
global dx dz gamma nu eta dt Bx Bz
nx=33;
nz=32;
nt=5000;
dt=0.05;  % 要满足 dt < min(dx,dz)/[sqrt(1.0+0.5*gamma*beta)*va]
dx=1.0; dz=2.0;
xl=dx*nx/2; zl=dz*(nz-1);
xx=-xl:dx:xl-dx;
zz=0:dz:zl;
[x,z]=meshgrid(xx,zz);
x=x';z=z';
gamma=1.66667;              %
eta=0.01;
nu=0.05;
beta=0.5;                                  % beta in x=Lx
rho0=1.0;             % rho0 -- mass density rho0 in x=Lx
B0=1.0;                                        % b0 -- B0
bl=3.0;                    % bl -- width of current sheet
va=sqrt(B0*B0/rho0);% Alfven velocity
p0=0.5*beta*B0*B0; % p0 -- pressure p0 in x=Lx
t0=0.5*beta*va*va; % t0 -- temperature T0 in x=Lx
cfl=1.5;  % Courant–Friedrichs–Lewy condition parameter
U=zeros(nx,nz,5);
% % 初始场分布
s=((1:nx)-nx)*dx/bl;
b=B0*tanh(s);
p=p0+0.5*(B0.^2-b.^2);
rho=p/t0;
Ay=B0.*bl.*log(cosh(s));
U(:,:,1)=repmat(rho',1,nz);
U(:,:,2)=repmat(p',1,nz);
U(:,:,5)=repmat(Ay',1,nz);
% % 初始场扰动
nw=1;% nw -- number of initial perturbation
abd=10; % abd -- perturbation width in x
am=[0.01,zeros(1,9)];% am -- prtb amplitude
s=((2:nx)-nx)*dx/abd;
for m=1:nw
    sii=exp(-s.*s).*am(m)*B0;
    sij=sin((2*m*((1:nz)-0.5*nz)/(nz-1)+0.5)*pi)*dz*(nz-1)/(2.0*m*pi);
end
U(:,:,5)=U(:,:,5)+[zeros(1,nz);sii'*sij];
%%
t=0;
rho0=U(:,:,1);p0=U(:,:,2);ux0=U(:,:,3);uz0=U(:,:,4);Ay0=U(:,:,5);[Bx0,Bz0]=calcBxz(t,U);
% figure(1)
% subplot(3,2,1);surf(x,z,rho0); xlabel('x');ylabel('z');title('\rho_0')
% subplot(3,2,2);surf(x,z,p0); xlabel('x');ylabel('z');title('p_0')
% subplot(3,2,3);surf(x,z,ux0); xlabel('x');ylabel('z');title('ux_0')
% subplot(3,2,4);surf(x,z,uz0); xlabel('x');ylabel('z');title('uz_0')
% subplot(3,2,5);surf(x,z,Ay0); xlabel('x');ylabel('z');title('Ay_0')
% subplot(3,2,6);quiver(x,z,ux0,uz0); xlabel('x');ylabel('z');title('u_0')
%% 创建新文件夹
% cd('E:\Vlasov\谢华生\code\chap5_mhd\mhd2d_f90\MHD2D New 3');%新文件夹路径
% dirname=['画图1'];%新的文件夹名
% Order=['mkdir ' dirname];%创建命令
% system(Order);%创建文件夹
%%
for it=1:nt
    t=t+dt;
    [Bx,Bz]=calcBxz(t,U);
    Bxm(it)=max(max(Bx));
    U=RK4(t,U);
    if mod(it,20)==0
        figure(2)
        subplot(3,2,1);surf(x,z,U(:,:,1)); xlabel('x');ylabel('z');title('\rho');axis tight
        subplot(3,2,2);surf(x,z,U(:,:,2)); xlabel('x');ylabel('z');title('p');axis tight
        subplot(3,2,3);surf(x,z,U(:,:,3)); xlabel('x');ylabel('z');title('ux');axis tight
        subplot(3,2,4);surf(x,z,U(:,:,4)); xlabel('x');ylabel('z');title('uz');axis tight
        subplot(3,2,5);surf(x,z,U(:,:,5)); xlabel('x');ylabel('z');title('Ay');axis tight
%         saveas(gcf,['E:\Vlasov\谢华生\code\chap5_mhd\mhd2d_f90\MHD2D New 3\画图1\','五个参量',num2str(1000+it),'.jpg']);
        figure(3)
        subplot(2,1,1);quiver(x,z,Bx,Bz,2); xlabel('x');ylabel('z');title('B');axis tight
        subplot(2,1,2);quiver(x,z,U(:,:,3),U(:,:,4),2); xlabel('x');ylabel('z');title('u');axis tight
%         saveas(gcf,['E:\Vlasov\谢华生\code\chap5_mhd\mhd2d_f90\MHD2D New 3\画图1\','向量场',num2str(1000+it),'.jpg']);        
        figure(4)
        subplot(2,2,1);pcolor(x,z,U(:,:,1)); shading('interp');
        xlabel('x');ylabel('z');title(['\rho(x,z,t), t=',num2str(t)])
        subplot(2,2,2);pcolor(x,z,U(:,:,3)); shading('interp');
        xlabel('x');ylabel('z');title(['u_x(x,z,t), t=',num2str(t)])
        subplot(2,2,3);pcolor(x,z,Bx); shading('interp');
        xlabel('x');ylabel('z');title(['B_x(x,z,t), t=',num2str(t)])
        subplot(2,2,4);
        plot([1:it]*dt,log(Bxm));xlabel('t');ylabel('max(Bx(t))');title(['max(Bx(t)), t=',num2str(t)])
        xlabel('t');ylabel('log(max(Bx))');title(['max(Bx(t)), t=',num2str(t)])
        xlim([0,nt*dt])
%         saveas(gcf,['E:\Vlasov\谢华生\code\chap5_mhd\mhd2d_f90\MHD2D New 3\画图1\','伪彩色图',num2str(1000+it),'.jpg']);
        drawnow
    end
    
end
function Uo=RK4(t,Ui)
 % U(:,:,i),其中 i= 1 2 3 4 5 --> rho p ux uz ay
global dt
k1=rightEq(t,Ui);
U1=Ui+0.5*dt*k1;
k2=rightEq(t,U1);
U2=Ui+0.5*dt*k2;
k3=rightEq(t,U2);
U3=Ui+dt*k3;
k4=rightEq(t,U3);
Uo=Ui+dt*(k1+2*k2+2*k3+k4)/6;
end
function k=rightEq(t,U)
 % U(:,:,i),其中 i= 1 2 3 4 5 --> rho p ux uz ay
global gamma nu eta dx dz Bx Bz
rho=U(:,:,1);
p=U(:,:,2);
ux=U(:,:,3);
uz=U(:,:,4);
Ay=U(:,:,5);
pt=p+0.5*(Bx.^2+Bz.^2);% pt=p+b^2/2
[nx,nz,~]=size(U);
ic=1:nx;
jc=1:nz;
ip=ic+1;ip(nx)=1;
im=ic-1;im(1)=nx;
jp=jc+1;jp(nz)=1;
jm=jc-1;jm(1)=nz;
k(:,:,1)=-ux.*(rho(ip,:)-rho(im,:))/(2*dx) ...
    -uz.*(rho(:,jp)-rho(:,jm))/(2*dz) ...
    -rho.*((ux(ip,:)-ux(im,:))/(2*dx) ...
    +(uz(:,jp)-uz(:,jm))/(2*dz));
k(:,:,2)=-ux.*(p(ip,:)-p(im,:))/(2*dx) ...
    -uz.*(p(:,jp)-p(:,jm))/(2*dz) ...
    -gamma.*p.*((ux(ip,:)-ux(im,:))/(2*dx) ...
    +(uz(:,jp)-uz(:,jm))/(2*dz));
k(:,:,3)=-ux.*(ux(ip,:)-ux(im,:))/(2*dx) ...
    -uz.*(ux(:,jp)-ux(:,jm))/(2*dz) ...
    -1./rho.*(pt(ip,:)-pt(im,:))/(2*dz) ... % 原程序有错误(pt(ipx,:)-pt(imx,:))/(2*dz)中应该是除以dx的
    +1./rho.*(Bx.*(Bx(ip,:)-Bx(im,:))/(2*dx)+Bz.*(Bx(:,jp)-Bx(:,jm))/(2*dz)) ...
    +nu./rho.*((ux(ip,:)+ux(im,:)-2*ux(ic,:))/(dx*dx)+(ux(:,jp)+ux(:,jm)-2*ux(:,jc))/(dz*dz));
k(:,:,4)=-ux.*(uz(ip,:)-uz(im,:))/(2*dx) ...
    -uz.*(uz(:,jp)-uz(:,jm))/(2*dz) ...
    -1./rho.*(pt(:,jp)-pt(:,jm))/(2*dz) ...
    +1./rho.*(Bx.*(Bz(ip,:)-Bz(im,:))/(2*dx)+Bz.*(Bz(:,jp)-Bz(:,jm))/(2*dz)) ...
    +nu./rho.*((uz(ip,:)+uz(im,:)-2*uz(ic,:))/(dx*dx)+(uz(:,jp)+uz(:,jm)-2*uz(:,jc))/(dz*dz));
k(:,:,5)=-ux.*(Ay(ip,:)-Ay(im,:))/(2*dx) ...
    -uz.*(Ay(:,jp)-Ay(:,jm))/(2*dz) ...
    +eta.*((Ay(ip,:)+Ay(im,:)-2*Ay(ic,:))/(dx*dx)+(Ay(:,jp)+Ay(:,jm)-2*Ay(:,jc))/(dz*dz));

%% x方向自由边界条件 % z方向周期性边界条件
%  % U(:,:,i),其中 i= 1 2 3 4 5 --> rho p ux uz ay
% % x direction boundary condition, ux(x=xd)=bz=0, p_x=uz_x=rho_x=0
% % x的下边界,i=1
% k(1,:,1)=k(2,:,1);
% k(1,:,2)=k(2,:,2);
% k(1,:,3)=k(2,:,3);% ->0
% k(1,:,4)=k(2,:,4);
% k(1,:,5)=k(2,:,5);% ->0
%% x的上边界, i=end
%% ux(end,:)=0;Bz(end,:)=0;
% k(end,:,1)=-ux(end,:).*(rho(end,:)-rho(end-1,:))/(1*dx) ...
%     -uz(end,:).*(rho(end,jp)-rho(end,jm))/(2*dz) ...
%     -rho(end,:).*((ux(end,:)-ux(end-1,:))/(1*dx) ...
%     +(uz(end,jp)-uz(end,jm))/(2*dz));
% k(end,:,2)=-ux(end,:).*(p(end,:)-p(end-1,:))/(1*dx) ...
%     -uz(end,:).*(p(end,jp)-p(end,jm))/(1*dz) ...
%     -gamma.*p(end,:).*((ux(end,:)-ux(end-1,:))/(1*dx) ...
%     +(uz(end,jp)-uz(end,jm))/(2*dz));
% k(end,:,3)=-ux(end,:).*(ux(end,:)-ux(end-1,:))/(1*dx) ...
%     -uz(end,:).*(ux(end,jp)-ux(end,jm))/(2*dz) ...
%     -1./rho(end,:).*(pt(end,:)-pt(end-1,:))/(1*dx) ...
%     +1./rho(end,:).*(Bx(end,:).*(Bx(end,:)-Bx(end-1,:))/(1*dx)+Bz(end,:).*(Bx(end,jp)-Bx(end,jm))/(2*dz)) ...
%     +nu./rho(end,:).*((2*ux(end,:)-2*ux(end-1,:))/(dx*dx)+(ux(end,jp)+ux(end,jm)-2*ux(end,jc))/(dz*dz));
% k(end,:,4)=-ux(end,:).*(uz(end,:)-uz(end-1,:))/(1*dx) ...
%     -uz(end,:).*(uz(end,jp)-uz(end,jm))/(2*dz) ...
%     -1./rho(end,:).*(pt(end,jp)-pt(end,jm))/(2*dz) ...
%     +1./rho(end,:).*(Bx(end,:).*(Bz(end,:)-Bz(end-1,:))/(1*dx)+Bz(end,:).*(Bz(end,jp)-Bz(end,jm))/(2*dz)) ...
%     +nu./rho(end,:).*((2*uz(end-1,:)-2*uz(end,:))/(dx*dx)+(uz(end,jp)+uz(end,jm)-2*uz(end,jc))/(dz*dz));
% k(end,:,5)=-ux(end,:).*(Ay(end,:)-Ay(end-1,:))/(1*dx) ...
%     -uz(end,:).*(Ay(end,jp)-Ay(end,jm))/(2*dz) ...
%     +eta.*((2*Ay(end-1,:)-2*Ay(end,:))/(dx*dx)+(Ay(end,jp)+Ay(end,jm)-2*Ay(end,jc))/(dz*dz));

%% 以下是参照谢华生书籍fortran程序的边界处理
% x的下边界, i=1
k(1,:,1)=k(2,:,1);
k(1,:,2)=k(2,:,2);
k(1,:,3)=0;
k(1,:,4)=k(2,:,4);
k(1,:,5)=0;
% x的上边界, i=end
k(end,:,1)=0 ...
    -uz(end,:).*(rho(end,jp)-rho(end,jm))/(2*dz) ...
    -rho(end,:).*((ux(end,:)-ux(end-1,:))/(1*dx) ...
    +(uz(end,jp)-uz(end,jm))/(2*dz));
k(end,:,2)=0 ...
    -uz(end,:).*(p(end,jp)-p(end,jm))/(1*dz) ...
    -gamma.*p(end,:).*((ux(end,:)-ux(end-1,:))/(1*dx) ...
    +(uz(end,jp)-uz(end,jm))/(2*dz));
k(end,:,3)=0;
k(end,:,4)=0 ...
    -uz(end,:).*(uz(end,jp)-uz(end,jm))/(2*dz) ...
    -1./rho(end,:).*(pt(end,jp)-pt(end,jm))/(2*dz) ...
    +1./rho(end,:).*(Bx(end,:).*(Bz(end,:)-Bz(end-1,:))/(1*dx)+0) ...
    +nu./rho(end,:).*((2*uz(end-1,:)-2*uz(end,:))/(dx*dx)+(uz(end,jp)+uz(end,jm)-2*uz(end,jc))/(dz*dz));
k(end,:,5)=0 ...
    -uz(end,:).*(Ay(end,jp)-Ay(end,jm))/(2*dz) ...
    +eta.*((2*Ay(end-1,:)-2*Ay(end,:))/(dx*dx)+(Ay(end,jp)+Ay(end,jm)-2*Ay(end,jc))/(dz*dz));
%%
end

function [Bx,Bz]=calcBxz(t,U)
%********************************************************************************
global dx dz
[nx,nz,~]=size(U);
ic=1:nx;
jc=1:nz;
ip=ic+1;ip(nx)=1;
im=ic-1;im(1)=nx;
jp=jc+1;jp(nz)=1;
jm=jc-1;jm(1)=nz;
Bx=-(U(:,jp,5)-U(:,jm,5))/(2*dz);      
Bz= (U(ip,:,5)-U(im,:,5))/(2*dx);
%% x方向自由边界条件 % z方向周期性边界条件
% x的下边界 i=1
Bx(1,:)=-(U(1,jp,5)-U(1,jm,5))/(2*dz);
Bz(1,:)=(U(2,:,5)-U(1,:,5))/(1*dx);

%% x的上边界 i=end
% Bx(end,:)=-(U(end,jp,5)-U(end,jm,5))/(2*dz);
% Bz(end,:)=(U(end,:,5)-U(end-1,:,5))/(1*dx);
% 以下是参照谢华生书籍fortran程序的边界处理
Bx(end,:)=-(U(end,jp,5)-U(end,jm,5))/(2*dz);
Bz(end,:)=0;
end