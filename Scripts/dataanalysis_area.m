clear all;clc;
%
%
%   Frans Ebersohn:  Simulation data
%
%%
nfiles = 50;
dt = 1.0e-11;
q= -1.602e-19;
m = 9.1e-31;
eps = 8.854e-12;
L = 0.1;
A_in = 1.0;
% 
%
%%
for i=0:nfiles
    filename1 = ['electronOutput_cField' num2str(i) '.dat'];
    
    A = importdata(filename1,'\t',3); 
    
     x1(:,i+1) = A.data(:,1);
     y1(:,i+1) = A.data(:,2);
     N1(:,i+1) = A.data(:,3);
     U1(:,i+1) = A.data(:,4);
     En1(:,i+1) = A.data(:,7);
     Phi1(:,i+1) = A.data(:,8);
     Bx1(:,i+1) = A.data(:,9); 
     By1(:,i+1) = A.data(:,10);
     Bz1(:,i+1) = A.data(:,11);     
     
     t(i+1) = (i+1)*dt; 
     
     filename1 = ['electronOutput_pField' num2str(i) '.dat'];
     
     A = importdata(filename1,'\t',3);
     
     xp1(:,i+1) = A.data(:,1);
     Ex1(:,i+1) = A.data(:,3);
     Ey1(:,i+1) = A.data(:,4);
     Ez1(:,i+1) = A.data(:,5);

end
%
%
%
%%
%
IT = nfiles+1;
xplot = x1(:,1);
xpplot = xp1(:,1);
Area = A_in*Bx1(1,1)./Bx1(:,1);
Exint = interp1(xpplot,Ex1(:,IT),xplot);
%
sizeof = size(N1);
gridpts = sizeof(1);
dx = L/(gridpts-2);
%
%
massflux = N1(:,IT).*U1(:,IT).*Area;
energy = 0.5*U1(:,IT).*U1(:,IT)+q/m*Phi1(:,IT);
%
momlhs = m*N1(gridpts-1,IT)*U1(gridpts-1,IT)^2*Area(gridpts-1)-m*N1(2,IT)*U1(2,IT)^2*Area(2);
momrhs = 0.0;
%
poislhs = Ex1(gridpts-1,IT)*(Area(gridpts-1)+Area(gridpts))/2-Ex1(1,IT)*(Area(2)+Area(1))/2;
poisrhs = 0.0;
%
momrhs = trapz(xplot(2:gridpts-1), q.*N1(2:(gridpts-1),IT).*Exint(2:(gridpts-1)).*Area(2:(gridpts-1)));    
poisrhs = trapz(xplot(2:(gridpts-1)), q/eps.*N1(2:(gridpts-1),IT).*Area(2:(gridpts-1)));
%    
%%
%
% figure;
% plot(xplot,Area,xplot,Bx1(:,IT));
figure;
plot(xplot(2:(gridpts-1)),massflux(2:(gridpts-1)));
axis([0 0.1 0 1e19]); xlabel('x');ylabel('mass flux');
figure;
plot(xplot(2:(gridpts-1)),energy(2:(gridpts-1)));
axis([0 0.1 0 2e15]);xlabel('x');ylabel('energy');
%
% figure;
% subplot(2,1,1);
% imagesc(t,xplot,N1-mean(N1(:,1))); colorbar;
% %axis([0 62.8 0 0.314159]);
% xlabel('t','interpreter','latex'); ylabel('x','interpreter','latex'); title('$$\delta n$$','interpreter','latex');
% subplot(2,1,2);
% imagesc(t,xplot,-Ex1); colorbar;
% %axis([0 62.8 0 0.314159]);
% xlabel('t','interpreter','latex'); ylabel('x','interpreter','latex'); title('$$qE_x/m$$','interpreter','latex');
% 
%%
% for i=0:nfiles
%     filename = ['elec1Particles' num2str(i) '.dat'];
%     
%     A = importdata(filename,'\t',3);   
%     
% %     x(i+1,:) = A.data(:,1);
% %     vx(i+1,:) = A.data(:,3);
% %     vy(i+1,:) = A.data(:,4);
% %     vz(i+1,:) = A.data(:,5);
% %     en(i+1,:) = A.data(:,6);   
%     
%     t(i+1) = i*dt;
%     crossvt = sqrt((pi/(t(i+1)))^2*1/4);    
%     h = plot(A.data(:,1), A.data(:,5),'.',0.0, crossvt,'ro');
%     axis([-0.50 0.50 0.0 0.05]);
%     set(h,{'MarkerSize'},{1;10},{'LineWidth'},{1;3});
%     xlabel('$$z$$','interpreter','latex');ylabel('$$v_{\perp}$$','interpreter','latex');
%     M(i+1) = getframe;
% % 
% end
% %
% %%
% movie(M);
% %
% myVideo = VideoWriter('myfile.avi');
% myVideo.FrameRate = 15;  % Default 30
% myVideo.Quality = 100;    % Default 75
% open(myVideo);
% writeVideo(myVideo, M);
% close(myVideo);
%


