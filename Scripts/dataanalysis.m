clear all;clc;
%
%
%   Frans Ebersohn:  Simulation data
%
%%
nfiles = 100;
dt = 1.0;
% 
%
%%
for i=0:nfiles
    filename1 = ['ARGONOutput_cField' num2str(i) '.dat'];
    
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
     
     filename1 = ['ARGONOutput_pField' num2str(i) '.dat'];
     
     A = importdata(filename1,'\t',3);
     
     xp1(:,i+1) = A.data(:,1);
     Ex1(:,i+1) = A.data(:,3);
     Ey1(:,i+1) = A.data(:,4);
     Ez1(:,i+1) = A.data(:,5);

end
%
% fftN1D = fft(N1(:,1)-mean(N1(:,1)));
% fftN2D = fft2(N1-mean(N1(:,1)));
% fftE1D = fft(Ex1(:,1));
% fftE2D = fft2(Ex1);
%
xplot = x1(:,1);
xpplot = xp1(:,1);
Exint = interp1(xpplot,Ex1(:,101),xplot);
%
% for i=0:(length(xplot)-1)
%     k(i+1) = 2*pi/(2*pi)*(i);
% end
% %
% for i=0:(length(xplot)-2)
%     ke(i+1) = 2*pi/(2*pi)*(i);
% end
%
%
% matsize = size(px1)
% for i=0:(length(t)-1)
%     deltav(i+1) =(sum(pv1(:,i+1).*pv1(:,i+1))/matsize(1) - (sum(pv1(:,i+1))/matsize(1)).^2);
% end
%
%%
% figure; 
% semilogy(t,(Nmax-mean(N1(:,1)))/mean(N1(:,1)),t,0.001*exp(.265/sqrt(5).*t));
% xlabel('t','interpreter','latex'); ylabel('$$n_{max}$$','interpreter','latex'); title('Maximum density in time','interpreter','latex');
% axis([0 100 1e-4 1.0]);
% figure; 
% plot(t,deltav);
% xlabel('t','interpreter','latex'); ylabel('$$\Delta v$$','interpreter','latex'); title('$$\Delta v$$ in time','interpreter','latex');
figure;
subplot(2,1,1);
imagesc(t,xplot,N1-mean(N1(:,1))); colorbar;
%axis([0 62.8 0 0.314159]);
xlabel('t','interpreter','latex'); ylabel('x','interpreter','latex'); title('$$\delta n$$','interpreter','latex');
subplot(2,1,2);
imagesc(t,xplot,-Ex1); colorbar;
%axis([0 62.8 0 0.314159]);
xlabel('t','interpreter','latex'); ylabel('x','interpreter','latex'); title('$$qE_x/m$$','interpreter','latex');
% %
% figure; 
% subplot(2,3,1);
% plot(px1(:,1),pv1(:,1),'.',px2(:,1),pv2(:,1),'.');
% axis([0 8*3.1415*sqrt(3) -3 6]);
% xlabel('x','interpreter','latex'); ylabel('v','interpreter','latex'); title('$$t=0$$','interpreter','latex');
% subplot(2,3,2);
% plot(px1(:,25),pv1(:,25),'.',px2(:,25),pv2(:,25),'.');
% axis([0 8*3.1415*sqrt(3) -3 6]);
% xlabel('x','interpreter','latex'); ylabel('v','interpreter','latex'); title('$$t=4\pi$$','interpreter','latex');
% subplot(2,3,3);
% plot(px1(:,50),pv1(:,50),'.',px2(:,50),pv2(:,50),'.');
% axis([0 8*3.1415*sqrt(3) -3 6]);
% xlabel('x','interpreter','latex'); ylabel('v','interpreter','latex'); title('$$t=8\pi$$','interpreter','latex');
% subplot(2,3,4);
% plot(px1(:,75),pv1(:,75),'.',px2(:,75),pv2(:,75),'.');
% axis([0 8*3.1415*sqrt(3) -3 6]);
% xlabel('x','interpreter','latex'); ylabel('v','interpreter','latex'); title('$$t=12\pi$$','interpreter','latex');
% subplot(2,3,5);
% plot(px1(:,100),pv1(:,100),'.',px2(:,100),pv2(:,100),'.');
% axis([0 8*3.1415*sqrt(3) -3 6]);
% xlabel('x','interpreter','latex'); ylabel('v','interpreter','latex'); title('$$t=16\pi$$','interpreter','latex');
% subplot(2,3,6);
% plot(px1(:,125),pv1(:,125),'.',px2(:,125),pv2(:,125),'.');
% axis([0 8*3.1415*sqrt(3) -3 6]);
% xlabel('x','interpreter','latex'); ylabel('v','interpreter','latex'); title('$$t=20\pi$$','interpreter','latex');
% figure;
% [hC hC] = contourf(omega,k,abs(fftN2D));
% set(hC,'LineStyle','none');
% colorbar;
% axis([0 3 0 10]);% 
% xlabel('$$\omega$$','interpreter','latex');ylabel('$$k_x$$','interpreter','latex');title('Phase Space');
% %
% figure; 
% plot(pv1(:,1),pv2(:,1),'.');
% xlabel('v_x');ylabel('v_y');title('Initial');
% figure; 
% plot(pv1(:,100),pv2(:,100),'.');
% xlabel('v_x');ylabel('v_y');title('Final');

% figure; 
% plot(t,pv1(1,:),'.-',t,pv2(1,:),'.-');
% xlabel('t(s)'); ylabel('v');legend('v_x','v_y');
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


