clear all;clc;
%
%
%   Frans Ebersohn:  Data Analysis
%
%%
nfiles = 100;
dt = 10;
% 
%
format = '%f \t %f \t %f \t %f \t %f \t %f';
%%
for i=0:nfiles
    filename = ['elec1Particles' num2str(i) '.dat'];
    
    A = importdata(filename,'\t',3);   
    
%     x(i+1,:) = A.data(:,1);
%     vx(i+1,:) = A.data(:,3);
%     vy(i+1,:) = A.data(:,4);
%     vz(i+1,:) = A.data(:,5);
%     en(i+1,:) = A.data(:,6);   
    
    t(i+1) = i*dt;
    crossvt = sqrt((pi/(t(i+1)))^2*1/4);    
    h = plot(A.data(:,1), A.data(:,5),'.',0.0, crossvt,'ro');
    axis([-0.50 0.50 0.0 0.05]);
    set(h,{'MarkerSize'},{1;10});
    xlabel('$$v_{z}$$','interpreter','latex');ylabel('$$v_{\perp}$$','interpreter','latex');
    M(i+1) = getframe;
% 
end
%
%%
movie(M);
%
%


