clear all;clc;
%
%
%   Frans Ebersohn:  VDF data
%
%%
nfiles = 7;
dt = 2;
% 
%
%%
for i=0:nfiles
    filename = ['electrontotalvdf' num2str(i) '.dat'];
    
    A = importdata(filename,'\t',2);   
    
     v(:,i+1) = A.data(:,1);
     f(:,i+1) = A.data(:,2);
     vx(:,i+1) = A.data(:,3);
     fx(:,i+1) = A.data(:,4);
     vy(:,i+1) = A.data(:,5);
     fy(:,i+1) = A.data(:,6); 
     vz(:,i+1) = A.data(:,7);
     fz(:,i+1) = A.data(:,8);  

end
%
fx_analytic = 1./(log(3)*(2-vx(:,1)));
fy_analytic = 1./(log(3)*sqrt(3+vx(:,1).*vx(:,1)));
%
%%
%
figure;
subplot(1,3,1);
plot(vx(:,1),fx(:,1),vz(:,1),fz(:,1),'g');
axis([-1.0 1.0 0.0 30.0]);
title('Initial'); 
xlabel('v');ylabel('f');
legend('fx','fz');
subplot(1,3,2);
plot(vx(:,2),fx(:,2),vy(:,2),fz(:,2), 'g');
axis([-1.0 1.0 0.0 2.0]);
title('After 1 Collision');
xlabel('v');ylabel('f');
legend('fx','fz','fx analytical', 'fz analytical');
subplot(1,3,3);
plot(vx(:,8),fx(:,8),vz(:,8),fz(:,8),'g');
axis([-1.0 1.0 0.0 2.0]);
title('After 6 Collisions');
xlabel('v');ylabel('f');
legend('fx','fz');