clear all;clc;
%
%
%   Frans Ebersohn:  VDF data
%
%%
nfiles = 80;
dt = 2;
% 
%
%%
for i=0:nfiles
    filename = ['ARGONtotalvdf' num2str(i) '.dat'];
    
    A = importdata(filename,'\t',2);   
    
     v(:,i+1) = A.data(:,1);
     f(:,i+1) = A.data(:,2);
     vz(:,i+1) = A.data(:,3);
     fz(:,i+1) = A.data(:,4);
     vr(:,i+1) = A.data(:,5);
     fr(:,i+1) = A.data(:,6); 
     vp(:,i+1) = A.data(:,7);
     fp(:,i+1) = A.data(:,8);  

end
%
f_analytic = 4*pi*v(:,10).*v(:,10)*(sqrt(1/(208*10*2*pi)))^3.*exp(-0.5.*v(:,10).^2/(208*10));
fz_analytic = (sqrt(1/(208*10*2*pi))).*exp(-0.5.*vz(:,10).^2/(208*10));
fr_analytic = (sqrt(1/(208*10*2*pi))).*exp(-0.5.*vr(:,10).^2/(208*10));;
fp_analytic = 2*pi.*vp(:,10)*(sqrt(1/(208*10*2*pi)))^2.*exp(-0.5.*vp(:,10).^2/(208*10));
%
%%
%
subplot(2,4,1);
plot(vz(:,1),fz(:,1),vp(:,1),fp(:,1));
axis([-1000.0 1000.0 0.0 0.05]);
title('Initial'); 
xlabel('v');ylabel('f');
legend('fz','fp');
subplot(2,4,2);
plot(vz(:,2),fz(:,2),vp(:,2),fp(:,2),vz(:,10),fz_analytic,vp(:,10),fp_analytic);
axis([-1000.0 1000.0 0.0 0.015]);
title('After 1 Collision');
xlabel('v');ylabel('f');
legend('fz','fp','fz (neutrals)','fp (neutrals)');
subplot(2,4,3);
plot(vz(:,11),fz(:,6),vp(:,6),fp(:,6)',vz(:,10),fz_analytic,vp(:,10),fp_analytic);
axis([-1000.0 1000.0 0.0 0.015]);
title('After 5 Collisions');
xlabel('v');ylabel('f');
legend('fz','fp','fz (neutrals)','fp (neutrals)');
subplot(2,4,4);
plot(vz(:,51),fz(:,51),vp(:,51),fp(:,51),vz(:,10),fz_analytic,vp(:,10),fp_analytic);
axis([-1000.0 1000.0 0.0 0.015]);
title('After 50 Collisions');
xlabel('v');ylabel('f');
legend('fz','fp','fz (neutrals)','fp (neutrals)');

subplot(2,4,5);
plot(v(:,1), f(:,1));
axis([0.0 1000.0 0.0 0.05]);
title('Initial'); 
xlabel('v');ylabel('f');
legend('f');
subplot(2,4,6);
plot(v(:,2), f(:,2),v(:,10),f_analytic,'r');
axis([0.0 1000.0 0.0 0.02]);
title('After 1 Collision');
xlabel('v');ylabel('f');
legend('f','f (neutrals)');
subplot(2,4,7);
plot(v(:,6), f(:,6),v(:,10),f_analytic,'r');
axis([0.0 1000.0 0.0 0.02]);
title('After 5 Collisions');
xlabel('v');ylabel('f');
legend('f','f (neutrals)');
subplot(2,4,8);
plot(v(:,51), f(:,51),v(:,10),f_analytic,'r');
axis([0.0 1000.0 0.0 0.02]);
title('After 50 Collisions');
xlabel('v');ylabel('f');
legend('f','f (neutrals)');

