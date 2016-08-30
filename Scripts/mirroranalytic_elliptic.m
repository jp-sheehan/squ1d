clear all
close all


x=0:0.001:0.5-0.001;
bz = 1+4.*x.^2;
bmax = 2;
bmin = 1;

rmax = bmax/bmin;
thetam = asin(1/sqrt(rmax));  %% maximum trapping 

for k=1:length(x)
    rrr = bz(k)/bmin;   %% r(z)
    theta1 = asin(1/sqrt(rrr))-1e-10;  %% theta_0 (z)
    nnn0  = -1/cos(thetam)*sqrt(2*rrr*cos(2*thetam)-2*rrr+4) ...
              /2/rrr/sqrt(1-(rrr-1)*tan(thetam)*tan(thetam)) ...
              * (ellipticE(thetam,rrr) - ellipticF(thetam,rrr));
    nnn1  = -1/cos(theta1)*sqrt(2*rrr*cos(2*theta1)-2*rrr+4) ...
              /2/rrr/sqrt(1-(rrr-1)*tan(theta1)*tan(theta1)) ...
              * (ellipticE(theta1,rrr) - ellipticF(theta1,rrr));
  
     theta1save(k) = theta1;
    nnn(k) = nnn1-nnn0;
end

figure(1)
plot(x,nnn)

nnn_non = -cos(theta1save) + cos(thetam);
hold on
plot(x,nnn_non,'r--')
hold off


figure(2)
plot(x,nnn_non,'r--')


fp=fopen('analytic_center.dat','w');
for k=1:length(x)
    fprintf(fp,'%15.5f %15.5e \n', x(k),nnn(k));
end
fclose(fp);