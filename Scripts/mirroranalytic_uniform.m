clear all
close all


x=0:0.001:0.5-0.001;
bz = 1+4.*x.^2;
bmax = 2;
bmin = 1;

rmax = bmax/bmin;
thetam = asin(1/sqrt(rmax));

for k=1:length(x)
    rrr = bz(k)/bmin;
    theta1 = asin(1/sqrt(rrr))-1e-10;
    
    nnn0  = (rrr-2)/2/rrr* ...
        1/cos(thetam)/sqrt(1-(rrr-1)*tan(thetam)*tan(thetam)) ...
        * sqrt((-rrr*cos(2*thetam)+(rrr-2))/(rrr-2)) ...
        *(sqrt((-rrr*cos(2*thetam)+(rrr-2))/(rrr-2)) -1);
    nnn1  = (rrr-2)/2/rrr* ...
        1/cos(theta1)/sqrt(1-(rrr-1)*tan(theta1)*tan(theta1)) ...
        * sqrt((-rrr*cos(2*theta1)+(rrr-2))/(rrr-2)) ...
        *(sqrt((-rrr*cos(2*theta1)+(rrr-2))/(rrr-2)) -1);
    
    theta1save(k) = theta1;
    nnn(k) = nnn1-nnn0;
end

figure(1)
plot(x,nnn)

nnn_non = -cos(theta1save) + cos(thetam);
hold on
plot(x,nnn_non,'r--')
hold off

fp=fopen('analytic_uniform.dat','w');
for k=1:length(x)
    fprintf(fp,'%15.5f %15.5e \n', x(k),nnn(k));
end
fclose(fp);

theta = 0:pi/400:pi/2;

y = sin(theta)./sqrt(1-(1.999-1).*tan(theta).^2);

y1 = sin(theta)./sqrt(1-(1.1-1).*tan(theta).^2);
y2 = sin(theta)./sqrt(1-(1.001-1)*tan(theta).^2);

thetapi = theta/(pi/2);


t1 = asin(1/sqrt(1.999))/(pi/2);
t2 = asin(1/sqrt(1.1))/(pi/2);
t3 = asin(1/sqrt(1.001))/(pi/2);


figure(2)
plot(thetapi,y,'k-',thetapi,y1,'r:',thetapi,y2,'b--')
xlabel('\theta, * \pi/2')
ylabel('Arbitrary')
legend('r=1.999 (marginally trapped)','r=1.1','r=1.001 (deeply trapped)')
hold on
tt1 = [t1 t1];
tt2 = [t2 t2];
tt3 = [t3 t3];
yy1 = [0 40];
yy2 = [0 40];
yy3 = [0 40];

plot(tt1,yy1,'k--', tt2,yy2, 'r--',tt3,yy3,'b--')
hold off
