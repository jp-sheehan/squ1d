clear all;clc;
%
%
%   Frans Ebersohn:  Simulation data
%
%%
% start = 0;
% nfiles = 4;
dt = 1.0;
% 
%
%%
filename1 = ['ARGONOutput_cField' num2str(start) '.dat'];
    
A = importdata(filename1,'\t',3); 

asize = size(A.data(:,3));
numelem = asize(1);

N1 = zeros(1,numelem);
U1 = zeros(1,numelem);
En1 = zeros(1,numelem);
Phi1 = zeros(1,numelem);
Ex1 = zeros(1,numelem-1);
Ey1 = zeros(1,numelem-1);
Ez1 = zeros(1,numelem-1);

for i=(start):(start+nfiles)
    filename1 = ['ARGONOutput_cField' num2str(i) '.dat'];
    
    A = importdata(filename1,'\t',3); 

    N1(:) = N1(:) + A.data(:,3);
    U1(:) = U1(:) + A.data(:,4);
    En1(:) = En1(:) + A.data(:,7);
    Phi1(:) = Phi1(:) + A.data(:,8);  
    xplot = A.data(:,1);

    filename1 = ['ARGONOutput_pField' num2str(i) '.dat'];
     
    A = importdata(filename1,'\t',3);
     
    Ex1(:) = Ex1(:) + A.data(:,3);
    Ey1(:) = Ey1(:) + A.data(:,4);
    Ez1(:) = Ez1(:) + A.data(:,5);
end
%
N1(:) = N1(:)/nfiles;
U1(:) = U1(:)/nfiles;
En1(:) =  En1(:)/nfiles;
Phi1(:) = Phi1(:)/nfiles;   
%     
Ex1(:) = Ex1(:)/nfiles;
Ey1(:) = Ey1(:)/nfiles;
Ez1(:) = Ez1(:)/nfiles;
%
%
%%
cout(:,1) = xplot(:);
cout(:,2) = N1(:);
cout(:,3) = U1(:);
cout(:,4) = En1(:);
cout(:,5) = Phi1(:);
%%
figure;
plot(xplot,Phi1);
%
save('c_average.dat','cout','-ascii');

