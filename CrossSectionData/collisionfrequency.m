clear all;clc;
%
%
%   Frans Ebersohn: Plot Collision Frequency
%
%%
%m=6.626e-26;
m = 9.1e-31;
n = 1e18;
% 
%
filename = 'collisions_eAr.txt'; 
A = importdata(filename,'\t',0);   
filename = 'electron_ARGON_crosssection_data.txt';
B = importdata(filename,'\t',1);
nu1 = sqrt(B.data(:,1)*2/m)*n.*B.data(:,2);
nu2 = sqrt(B.data(:,1)*2/m)*n.*B.data(:,3);
nu3 = sqrt(B.data(:,1)*2/m)*n.*B.data(:,4);
%

loglog(A(:,1),A(:,2),'o',B.data(:,1)/(1.6e-19),nu1,A(:,1),A(:,3),'o',B.data(:,1)/(1.6e-19),nu2,A(:,1),A(:,4),'o',B.data(:,1)/(1.6e-19),nu3);
xlabel('Energy (eV)');ylabel('Collision Frequency (1/s)'); title('Electron-Argon Collision Frequency');
legend('Elastic: Simulation','Elastic: Analytical', 'Inelastic: Simulation', 'Inelastic: Analytical','Ionization: Simulation', 'Ionization: Analytical');
% loglog(A(:,1),A(:,2),'o',B.data(:,1)/(1.6e-19),nu1,A(:,1),A(:,3),'o',B.data(:,1)/(1.6e-19),nu2);
% xlabel('Energy (eV)');ylabel('Collision Frequency (1/s)'); title('Argon-Argon Collision Frequency');
%legend('Elastic: Simulation','Elastic: Analytical', 'CEX: Simulation', 'CEX: Analytical');