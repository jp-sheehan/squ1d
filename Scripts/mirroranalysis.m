clear all;clc;
%
%
%   Frans Ebersohn:  Mirror Data Analysis
%
%%
nfiles = 10;
dt = .314;
%
format = '%f \t %f \t %f \t %f';
%%
for i=0:nfiles
    filename = ['elec1Particles' num2str(i) '.dat'];
    fid = fopen(filename,'r');
    A = importdata(filename,'\t',3);
    fclose(fid);    
    
    x(:,i+1) = A.data(:,1);
    y(:,i+1) = A.data(:,2);
    v1(:,i+1) = A.data(:,3);
    v2(:,i+1) = A.data(:,4);   
    v3(:,i+1) = A.data(:,5);
    en(:,i+1) = A.data(:,6);
    
end