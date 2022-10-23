clc;
clear all;
warning off;
cd('C:\Users\IOT Lab\Desktop');%change directory
A=readtable('YBUS.xlsx');%open excel file as table
A=table2array(A);%convert table into array
for i=1:length(A) %Generate impedance matrix
    z(A(i,1),A(i,2))=A(i,3)+A(i,4)*j;
    z(A(i,2),A(i,1))=A(i,3)+A(i,4)*j;
end
% disp (z)
for i=1:length(z) %set the unassigned impedances to infinite
    for j=1:1:length(z)
        if z(i,j)==0
            z(i,j)=inf;
        end
    end
end
y=1./z; %Generate admittance matrix
for i=1:length(z) %generate Y-BUS Matrix
    for j=1:1:length(z)
        if i==j
            Y(i,j)=sum(y(i,:));
        else
            Y(i,j)=-y(i,j);
        end
    end
end
fprintf('Y BUS Matrix is\n')  %print 
disp(Y) %display Y bus 

