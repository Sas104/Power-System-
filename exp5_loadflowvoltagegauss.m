%this code is example of load flow analysis
%It determines bus voltage with gauss seidel method
%Power system I sessional, EEE 3112, RUET
%Written by Sabbir Ahmed Sumon, RUET EEE-18
%On 07 November 2022
%Course Teacher Md. Rashidul Islam, 
%Assistant Professor,EEE, RUET, Bangladesh
clc;
clear all;
warning off;
cd('C:\Users\Asus\OneDrive\Desktop\reports_lab\PowerSystemISessional');% change current directory
A=readtable('gauss.xlsx');%open excel file as table
A=table2array(A);%convert table into array
sz=size(A);
for i=1:sz(1) %Generate impedance matrix
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
% fprintf('Y BUS Matrix is\n')  %print 
% disp(Y);
B=readtable('Gauss2.xlsx');%open excel file as table
B=table2array(B);%convert table into array
a=size(B); %compute size of given matrix
for i=1:a(1) 
    V(B(i,1))=B(i,2);
    Pg(B(i,1))=B(i,3);
    Qg(B(i,1))=B(i,4);
    Pl(B(i,1))=B(i,5);
    Ql(B(i,1))=B(i,6);
end
%e=input('Enter accuracy\n'); 
e=0.00000001;
r=10;
for i=1:a(1)
    if Pg(i)>0
        Vconstant(i)=V(i);
        Pgconstant(i)=Pg(i);
    end
end
for j=1:r
    for i=2:a(1)
        Vprev(i)=V(i);
        if Pg(i)>0
            I=0;
            for k=1:a(1)
                I=I+(Y(i,k)*V(k));
            end
            I(i)=I;
            S(i)=V(i)*conj(I(i));
            S(i)=S(i)+Pg(i)-Pl(i)-real(S(i));
            VY=0;
            for k=1:a(1)
                VY=VY+(Y(i,k)*V(k));
            end 
            VY=VY-V(i)*Y(i,i);
            V(i)=(1/Y(i,i))*(conj(S(i))/conj(V(i))-VY);
            Vangle=angle(V(i));
            Vr=Vconstant(i);
            V(i)=complex(Vr*cos(Vangle),Vr*sin(Vangle));
        else
            S(i)=complex(Pg(i)-Pl(i),(Qg(i)-Ql(i)));
            VY=0;
            for k=1:a(1)
                VY=VY+(Y(i,k)*V(k));
            end 
            VY=VY-V(i)*Y(i,i);
            V(i)=(1/Y(i,i))*(conj(S(i))/conj(V(i))-VY);
        end
    end
%     fprintf('%u iteration\n',j+1)
%     disp(V)
    for j=2:a(1)
        terminate=0;
        if (Vprev(i)-V(i))<e
            terminate=1*terminate
        else
            terminate=0;
        end
    end
    if terminate==1
        break;
    end
end
Vpolar=abs(V);
Vangle=angle(V)*180/pi;
for i=1:a(1)
    fprintf('V%u=%f<%f\n',i,Vpolar(i),Vangle(i));
end
