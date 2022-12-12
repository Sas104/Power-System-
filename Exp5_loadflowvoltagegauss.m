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
cd('C:\Users\Asus\OneDrive\Desktop\reports_lab\PowerSystemISessional');
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
%e=input('Enter tolerance\n'); 
e=0.00000001; %tolerance
r=100;   %maximum number of iteration
for i=1:a(1)
    if Pg(i)>0
        Vconstant(i)=V(i);
    end
end
for j=1:r
    for i=2:a(1)
        Vprev(i)=V(i); %preserving the previous results
        if Pg(i)>0  %for generation bus
            I=0;
            for k=1:a(1) %calculating I with generalized equation
                I=I+(Y(i,k)*V(k));
            end
            I(i)=I;
            S(i)=V(i)*conj(I(i)); %calculating SI
            S(i)=S(i)+Pg(i)-Pl(i)-real(S(i)); %keeping real part of SI unchanged
        else %for load bus
            S(i)=complex(Pg(i)-Pl(i),(Qg(i)-Ql(i)));
        end
            VY=0;
            for k=1:a(1)        %calculating sums of V*Y
                VY=VY+(Y(i,k)*V(k));
            end 
            VY=VY-V(i)*Y(i,i); %subtracting current V*Y from total
            V(i)=(1/Y(i,i))*(conj(S(i))/conj(V(i))-VY); %calculating bus voltage
            if Pg(i)>0  %to keep magnitude of geeration bus voltage constant 
                Vangle=angle(V(i));
                Vr=Vconstant(i);
                V(i)=complex(Vr*cos(Vangle),Vr*sin(Vangle));
            end
    end
%     fprintf('%u iteration\n',j+1); %uncomment to see all iterations.
%     disp(V);
    terminate=1;
    for j=2:a(1) %tolerance termination condition
        if (Vprev(i)-V(i))<e
            terminate=1*terminate;
        else
            terminate=0;
        end
    end
    if terminate==1
        break;
    end
end
Vpolar=abs(V); %calculating manitude
Vangle=angle(V)*180/pi; %calculating Polar angle
for i=1:a(1) %printing the results in polar form 
    fprintf('V%u=%f<%f\n',i,Vpolar(i),Vangle(i));
end
