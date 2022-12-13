clc;
clear all;
warning off;
cd('C:\Users\Asus\OneDrive\Desktop\reports_lab\PowerSystemISessional\Exp6');
% A=readtable('Example2Z.xlsx');%open excel file as table
A=readtable('gauss.xlsx');%open excel file as table
A=table2array(A);%convert table into array
sz=size(A);
for i=1:sz(1) %Generate impedance matrix
    z(A(i,1),A(i,2))=A(i,3)+A(i,4)*j;
    z(A(i,2),A(i,1))=A(i,3)+A(i,4)*j;
end
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
% B=readtable('Example2Bus.xlsx');%open excel file as table
B=readtable('Gauss2.xlsx');%open excel file as table
B=table2array(B);%convert table into array
a=size(B); %compute size of given matrix
n=a(1);
for i=1:n 
    V(B(i,1))=B(i,2);
    Pg(B(i,1))=B(i,3);
    Qg(B(i,1))=B(i,4);
    Pl(B(i,1))=B(i,5);
    Ql(B(i,1))=B(i,6);
end
%e=input('Enter tolerance\n'); 
e=0.000001; %tolerance
r=10;   %maximum number of iteration
Pg(isnan(Pg))=0;
Qg(isnan(Qg))=0;
Pl(isnan(Pl))=0;
Ql(isnan(Ql))=0;
for i=1:n
    if Pg(i)>0 %generation bus voltage is constant
        Vconstant(i)=V(i);
    else 
        Vconstant(i)=0;
    end
end
lamda=angle(V);
for j=1:r
    Vprev=V;
    lamdaprev=lamda;
    fprintf('Iteration %u\n',j);
    for i=2:n 
        current=0;
        for k=1:n %calculating I with generalized equation
            current=current+(Y(i,k)*V(k));
        end
        I(i)=current;
        S(i)=V(i)*conj(I(i)); %calculating S(i)
        f(i)=real(S(i)); %taking real part of S as f
        g(i)=imag(S(i)); %taking imaginary part of S as g
        delP(i)=Pg(i)-Pl(i)-f(i); %calculating delP
        delQ(i)=Qg(i)-Ql(i)-g(i); %calculating delQ
    end
    delP(:,1)=[];  %eliminating 1st row of delP
    for i=1:n    
        if Vconstant(i)~=0 %generation bus Q 
            delQ(:,i)=[]; %eliminates corresponding rows of generation bus
        end
    end
    delQ(:,1)=[];
    PQ=[delP';delQ'];
    for i=2:n
        for k=1:n
            if i~=k %for non diagonal elements
                H(i,k)=(-1)*abs(V(i))*abs(V(k))*abs(Y(i,k))*sin(angle(Y(i,k))-lamda(i)+lamda(k));
                L(i,k)=abs(V(i))*abs(Y(i,k))*cos(angle(Y(i,k))-lamda(i)+lamda(k));
                M(i,k)=(-1)*abs(V(i))*abs(V(k))*abs(Y(i,k))*cos(angle(Y(i,k))-lamda(i)+lamda(k));
                N(i,k)=(-1)*abs(V(i))*abs(Y(i,k))*sin(angle(Y(i,k))-lamda(i)+lamda(k));
            else  %Setting diagonal elements to zero for now
                    H(i,i)=0;
                    L(i,i)=0;
                    M(i,i)=0;
                    N(i,i)=0;
            end
        end
    end
    for i=2:n      
        for k=1:n %determinig the sum of non-diagonal terms of each row
            if i~=k
                H(i,i)=H(i,i)+H(i,k);
                L(i,i)=L(i,i)+L(i,k);
                M(i,i)=M(i,i)+M(i,k);
                N(i,i)=N(i,i)+N(i,k);
            end
        end
        H(i,i)=H(i,i)*(-1); %determining final values of diagonal elements
        L(i,i)=2*abs(V(i))*abs(Y(i,i))*cos(angle(Y(i,i)))+L(i,i);
        M(i,i)=M(i,i)*(-1);
        N(i,i)=-2*abs(V(i))*abs(Y(i,i))*sin(angle(Y(i,i)))+N(i,i);
    end
    for i=1:n    
        if Vconstant(i)~=0 %LMN are edited for generation bus 
            L(:,i)=[];  %eliminating rows of L refering generation bus
            M(i,:)=[];  %eliminating columns of M refering generation bus
            N(i,:)=[];  %eliminating columns of N refering generation bus
            N(:,i)=[];  %eliminating rows of N refering generation bus
        end
    end
    H(1,:)=[];  %eliminating 1st column of H
    H(:,1)=[];  %eliminating 1st rows of H
    L(1,:)=[];  %eliminating 1st column  of L
    L(:,1)=[];  %eliminating 1st rows of L
    M(1,:)=[];  %eliminating 1st column  of M
    M(:,1)=[];  %eliminating 1st row of M
    N(1,:)=[];  %eliminating 1st column  of N
    N(:,1)=[];  %eliminating 1st row of N
    J=[H L;M N];  %forming Jaccobian matrix
    dlamdadv=inv(J)*PQ;  %calucating dellamda-delv
    lamda(1)=0;  %setting slack bus lamda to 1
    g=0; %for tracking number of generation bus
    for i=2:n
        lamda(i)=lamda(i)+dlamdadv(i-1);
        if Vconstant(i)~=0
            V(i)=Vconstant(i);
            g=g+1;   
        else
            V(i)=abs(V(i))+dlamdadv(n+i-g-2);
        end
    end
    for i=1:n
        fprintf("V(%u)=%.6f<%.6f\n",i,V(i),(lamda(i)*180/pi));
        V(i)=complex(abs(V(i))*cos(lamda(i)),abs(V(i))*sin(lamda(i)));
    end
    terminate1=1;
    terminate2=1;
    for i=2:a(1) %tolerance termination condition
        if (Vprev(i)-V(i))<e
            terminate1=1*terminate1;
        else
            terminate1=0;
        end
    end
    for i=2:a(1) %tolerance termination condition
        if (lamdaprev(i)-lamda(i))<e
            terminate2=1*terminate2;
        else
            terminate2=0;
        end
    end
    if terminate1==1 && terminate2==1
        break;
    end
end
fprintf("\nFinal result\n");
for i=1:n
    fprintf("V(%u)=%.6f<%.6f\n",i,abs(V(i)),(lamda(i)*180/pi));
end