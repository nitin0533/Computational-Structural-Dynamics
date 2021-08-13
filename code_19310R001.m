
clc
clear all
close all
format long
load('SubStructure1_Matlab.mat');
load('SubStructure2_Matlab.mat');
load('SubStructure3_Matlab.mat');
load('SubStructure4_Matlab.mat');     
% creating identity matrices of the order of internal and boundary degrees
% of freedom of the individual substructures.
I11=sparse(eye(1256));
I22=sparse(eye(440));
I33=sparse(eye(840));
I44=sparse(eye(440));
I12=sparse(eye(22));
I13=sparse(eye(42));
I14=sparse(eye(22));
%% Assembling the global mass and stiffness matrix by using the interface continuity conditions
% Cosntructing the tranformation matrix
row12=length(I11)+length(I12)+length(I13)+length(I14)+length(I22)+length(I12)+length(I33)...
      +length(I13)+length(I44)+length(I14);
col12=length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+length(I13)+length(I14);
T=sparse(row12,col12);

T([1:length(I11)],[1:length(I11)])=I11;

tempr=length(I11)+1:length(I11)+length(I12);
tempc=length(I11)+length(I22)+length(I33)+length(I44)+1:...
    length(I11)+length(I22)+length(I33)+length(I44)+length(I12);
T([tempr],[tempc])=I12;

tempr=length(I12)+length(I11)+1:length(I12)+length(I11)+length(I13);
tempc=length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+1:...
    length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+length(I13);
T([tempr],[tempc])=I13;

tempr=length(I12)+length(I11)+length(I13)+1:length(I12)+length(I11)+length(I13)+length(I14);
tempc=length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+length(I13)+1:...
    length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+length(I13)+length(I14);
T([tempr],[tempc])=I14;

tempr=length(I12)+length(I11)+length(I13)+length(I14)+1:...
    length(I12)+length(I11)+length(I13)+length(I14)+length(I22);
tempc=length(I11)+1:length(I11)+length(I22);
T([tempr],[tempc])=I22;

tempr=length(I12)+length(I11)+length(I13)+length(I14)+length(I22)+1:...
    length(I12)+length(I11)+length(I13)+length(I14)+length(I22)+length(I12);
tempc=length(I11)+length(I22)+length(I33)+length(I44)+1:...
    length(I11)+length(I22)+length(I33)+length(I44)+length(I12);
T([tempr],[tempc])=I12;

tempr=length(I12)+length(I11)+length(I13)+length(I14)+length(I22)+length(I12)+1:...
    length(I12)+length(I11)+length(I13)+length(I14)+length(I22)+length(I12)+length(I33);
tempc=length(I11)+length(I22)+1:length(I11)+length(I22)+length(I33);
T([tempr],[tempc])=I33;

tempr=length(I12)+length(I11)+length(I13)+length(I14)+length(I22)+length(I12)+length(I33)+1:...
      length(I12)+length(I11)+length(I13)+length(I14)+length(I22)+length(I12)+length(I33)+length(I13);
tempc=length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+1:...
    length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+length(I13);
T([tempr],[tempc])=I13;

tempr=length(I12)+length(I11)+length(I13)+length(I14)+length(I22)+length(I12)+length(I33)+length(I13)+1:...
      length(I12)+length(I11)+length(I13)+length(I14)+length(I22)+length(I12)+length(I33)+length(I13)+length(I44);
tempc=length(I11)+length(I22)+length(I33)+1:length(I11)+length(I22)+length(I33)+length(I44);
T([tempr],[tempc])=I44;

tempr=length(I12)+length(I11)+length(I13)+length(I14)+length(I22)+length(I12)+length(I33)+length(I13)+length(I44)+1:...
      length(I12)+length(I11)+length(I13)+length(I14)+length(I22)+length(I12)+length(I33)+length(I13)+length(I44)+length(I14);
tempc=length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+length(I13)+1:...
    length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+length(I13)+length(I14);
T([tempr],[tempc])=I14;

M_ass=[M1 zeros(size(M1,1),size(M2,2)) zeros(size(M1,1),size(M3,2)) zeros(size(M1,1),size(M4,2));...
      zeros(size(M2,1),size(M1,2)) M2 zeros(size(M2,1),size(M3,2)) zeros(size(M2,1),size(M4,2));...
      zeros(size(M3,1),size(M1,2)) zeros(size(M3,1),size(M2,2)) M3 zeros(size(M3,1),size(M4,2));...
      zeros(size(M4,1),size(M1,2)) zeros(size(M4,1),size(M2,2)) zeros(size(M4,1),size(M3,2)) M4];

K_ass=[K1 zeros(size(K1,1),size(K2,2)) zeros(size(K1,1),size(K3,2)) zeros(size(K1,1),size(K4,2));...
      zeros(size(K2,1),size(K1,2)) K2 zeros(size(K2,1),size(K3,2)) zeros(size(K2,1),size(K4,2));...
      zeros(size(K3,1),size(K1,2)) zeros(size(K3,1),size(K2,2)) K3 zeros(size(K3,1),size(K4,2));...
      zeros(size(K4,1),size(K1,2)) zeros(size(K4,1),size(K2,2)) zeros(size(K4,1),size(K3,2)) K4];
% 
 M1=T'*M_ass*T;
 K1=T'*K_ass*T;
%% Rearranging the global mass and stiffness matrix to the form mentioned in Pierre's paper.
row12=length(I11)+length(I12)+length(I13)+length(I14)+length(I22)+length(I33)...
      +length(I44);
col12=length(I11)+length(I12)+length(I13)+length(I14)+length(I22)+length(I33)...
     +length(I44);
T=sparse(row12,col12);

tempr=1:length(I11);
tempc=length(I12)+length(I13)+length(I14)+1:length(I12)+length(I13)+length(I14)+length(I11);
T([tempr],[tempc])=I11;

tempr=length(I11)+1:length(I11)+length(I22);
tempc=length(I12)+length(I13)+length(I14)+length(I11)+1:...
   length(I12)+length(I13)+length(I14)+length(I11)+length(I22);
T([tempr],[tempc])=I22;

tempr=length(I11)+length(I22)+1:length(I11)+length(I22)+length(I33);
tempc=length(I12)+length(I13)+length(I14)+length(I11)+length(I22)+1:...
    length(I12)+length(I13)+length(I14)+length(I11)+length(I22)+length(I33);
T([tempr],[tempc])=I33;

tempr=length(I11)+length(I22)+length(I33)+1:...
    length(I11)+length(I22)+length(I33)+length(I44);
tempc=length(I12)+length(I13)+length(I14)+length(I11)+length(I22)+length(I33)+1:...
       length(I12)+length(I13)+length(I14)+length(I11)+length(I22)+length(I33)+length(I44);
T([tempr],[tempc])=I44;

tempr=length(I11)+length(I22)+length(I33)+length(I44)+1:...
     length(I11)+length(I22)+length(I33)+length(I44)+length(I12);
tempc=1:length(I12);
T([tempr],[tempc])=I12;

tempr=length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+1:...
     length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+length(I13);
tempc=length(I12)+1:length(I12)+length(I13);
T([tempr],[tempc])=I13;

tempr=length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+length(I13)+1:...
    length(I11)+length(I22)+length(I33)+length(I44)+length(I12)+length(I13)+length(I14);
tempc=length(I12)+length(I13)+1:length(I12)+length(I13)+length(I14);
T([tempr],[tempc])=I14;

M=T'*M1*T;
K=T'*K1*T;
[V,D]=eigs(K,M,10,'sm');
for i=1:10;
    w(i)=sqrt(D(i,i));                %determining w_{actual}
end
%%
b=20;  %number of fixed interface modes to be considered 

k1=K([86+1:86+1256],[86+1:86+1256]);
m1=M([86+1:86+1256],[86+1:86+1256]);
[v1,d1]=eigs(k1,m1,b,'sm');

k2=K([86+1256+1:86+1256+440],[86+1256+1:86+1256+440]);
m2=M([86+1256+1:86+1256+440],[86+1256+1:86+1256+440]);
[v2,d2]=eigs(k2,m2,b,'sm');

k3=K([86+1256+440+1:86+1256+440+840],[86+1256+440+1:86+1256+440+840]);
m3=M([86+1256+440+1:86+1256+440+840],[86+1256+440+1:86+1256+440+840]);
[v3,d3]=eigs(k3,m3,b,'sm');

k4=K([86+1256+440+840+1:86+1256+440+840+440],[86+1256+440+840+1:86+1256+440+840+440]);
m4=M([86+1256+440+840+1:86+1256+440+840+440],[86+1256+440+840+1:86+1256+440+840+440]);
[v4,d4]=eigs(k4,m4,b,'sm');

k12_1=K([86+1:86+1256],[1:22]);
k13_1=K([86+1:86+1256],[22+1:22+42]);
k14_1=K([86+1:86+1256],[22+42+1:22+42+22]);

k12_2=K([86+1256+1:86+1256+440],[1:22]);
k13_2=K([86+1256+1:86+1256+440],[22+1:22+42]);
k14_2=K([86+1256+1:86+1256+440],[22+42+1:22+42+22]);

k12_3=K([86+1256+440+1:86+1256+440+840],[1:22]);
k13_3=K([86+1256+440+1:86+1256+440+840],[22+1:22+42]);
k14_3=K([86+1256+440+1:86+1256+440+840],[22+42+1:22+42+22]);

k12_4=K([86+1256+440+840+1:86+1256+440+840+440],[1:22]);
k13_4=K([86+1256+440+840+1:86+1256+440+840+440],[22+1:22+42]);
k14_4=K([86+1256+440+840+1:86+1256+440+840+440],[22+42+1:22+42+22]);

si1_12=-k1\k12_1;
si1_13=-k1\k13_1;  
si1_14=-k1\k14_1;

si2_12=-k2\k12_2;
si2_13=-k2\k13_2;
si2_14=-k2\k14_2;

si3_12=-k3\k12_3;
si3_13=-k3\k13_3;
si3_14=-k3\k14_3;

si4_12=-k4\k12_4;
si4_13=-k4\k13_4;
si4_14=-k4\k14_4;


Ib1=eye(22);
Ib2=eye(42);
Ib3=eye(22);

I=[Ib1 zeros(size(Ib1,1),size(Ib2,2)) zeros(size(Ib1,1),size(Ib3,2));...
   zeros(size(Ib2,1),size(Ib1,2)) Ib2 zeros(size(Ib2,1),size(Ib3,2));...
   zeros(size(Ib3,1),size(Ib1,2)) zeros(size(Ib3,1),size(Ib2,2)) Ib3];

psic= [si1_12 si1_13 si1_14;...
       si2_12 si2_13 si2_14;...
       si3_12 si3_13 si3_14;...
       si4_12 si4_13 si4_14];

T1=[I;psic];

T2=[zeros(size(I,1),size(v1,2));...                  
    v1;...
    zeros(size(si2_12,1),size(v1,2));...
    zeros(size(si3_13,1),size(v1,2));...
    zeros(size(si4_14,1),size(v1,2))];

T3=[zeros(size(I,1),size(v2,2));...
    zeros(size(si1_12,1),size(v2,2));...
    v2;...
    zeros(size(si3_13,1),size(v2,2));...
    zeros(size(si4_14,1),size(v2,2))];

T4=[zeros(size(I,1),size(v3,2));...
    zeros(size(si1_12,1),size(v3,2));...
    zeros(size(si2_12,1),size(v3,2));...
    v3;...
    zeros(size(si4_14,1),size(v3,2))];

T5=[zeros(size(I,1),size(v4,2));
    zeros(size(si1_12,1),size(v4,2));...
    zeros(size(si2_12,1),size(v4,2));...
    zeros(size(si3_13,1),size(v4,2));...
    v4];

Rcb=[T1 T2 T3 T4 T5];

M_cb=Rcb'*M*Rcb;
K_cb=Rcb'*K*Rcb;

[V_cb,D_cb]=eigs(K_cb,M_cb,10,'sm');              

for i=1:10;
    w_cb(i)=sqrt(D_cb(i,i));               %determining the frequencies by Craig-Bampton method
end

%% Using Characeristic constraint mode method
temp=[5:5:60];
error=zeros(length(w_cb),length(temp));
ws=zeros(10,3);
ws(:,1)=real(w');
ws(:,2)=real(w_cb');
mc=M_cb([1:86],[1:86]);
kc=K_cb([1:86],[1:86]);

for j=1:length(temp)
    
[vc,dc]=eigs(kc,mc,temp(j),'sm');

Tcc=[vc zeros(size(vc,1),b) zeros(size(vc,1),b) zeros(size(vc,1),b) zeros(size(vc,1),b);...
    zeros(b,size(vc,2)) eye(b) zeros(b,b) zeros(b,b) zeros(b,b);...
    zeros(b,size(vc,2)) zeros(b,b)  eye(b) zeros(b,b) zeros(b,b);...
    zeros(b,size(vc,2)) zeros(b,b) zeros(b,b) eye(b) zeros(b,b);...
    zeros(b,size(vc,2)) zeros(b,b) zeros(b,b) zeros(b,b) eye(b)];

M_cc=Tcc'*M_cb*Tcc;
K_cc=Tcc'*K_cb*Tcc;
[vcc,dcc]=eigs(K_cc,M_cc,10,'sm');
for i=1:10;
    w_cc(i,j)=sqrt(dcc(i,i));  %determining the frequencies with CC method.
end
end
%% error analysis and plots
for  i=1:12
    error(:,i)=real(w)-real(w_cc(:,i)');  
end
for i=4:10
error(i,:)=error(i,:)./real(w(i))*100;
end

figure(1)
plot(temp,error(9,:),'LineWidth',2,'Marker','*')
grid on
hold on
plot(temp,error(10,:),'LineWidth',2,'Marker','+')
xlabel('Number of Charaterisitc constraint modes')
ylabel('Percentage error (%)')
title('Error Analysis with respect to the exact frequencies (w_{cc} vs w_{actual})')
legend('9th Natural frequency','10th Natural frequency')

figure(2)
plot(temp,error(8,:),'LineWidth',2,'Marker','*')
grid on
 hold on
plot(temp,error(7,:),'LineWidth',2,'Marker','+')
xlabel('Number of Charaterisitc constraint modes')
ylabel('Percentage error (%)')
title('Error Analysis with respect to the exact frequencies (w_{cc} vs w_{actual})')
legend('8th Natural frequency','7th Natural frequency')

figure(3)
plot(temp,error(6,:),'LineWidth',2,'Marker','*')
grid on
xlabel('Num-ber of Charaterisitc constraint modes')
ylabel('Percentage error (%)')
title('Error Analysis with respect to the exact frequencies (w_{cc} vs w_{actual})')
legend('6th Natural frequency')
figure(4)
plot(temp,error(5,:),'LineWidth',2,'Marker','*')
grid on
xlabel('Number of Charaterisitc constraint modes')
ylabel('Percentage error (%)')
title('Error Analysis w.r.t the exact frequencies (w_{cc} vs w_{actual})')
hold on
plot(temp,error(4,:),'LineWidth',2,'Marker','+')
legend('5th Natural frequency','4th Natural frequency')

for  i=1:12
    error(:,i)=real(w_cb)-real(w_cc(:,i)');  
end

for i=4:10
error(i,:)=error(i,:)./real(w_cb(i))*100;
end


figure(5)
plot(temp,error(9,:),'LineWidth',2,'Marker','*')
grid on
hold on
plot(temp,error(10,:),'LineWidth',2,'Marker','+')
xlabel('Number of Charaterisitc constraint modes')
ylabel('Percentage error (%)')
title('Error Analysis with respect to the exact frequencies (w_{cc} vs w_{cb})')
legend('9th Natural frequency','10th Natural frequency')

figure(6)
plot(temp,error(8,:),'LineWidth',2,'Marker','*')
grid on
 hold on
plot(temp,error(7,:),'LineWidth',2,'Marker','+')
xlabel('Number of Charaterisitc constraint modes')
ylabel('Percentage error (%)')
title('Error Analysis with respect to the exact frequencies (w_{cc} vs w_{cb})')
legend('8th Natural frequency','7th Natural frequency')

figure(7)
plot(temp,error(6,:),'LineWidth',2,'Marker','*')
grid on
xlabel('Num-ber of Charaterisitc constraint modes')
ylabel('Percentage error (%)')
title('Error Analysis with respect to the exact frequencies (w_{cc} vs w_{cb})')
legend('6th Natural frequency')

figure(8)
plot(temp,error(5,:),'LineWidth',2,'Marker','*')
grid on
xlabel('Number of Charaterisitc constraint modes')
ylabel('Percentage error (%)')
title('Error Analysis with respect to the exact frequencies (w_{cc} vs w_{cb})')
hold on
plot(temp,error(4,:),'LineWidth',2,'Marker','+')
legend('5th Natural frequency','4th Natural frequency')

ws(:,3)=real(w_cc(:,12));
disp('     Actual                 CB               CC')   
disp(ws)

































































