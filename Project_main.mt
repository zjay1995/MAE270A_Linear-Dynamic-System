
clc; clear all; close all
tic
% ns = 7 for Task 2
for iii = 1:4
%% ns 
ns      = [6,7,10,20];
ns      = ns(iii);

load u1_impulse.mat
y11     = u1_impulse.Y(3).Data;
y21     = u1_impulse.Y(4).Data;
u1      = u1_impulse.Y(1).Data; %%% note that the pulse magnitude is 5
[m,mi]  = max(u1>0); %%% find index where pulse occurs
load u2_impulse.mat
y12     = u2_impulse.Y(3).Data;
y22     = u2_impulse.Y(4).Data;
u2      = u2_impulse.Y(2).Data;

%%% remove any offsets in output data using data prior to pulse application
y11     = y11 - mean(y11([1:mi-1]));
y12     = y12 - mean(y12([1:mi-1]));
y21     = y21 - mean(y21([1:mi-1]));
y22     = y22 - mean(y22([1:mi-1]));
%%% rescale IO data so that impulse input has magnitude 1
y11     = y11/max(u1);
y12     = y12/max(u2);
y21     = y21/max(u1);
y22     = y22/max(u2);
u1      = u1/max(u1);
u2      = u2/max(u2);
ts      = 1/40; %%%% sample period
N       = length(y11); %%%% length of data sets
t       = [0:N-1]*ts - 1;
%% Task1_p1: Construct Hankel Matrix
cl      = 100;
bi      = 41;
H       = zeros(cl*2);
for i = 1:cl
    for j = 1:cl
        k   =i+j-1;
        H(2*i-1,j*2-1)=y11(k+bi); %hk(1,1)
        H(2*i-1,j*2)=y12(k+bi);   %hk(1,2)
        H(2*i,j*2-1)=y21(k+bi); %hk(2,1)
        H(2*i,j*2)=y22(k+bi); %hk(2,2)
    end
end 
x       = 1:cl;
[U,S,V] =svd(H);
[T d]   = eig(H);

s       = diag(S);



[U,S,V] =svd(H); 
Si      = zeros(cl*2);
for i = 1:ns
    Si(i,i) = S(i,i);
end 
O       = U;
C       = Si*V';
Hh      = zeros(cl*2);

for i = 1:cl
    for j = 1:cl
        k   =i+j;
        Hh(2*i-1,j*2-1)=y11(k+bi); %hk(1,1)
        Hh(2*i-1,j*2)=y12(k+bi);   %hk(1,2)
        Hh(2*i,j*2-1)=y21(k+bi); %hk(2,1)
        Hh(2*i,j*2)=y22(k+bi); %hk(2,2)
    end
end 
U       = U(:,1:ns);
V       = V(:,1:ns);
Si      = Si(1:ns,1:ns);

A       = U'*Hh*(V*inv(Si));
B       = C(1:ns,1:2);
C       = O(1:2,1:ns);
D       = zeros (2);

h = zeros(2,2*cl);
    for k = 1:cl
        a           = C*A^(k-1)*B;
        h(1,2*k-1)  = a(1,1);
        h(2,2*k-1)  = a(2,1);
        h(1,2*k)    = a(1,2);
        h(2,2*k)    = a(2,2);

        qa(k+bi)    = a(1,1);
        qb(k+bi)    = a(2,1);
        qc(k+bi)    = a(1,2);
        qd(k+bi)    = a(2,2);
    end
    
%Stability Check
scheck      = max(abs(eig(A)));    
%% Task1_p2: Plot Response   

plot1(u1,y11,y21,y12,y22,u2,qa,qb,qc,qd,t)

%% Task1_p3: 
I   = eye(ns);
%sample period
ts  = 1/40; %s
%Nyquist frequency
is  = 205;

w   = 20; %hz
w   = w*2*pi; %rad/s
w   = linspace(0,w,is);
% model frequency resposne
t   = 0;

for i = 1:is %size(100) 
    fr_s    = (C*inv(exp(1j*w(i)*ts)*I-A)*B); %+D Transfer Function
    a11(i)  = fr_s(1,1);
    a12(i)  = fr_s(1,2);
    a21(i)  = fr_s(2,1);
    a22(i)  = fr_s(2,2);
    frs_s{i}   = svd(fr_s);
end
        mag_a11 = abs(a11);mag_a12 = abs(a12);mag_a21 = abs(a21);mag_a22 = abs(a22);
        ang_a11 = angle(a11);ang_a12 = angle(a12); ang_a21 = angle(a21); ang_a22 = angle(a22);
       
        y11f    = fft(y11)./fft(u1);
        N       = length(y11f);
        om1     = [0:N-1]/(ts*N);
        y21f    = fft(y21)./fft(u1);
        N       = length(y21f);
        om2     = [0:N-1]/(ts*N);
        y12f    = fft(y12)./fft(u2);
        N       = length(y12f);
        om3     = [0:N-1]/(ts*N);
        y22f    = fft(y22)./fft(u2);
        N       = length(y22f);
        om4     = [0:N-1]/(ts*N);
        
         %plot2(mag_a11,mag_a12,mag_a21,mag_a22,y11f,y12f,y21f,y22f,ang_a11,ang_a12,ang_a21,ang_a22)
end
for i = 1:is
    frs_e{i} = svd ([y11f(i)    y12f(i);
                     y21f(i)    y22f(i)]);
end 
%% Task 2_p1:
%confirm rank(S) = 9
for al = 1:100
    S1 = [   al*eye(ns)-A       -B;
              C                D];
    rank(S1);
end 

S2  =[  A    B; 
        -C   -D];
I   = [eye(size(A))  zeros(ns,2);
     zeros(2,ns) zeros(2,2)];
[v2 d2] = eig(S2,I);
%% Task 2_p2:
figure
plot(eig(A) ,'x')
hold on
plot(diag(d2),'o')
xlim([-3 2])
xlabel('real')
ylabel('complex')
axis equal
grid on
circle(1)
title('Eigenvalues and Transmission zeros of the 7-State Model')
legend('Eigenvalues','Transmission Zeros');
%% Task 2_p3:
exp(eig(A)*ts)

%% Task2_p4:


for jj = 1:4

B1  = B(:,1);
B2  = B(:,2);
C1  = C(1,:);
C2  = C(2,:);
D1  = [0];


%     x   = 1:cl;
%     [U,S,V]=svd(H);
%     [T d] = eig(H);
%     Si  =diag(S);
%     Si  = diag(Si(1:ns));
%     U   = U(:,1:ns);
%     V   = V(:,1:ns);
%     O   = U;
%     C   = Si*V';

% Hh = zeros(cl*2);
% for i = 1:cl
%     for j = 1:cl
%         k=i+j;
%         Hh(2*i-1,j*2-1)=y(k+bi); %hk(1,1)
%         Hh(2*i-1,j*2)=y(k+bi);   %hk(1,2)
%         Hh(2*i,j*2-1)=y(k+bi); %hk(2,1)
%         Hh(2*i,j*2)=y(k+bi); %hk(2,2)
%     end
% end 
% Hh = downsample(Hh,2);
% Hh = downsample(Hh',2);
% Hh = Hh';
% Define New A,B,C, and D

    H   = zeros(2*cl);
    cl = 100;
switch jj
    case 1
        S2=[ A      B1; 
             -C1    -D1];
        for i = 1:cl
        for j   = 1:cl
            k   = i+j-1;
            H(i,j)  =y11(k+bi);
        end
        end 
    case 2
        S2=[ A      B2;
            -C1     -D1];
        for i = 1:cl
        for j   = 1:cl
            k   = i+j-1;
            H(i,j)  =y12(k+bi);
        end
        end 
    case 3
        S2 =[A      B1;
             -C2    -D1];
        for i = 1:cl
        for j   = 1:cl
            k   = i+j-1;
            H(i,j)  =y21(k+bi);
        end
        end 
    case 4
        S2 = [A     B2;
             -C2    -D1];
        for i = 1:cl
        for j   = 1:cl
            k   = i+j-1;
            H(i,j)  =y22(k+bi);
        end
        end 
end 

subplot(2,2,jj)
plot (svd(H),'*')
title('Hankel Singular Value for Each Channel')
grid on
xlim([0 40])
hold on


I = [eye(size(A))  zeros(ns,1);
     zeros(1,ns)    zeros(1,1)];
[v2 d2] = eig(S2,I);

subplot(2,2,jj)
plot(eig(A),'x','LineWidth',1);
hold on
plot(diag(d2),'o','LineWidth',1);
xlim([-3 2])
ylim([-1.5 1.5])
xlabel('real')
ylabel('complex')
title(['plot ' num2str(jj) ', when ns = ' num2str(ns)])
axis equal
grid on
circle(1);
hold on
end 
%% Task4_p1
load u_rand.mat
y1 = u_rand.Y(3).Data;
y2 = u_rand.Y(4).Data;
u1 = u_rand.Y(1).Data;
u2 = u_rand.Y(2).Data;

N = length(y1);
t = [0:N-1]*ts - 1;

p = 12000;
u = [u1; u2];
y = [y1; y2];

mean(u1);
mean(u2);

%% Task4_p2:
% a = zeros(1,401);b = zeros(1,401);c = zeros(1,401);d = zeros(1,401);
% 
% for k = 1:401
%     sm = zeros(2);
%     for q = 202:2*p-202
%           sm = sm + u(:,q+k-201)*(u(:,q))' ;   
%     end
%     matr = 1/(2*p)*sm;
%     a(k) = matr(1,1);b(k) = matr(1,2);c(k) = matr(2,1);d(k) = matr(2,2);
% end
% 
% figure
% x = linspace(-5,5,401);
% subplot(2,2,1)
% plot (x,a,'*','LineWidth',1)
% ylim([-1 5])
% title('Ruu Estimation for u11 Channel')
% grid on
% subplot(2,2,2)
% plot (x,b,'*','LineWidth',1)
% ylim([-1 5])
% title('Ruu Estimation for u12 Channel')
% grid on
% subplot(2,2,3)
% plot (x,c,'*','LineWidth',1)
% ylim([-1 5])
% title('Ruu Estimation for u21 Channel')
% grid on
% subplot(2,2,4)
% plot (x,d,'*','LineWidth',1)
% ylim([-1 5])
% title('Ruu Estimation for u22 Channel')
% grid on

%% Task4_p3
Ruu = [0 0 ; 0 0];
for q = 1:24001
    Ruu = (Ruu + u(:,q)*(u(:,q))');
end 
(1/(2*p))*Ruu

%% Task4_p4
% a = zeros(1,89);b = zeros(1,89);c = zeros(1,89);d = zeros(1,89);
% 
% for k = 1:89
%     sm = zeros(2);
%     for q = 89:2*p-89
%           sm = sm + y(:,q+k-8)*(u(:,q))' ;   
%     end
%     matr = 1/(2*p)*sm;
% 
%     a(k) = matr(1,1)/var(u1);b(k) = matr(1,2)/var(u2);c(k) = matr(2,1)/var(u1);d(k) = matr(2,2)/var(u2);
% end
% 
% x = linspace(-0.2,2,89);
% figure
% subplot(4,1,1)
% plot (x,a,'*')
% hold on
% grid on
% 
% subplot(4,1,2)
% plot (x,c,'*')
% hold on
% grid on
% 
% subplot(4,1,3)
% plot (x,b,'*')
% hold on
% grid on
% 
% subplot(4,1,4)
% plot (x,d,'*')
% hold on
% grid on

%% Task5_p1

ns  = 7;
cl  = 100;
bi  = 41;
H   = zeros(cl*2);
for i = 1:cl
    for j = 1:cl
                        k=i+j-1;
                        H(2*i-1,j*2-1)=y11(k+bi); %hk(1,1)
                        H(2*i-1,j*2)=y12(k+bi);   %hk(1,2)
                        H(2*i,j*2-1)=y21(k+bi); %hk(2,1)
                        H(2*i,j*2)=y22(k+bi);
    end
end

x       = 1:cl;
[U,S,V] =svd(H);
[T d]   = eig(H);
Si      = zeros(cl*2);
for i = 1:ns
    Si(i,i) = S(i,i);
end 
O = U;
C = Si*V';

Hh = zeros(cl*2);
for i = 1:cl
    for j = 1:cl
                        k=i+j;
                        Hh(2*i-1,j*2-1)=y11(k+bi); %hk(1,1)
                        Hh(2*i-1,j*2)=y12(k+bi);   %hk(1,2)
                        Hh(2*i,j*2-1)=y21(k+bi); %hk(2,1)
                        Hh(2*i,j*2)=y22(k+bi); %hk(2,2)
    end
end 
U       = U(:,1:ns);
V       = V(:,1:ns);
Si      = Si(1:ns,1:ns);

A       = U'*Hh*(V*inv(Si));
B       = C(1:ns,1:2);
C       = O(1:2,1:ns);
D       = zeros(2);

load u_rand.mat
y1 = u_rand.Y(3).Data/2;
y2 = u_rand.Y(4).Data/2;
u1 = u_rand.Y(1).Data/2;
u2 = u_rand.Y(2).Data/2;

ts = 1/40;
N = length(y1);
t = [0:N-1]*ts - 1;

p = 12000;
u = [u1; u2];
y = [y1; y2];

%%Problem 1
Ruu     = [0 0 ; 0 0];
for q = 1:24001
    Ruu = Ruu + y(:,q)*(y(:,q))';
end 
y_rms   = sqrt(diag((1/(2*p))*Ruu));
RMS     = norm(y_rms);

%% Task5_p2
sys     = ss(A,B,C,D,ts);
Gc      = gram(sys,'c');
Go      = gram(sys,'o');
PH2_1   = norm(sqrt(diag(B'*Go*B)));
PH2_2   = norm(sqrt(diag(C*Gc*C')));

%% Task5_p3
%experimental norm
 PH2_3  = 0;
for i = 1:401
    n = norm([y11(i) y12(i);y21(i) y22(i)],'fro');
    PH2_3 = PH2_3+n^2;
end 
PH2_3   = sqrt(PH2_3);
PH2_3   = norm(PH2_3);

%% Task6_1
ll  = 0.0000; %lower limit
ul  = 5; %upper limit
tolerance = 1e-5;

[gam_c,eigd]= hinfnormc(ll,ul,tolerance,A,B,C,D);
[gam_d,fre] = hinfnormd(ll,ul,tolerance,A,B,C,D,ts);

%Hankel singular value
max(eig(Gc*Go));

%% Task6_4
% Singular value of the frequency response.
    frs_s   = cell2mat(frs_s);
    frs_e   = cell2mat(frs_e);
    figure
    xaxis   = linspace(1,20*2*pi,205);
    plot(xaxis,frs_s(1,:),'*')
    hold on
    plot(xaxis,frs_e(1,:),'*')
    title('Singular Value of Frequency Response')
    grid on
    legend('Empirical Data','Simulated Data')
    xlabel('Frequency (Hz)')
    ylabel('Value')

toc
