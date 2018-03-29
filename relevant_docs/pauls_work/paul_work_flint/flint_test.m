% Example program to run regularised non-negative least squares
% Copyright Victoria University of Wellington
% P Teal, Thursday, 15 August 2013
% Tuesday, 6 May 2014

clear all


N1 = 50;       % number of data points in each dimension
N2 = 10000;

Nx = 100;      % number of bins in relaxation time grids
Ny = 101;      

tau1min = 1e-4;
tau1max = 10;
deltatau2 = 3.5e-4;

T1 = logspace(-2,1,Nx);
T2 = logspace(-2,1,Ny);

[T2a,T1a] = meshgrid(log10(T2),log10(T1));

tau2 = (1:N2)'*deltatau2;  
% measurement time array (note, no t=0 data available)

tau1 = logspace(log10(tau1min),log10(tau1max),N1)';

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation data
K1 = 1-2*exp(-tau1 *(1./T1) );  % T1 relaxation data

if 1==1
  % Generate some synthetic data
  Ftrue = zeros(Nx,Ny);
  
  centre1 = [-0.5 -1];
  radius1 = [0.35 0.75];
  centre2 = [-1.5 0.3];
  radius2 = 0.2;
  centre3 = [0.5 0.3];
  radius3 = 0.2;
  
  dist1 = sqrt( (T2a-centre1(1)).^2 + (T1a-centre1(2)).^2 );
  dist2 = sqrt( (T2a-centre2(1)).^2 + (T1a-centre2(2)).^2 );
  dist3 = sqrt( (T2a-centre3(1)).^2 + (T1a-centre3(2)).^2 );
  ang1 = atan2( T1a-centre1(2), T2a-centre1(1));
  
  f1 = find(radius1(1) < dist1 & dist1<radius1(2) & ang1<=0);
  Ftrue(f1) = 0.4;
  f2 = find(dist2<radius2);
  Ftrue(f2) = 0.18;
  f3 = find(dist3<radius3);
  Ftrue(f3) = 0.58;
  
  %Ftrue = rand(size(Ftrue));
  %Ftrue = max(0,randn(size(Ftrue))-2);
end

Mmeas_cln = K1*Ftrue*K2';   % without noise
sigma = max(Mmeas_cln(:))*1e-2; sigma = sigma*1e-1;
noise = sigma*randn(N1,N2);
Mmeas = Mmeas_cln + noise;  % with noise

%Mmeas = Mmeas * 1e6;

reg = 1e-3;
alpha = sigma*reg;
alpha = alpha*1e5;

alpha = 1e-1;

%[Sflint,resa,sss] = flint1(K1,K2,Mmeas,alpha);
[Sflint,resa] = flint1(K1,K2,Mmeas,alpha);
%Sflint = flint1(K1,K2,Mmeas,alpha);


figure(1);
mesh(T2a,T1a,Ftrue);
xlabel('y (T2)');
ylabel('x (T1)');
view(330,75);

figure(2);
mesh(T2a,T1a,Sflint);
xlabel('y (T2)');
ylabel('x (T1)');
view(330,75);

figure(3);
semilogy(1:length(resa),resa,'+',1:length(resa),resa-min(resa));
zoom on;
grid on;

%figure(4);
%semilogy(1:length(sss),sss./resa);
%zoom on;
%grid on;
