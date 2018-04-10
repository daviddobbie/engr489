%% David Dobbie
% Victoria University of Wellington
% Recreating paper 2 (Mellin Transform of CPMG data)
%
% 
% L. Venkataramanan et al/ journal of Magnetic Resonance 206 (2010) 20-31

%Aim: Computing moments of transverse T2 relaxation time from measured CPMG
%       data. Will implement in one dimension.


close all
clc
clf
clear all

%% init variables
% number of data points in each dimension
N2 = 1000;

% number of bins in relaxation time grids
Ny = 101;      

%sets how we compress the data points.
% NOTE that:
%    N2 > trunc2 >= Ny
%    N1 > trunc1 >= Nx
trunc2=25;

deltatau2 = 3.5e-4;
T2 = logspace(-3,1,Ny); %form T2 domain, use log since will be small

%forms measurement arrays, time tau1 and tau2 domains
tau2 = (1:N2)'*deltatau2;  

K2 = exp(-tau2 * (1./T2) );     % simple T2 relaxation data
