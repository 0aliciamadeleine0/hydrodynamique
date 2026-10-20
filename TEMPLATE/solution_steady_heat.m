clear all; close all;
addpath ../COMMON/;

% Visual aspects
set(0,'defaultAxesFontSize',20);
set(0,'defaultAxesLineWidth',2);

% Physical parameters
kappa=0.005;
amplitude_s=0.0; % 2.0
Thole=5;
Tborder=2;
Tcolormap=[0 5];

% Read the informations from the geometry
x=readarray('x.dat');
y=readarray('y.dat');
mesh_info=load('meshinfo.dat');

% Read matrices from the Finite Element Formulation
% Read mass matrix
M=readmat('M.dat');
[n,m]=size(M);
%pcolor(M);
%
K=readmat('K.dat');
%
% Get the position of the Dirichlet masks
BC_Dir_hole=readarray('BCdir5.dat');
ii_Dir_hole=find(abs(BC_Dir_hole)>0);
BC_Dir_borders=readarray('BCdir1234.dat');
ii_Dir_borders=find(abs(BC_Dir_borders)>0);
BC_Dir=readarray('BCdir.dat');
ii_int=find(BC_Dir==0);
%
% Build corresponding matrices
A=kappa*K;
Ared=A(ii_int,ii_int);
s=amplitude_s*ones(n,1);

% SOLVE
Td=zeros(n,1);
Td(ii_Dir_hole)=Thole;
Td(ii_Dir_borders)=Tborder;
%
C=M*s-A*Td;
Cred=C(ii_int);
Tred=Ared\Cred;
%
T=Td;
T(ii_int)=Tred;
%

%FIG=visu_2D(x,y,mesh_info,T,'black','on');
%FIG=visu_2D(x,y,mesh_info,T,'none','off');
FIG=visu_2D(x,y,mesh_info,T,'none','on');
caxis(Tcolormap);
%
savename=['steady_solution'];
exportgraphics(FIG,[savename '.png']);


