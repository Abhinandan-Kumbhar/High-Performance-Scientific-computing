%%%%%%%%%%%%%%%%%%%%%% Guyan Reduction %%%%%%%%%%%%%%%%%%%%
format compact;format short;clc;clear all;
load("Project1_2_GMAM_GMTA_Matlab.mat");
%% Stiffness inertia ratio %%
disp('Stiffness inertia ratio')
Y = diag(K)./diag(M);
% sorting in increasing order %
[y I] = sort(Y);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Number of master degress of freedom %%
mdof = 500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I(mdof)=I(441);
%% Constructing Permutation matrix for Partioning of Mass stiffness matrix %%
P = zeros(size(M));
for i = 1:mdof
    P(I(i),i) = 1;
end
for j = (mdof+1):max(size(M))
    P(I(j),j) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w=30000;
t=0:0.01:1;
L=zeros(length(M),length(t));
for i=1:length(t)
    L(3781,i)=100000*sin(w*t(i));
    L(3782,i)=100000*sin(w*t(i));
end
%% Constructing Partioned matrix %%
Mr = P'*M*P;
Kr = P'*K*P;
Lr = P'*L;
Mmm = Mr(1:mdof,1:mdof);
Mms = Mr(1:mdof,(mdof+1):max(size(M)));
Msm = Mr((mdof+1):max(size(M)),1:mdof);
Mss = Mr((mdof+1):max(size(M)),(mdof+1):max(size(M)));
Kmm = Kr(1:mdof,1:mdof);
Kms = Kr(1:mdof,(mdof+1):max(size(M)));
Ksm = Kr((mdof+1):max(size(M)),1:mdof);
Kss = Kr((mdof+1):max(size(M)),(mdof+1):max(size(M)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Constructing Transformation matrix %%
RG = -Kss\Ksm;
TG = [eye(mdof);RG];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Guyan reduction of system %%
MG = TG'*Mr*TG;
KG = TG'*Kr*TG;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LG=TG'*Lr;
modes=10;

%% Approximate Natural-Frequencies %%
[Vapprox,Eapprox] = eigs(KG,MG,10,'smallestabs');
omegaGR = sqrt(diag(Eapprox));
disp('Natural frequencies using Guyan condensation')
disp(diag(Eapprox'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exact Natural-Frequencies %%
Eexact = eigs(K,M,modes,'smallestabs','Tolerance',1e-15);
omegaE = sqrt(Eexact);
disp('Exac Natural frequencies')
disp(Eexact')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Error in approximation %%
e = abs(omegaGR-omegaE);
percente = abs(e./omegaE)*100;
disp('Absolute error in Natural frequencies')
disp(e)
disp('Percent error in Natural frequencies')
disp(percente)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% to determine position of load vactors in condensed matrix
for i=1:length(I)
    if I(i)==3782
        disp('3782 present in ')
        disp(i)
    elseif I(i)==3781
    disp('3781 present in')
    disp(i)
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
