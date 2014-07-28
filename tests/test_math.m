clear all; clc; close all;

% format short e
format long e

disp('*********************************************************************');
disp('***************************** 3D Tests ******************************');
disp('*********************************************************************');
rng('default');
V_1 = randn(1,3)
V_2 = randn(1,3)
M_1 = randn(3,3)
M_2 = randn(3,3)
M_2_DET = det(M_2)
M_2 = M_2'
M_3 = M_1 * M_2
V_3 = V_1 + V_2
V_4 = V_1 - V_2
V_5 = V_1 .* V_2
VOT = dot(V_1, V_2)
M_4 = M_3^(-1)

%% Try eigenvector and Eigenvalue calculation
rng('default');
M_1 = randn(3,3);
M_1 = M_1*M_1'
[EigenVectors,eigenValues] = eig(M_1);  %% the eigen vectors are colums in the return matrix
eigenValues = [eigenValues(1,1) eigenValues(2,2) eigenValues(3,3)]
eigenVector0 = [EigenVectors(1,1) EigenVectors(2,1) EigenVectors(3,1)]
eigenVector1 = [EigenVectors(1,2) EigenVectors(2,2) EigenVectors(3,2)]
eigenVector2 = [EigenVectors(1,3) EigenVectors(2,3) EigenVectors(3,3)]

BIG_N = 10;
inA = zeros(BIG_N, BIG_N);
k1 = 1;
k2 = 1;
for i = 1:BIG_N
    for j = i:BIG_N
        inA(i,j) = k1;
        nextK = k1 + k2;
        k1 = k2;
        k2 = nextK;
        inA(j,i) = inA(i,j);
    end
end
inA
[EigenVectors,eigenValues] = eig(inA);  %% the eigen vectors are colums in the return matrix
for i = 1:length(inA)
    eigs(i) = eigenValues(i,i);
end
eigs

%% 4D Calculations

disp('*********************************************************************');
disp('***************************** 4D Tests ******************************');
disp('*********************************************************************');
rng('default');
V4D_1 = randn(1,4)
V4D_2 = randn(1,4)
M4D_1 = randn(4,4)
M4D_2 = randn(4,4)
M4D_2_DET = det(M4D_2)
M4D_2 = M4D_2'
M4D_3 = M4D_1 * M4D_2
V4D_3 = V4D_1 + V4D_2
V4D_4 = V4D_1 - V4D_2
V4D_5 = V4D_1 .* V4D_2
V4DOT = dot(V4D_1, V4D_2)
M4D_4 = M4D_3^(-1)

%% 2D Calculations

disp('*********************************************************************');
disp('***************************** 2D Tests ******************************');
disp('*********************************************************************');
rng('default');
V_1 = randn(1,2)
V_2 = randn(1,2)
M_1 = randn(2,2)
M_2 = randn(2,2)
M_2_DET = det(M_2)
M_2 = M_2'
M_3 = M_1 * M_2
V_3 = V_1 + V_2
V_4 = V_1 - V_2
V_5 = V_1 .* V_2
VOT = dot(V_1, V_2)
M_4 = M_3^(-1)

disp('*********************************************************************');
disp('************************* Quaternion Tests **************************');
disp('*********************************************************************');
rng('default');
axis = [randn(1), randn(1), randn(1)];
length = sqrt(sum(axis.*axis));
axis = axis./length
angle = randn(1)
Q1 = spin_calc('EVtoQ',[axis (360*angle/(2*pi))],eps,0) 
RotMat1 = spin_calc('QtoDCM',Q1,eps,1)' 
RotMat2 = RotMat1;
RotMat2(4,1) = 0; RotMat2(4,2) = 0; RotMat2(4,3) = 0; RotMat2(4,4) = 1;
RotMat2(1,4) = 0; RotMat2(2,4) = 0; RotMat2(3,4) = 0
Q2 = Q1
Q3 = Q1