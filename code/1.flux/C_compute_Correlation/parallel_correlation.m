%CORR1
combination='C951-1000'
%CORR2
combination='C901-950'
%CORR3
combination='C851-900'
%CORR4
combination='C801-850'
%CORR5
combination='C751-800'
%CORR6
combination='C701-750'
%CORR7
combination='C651-700'
%CORR8
combination='C601-650'
%CORR9
combination='C551-600'
%CORR10
combination='C51-100'

P1_data1=[];P2_data1=[];P3_data1=[];P4_data1=[];P5_data1=[];P6_data1=[];P7_data1=[];P8_data1=[];P9_data1=[];P10_data1=[];
R1_data1=[];R2_data1=[];R3_data1=[];R4_data1=[];R5_data1=[];R6_data1=[];R7_data1=[];R8_data1=[];R9_data1=[];R10_data1=[];
P1_data1_secretor=[];P2_data1_secretor=[];P3_data1_secretor=[];P4_data1_secretor=[];P5_data1_secretor=[];P6_data1_secretor=[];P7_data1_secretor=[];P8_data1_secretor=[];P9_data1_secretor=[];P10_data1_secretor=[];
R1_data1_secretor=[];R2_data1_secretor=[];R3_data1_secretor=[];R4_data1_secretor=[];R5_data1_secretor=[];R6_data1_secretor=[];R7_data1_secretor=[];R8_data1_secretor=[];R9_data1_secretor=[];R10_data1_secretor=[];
P1_data2=[];P2_data2=[];P3_data2=[];P4_data2=[];P5_data2=[];P6_data2=[];P7_data2=[];P8_data2=[];P9_data2=[];P10_data2=[];
R1_data2=[];R2_data2=[];R3_data2=[];R4_data2=[];R5_data2=[];R6_data2=[];R7_data2=[];R8_data2=[];R9_data2=[];R10_data2=[];
P1_data2_secretor=[];P2_data2_secretor=[];P3_data2_secretor=[];P4_data2_secretor=[];P5_data2_secretor=[];P6_data2_secretor=[];P7_data2_secretor=[];P8_data2_secretor=[];P9_data2_secretor=[];P10_data2_secretor=[];
R1_data2_secretor=[];R2_data2_secretor=[];R3_data2_secretor=[];R4_data2_secretor=[];R5_data2_secretor=[];R6_data2_secretor=[];R7_data2_secretor=[];R8_data2_secretor=[];R9_data2_secretor=[];R10_data2_secretor=[];

load(['../../data/models_split/solCombi_',combination])
load('Network_Red_ext');

Impossible_model = [];
tic



for z=1:size(solCombi,2)
[P1_D1 P2_D1 P3_D1 P4_D1 P5_D1 P6_D1 P7_D1 P8_D1 P9_D1 P10_D1 ...
 P1_D1_S P2_D1_S P3_D1_S P4_D1_S P5_D1_S P6_D1_S P7_D1_S P8_D1_S P9_D1_S P10_D1_S... 
 P1_D2 P2_D2 P3_D2 P4_D2 P5_D2 P6_D2 P7_D2 P8_D2 P9_D2 P10_D2... 
 P1_D2_S P2_D2_S P3_D2_S P4_D2_S P5_D2_S P6_D2_S P7_D2_S P8_D2_S P9_D2_S P10_D2_S ...
 R1_D1 R2_D1 R3_D1 R4_D1 R5_D1 R6_D1 R7_D1 R8_D1 R9_D1 R10_D1 ...
 R1_D1_S R2_D1_S R3_D1_S R4_D1_S R5_D1_S R6_D1_S R7_D1_S R8_D1_S R9_D1_S R10_D1_S ...
 R1_D2 R2_D2 R3_D2 R4_D2 R5_D2 R6_D2 R7_D2 R8_D2 R9_D2 R10_D2 ...
 R1_D2_S R2_D2_S R3_D2_S R4_D2_S R5_D2_S R6_D2_S R7_D2_S R8_D2_S R9_D2_S R10_D2_S Imp]=computeCorr(z,solCombi(:,z),Network_Red,false);


P1_data1(:,z)= P1_D1;
P2_data1(:,z)= P2_D1;
P3_data1(:,z)= P3_D1;
P4_data1(:,z)= P4_D1;
P5_data1(:,z)= P5_D1;
P6_data1(:,z)= P6_D1;
P7_data1(:,z)= P7_D1;
P8_data1(:,z)= P8_D1;
P9_data1(:,z)= P9_D1;
P10_data1(:,z)= P10_D1;

P1_data1_secretor(:,z)= P1_D1_S;
P2_data1_secretor(:,z)= P2_D1_S;
P3_data1_secretor(:,z)= P3_D1_S;
P4_data1_secretor(:,z)= P4_D1_S;
P5_data1_secretor(:,z)= P5_D1_S;
P6_data1_secretor(:,z)= P6_D1_S;
P7_data1_secretor(:,z)= P7_D1_S;
P8_data1_secretor(:,z)= P8_D1_S;
P9_data1_secretor(:,z)= P9_D1_S;
P10_data1_secretor(:,z)= P10_D1_S;

P1_data2(:,z)= P1_D2;
P2_data2(:,z)= P2_D2;
P3_data2(:,z)= P3_D2;
P4_data2(:,z)= P4_D2;
P5_data2(:,z)= P5_D2;
P6_data2(:,z)= P6_D2;
P7_data2(:,z)= P7_D2;
P8_data2(:,z)= P8_D2;
P9_data2(:,z)= P9_D2;
P10_data2(:,z)= P10_D2;

P1_data2_secretor(:,z)= P1_D2_S;
P2_data2_secretor(:,z)= P2_D2_S;
P3_data2_secretor(:,z)= P3_D2_S;
P4_data2_secretor(:,z)= P4_D2_S;
P5_data2_secretor(:,z)= P5_D2_S;
P6_data2_secretor(:,z)= P6_D2_S;
P7_data2_secretor(:,z)= P7_D2_S;
P8_data2_secretor(:,z)= P8_D2_S;
P9_data2_secretor(:,z)= P9_D2_S;
P10_data2_secretor(:,z)= P10_D2_S;

R1_data1(:,z)= R1_D1;
R2_data1(:,z)= R2_D1;
R3_data1(:,z)= R3_D1;
R4_data1(:,z)= R4_D1;
R5_data1(:,z)= R5_D1;
R6_data1(:,z)= R6_D1;
R7_data1(:,z)= R7_D1;
R8_data1(:,z)= R8_D1;
R9_data1(:,z)= R9_D1;
R10_data1(:,z)= R10_D1;

R1_data1_secretor(:,z)= R1_D1_S;
R2_data1_secretor(:,z)= R2_D1_S;
R3_data1_secretor(:,z)= R3_D1_S;
R4_data1_secretor(:,z)= R4_D1_S;
R5_data1_secretor(:,z)= R5_D1_S;
R6_data1_secretor(:,z)= R6_D1_S;
R7_data1_secretor(:,z)= R7_D1_S;
R8_data1_secretor(:,z)= R8_D1_S;
R9_data1_secretor(:,z)= R9_D1_S;
R10_data1_secretor(:,z)= R10_D1_S;

R1_data2(:,z)= R1_D2;
R2_data2(:,z)= R2_D2;
R3_data2(:,z)= R3_D2;
R4_data2(:,z)= R4_D2;
R5_data2(:,z)= R5_D2;
R6_data2(:,z)= R6_D2;
R7_data2(:,z)= R7_D2;
R8_data2(:,z)= R8_D2;
R9_data2(:,z)= R9_D2;
R10_data2(:,z)= R10_D2;

R1_data2_secretor(:,z)= R1_D2_S;
R2_data2_secretor(:,z)= R2_D2_S;
R3_data2_secretor(:,z)= R3_D2_S;
R4_data2_secretor(:,z)= R4_D2_S;
R5_data2_secretor(:,z)= R5_D2_S;
R6_data2_secretor(:,z)= R6_D2_S;
R7_data2_secretor(:,z)= R7_D2_S;
R8_data2_secretor(:,z)= R8_D2_S;
R9_data2_secretor(:,z)= R9_D2_S;
R10_data2_secretor(:,z)= R10_D2_S;

Impossible_model = [Impossible_model Imp];

end
toc

save (['P_data1_',combination], 'P1_data1','P2_data1','P3_data1','P4_data1','P5_data1','P6_data1','P7_data1','P8_data1','P9_data1','P10_data1')
save (['P_data1_secretor_',combination], 'P1_data1_secretor','P2_data1_secretor','P3_data1_secretor','P4_data1_secretor','P5_data1_secretor','P6_data1_secretor','P7_data1_secretor','P8_data1_secretor','P9_data1_secretor','P10_data1_secretor') 
save (['P_data2_',combination], 'P1_data2','P2_data2','P3_data2','P4_data2','P5_data2','P6_data2','P7_data2','P8_data2','P9_data2','P10_data2') 
save (['P_data2_secretor_',combination], 'P1_data2_secretor','P2_data2_secretor','P3_data2_secretor','P4_data2_secretor','P5_data2_secretor','P6_data2_secretor','P7_data2_secretor','P8_data2_secretor','P9_data2_secretor','P10_data2_secretor') 

save (['R_data1_',combination], 'R1_data1','R2_data1','R3_data1','R4_data1','R5_data1','R6_data1','R7_data1','R8_data1','R9_data1','R10_data1') 
save (['R_data1_secretor_',combination], 'R1_data1_secretor','R2_data1_secretor','R3_data1_secretor','R4_data1_secretor','R5_data1_secretor','R6_data1_secretor','R7_data1_secretor','R8_data1_secretor','R9_data1_secretor','R10_data1_secretor') 
save (['R_data2_',combination], 'R1_data2','R2_data2','R3_data2','R4_data2','R5_data2','R6_data2','R7_data2','R8_data2','R9_data2','R10_data2') 
save (['R_data2_secretor_',combination], 'R1_data2_secretor','R2_data2_secretor','R3_data2_secretor','R4_data2_secretor','R5_data2_secretor','R6_data2_secretor','R7_data2_secretor','R8_data2_secretor','R9_data2_secretor','R10_data2_secretor') 

save (['model_impossible_',combination], 'Impossible_model')
