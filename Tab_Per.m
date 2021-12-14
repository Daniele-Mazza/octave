clear all;clc;format short;format compact;
d = 0.18; % clearance

figure(1,'position',[100 100 1000 800]);
axis([0 18 0 18 0 10]);axis on;grid on;grid minor on;xticks(20);yticks(20);

AtNo = [0 18;17 18;0 17;1 17;12 17;13 17;14 17;15 17;16 17;17 17;... % x(i),y(i)
        0 16;1 16;12 16;13 16;14 16;15 16;16 16;17 16;...
        0 15;1 15;2 15;3 15;4 15;5 15;6 15;7 15;8 15;9 15;10 15;11 15;12 15;13 15;14 15;15 15;16 15;17 15;...
        0 14;1 14;2 14;3 14;4 14;5 14;6 14;7 14;8 14;9 14;10 14;11 14;12 14;13 14;14 14;15 14;16 14;17 14;...
        0 13;1 13;2 10;3 10;4 10;5 10;6 10;7 10;8 10;9 10;10 10;11 10;12 10;13 10;14 10;15 10;...
        2 13;3 13;4 13;5 13;6 13;7 13;8 13;9 13;10 13;11 13;12 13;13 13;14 13;15 13;16 13;17 13;...
        0 12;1 12;2 9;3 9;4 9;5 9;6 9;7 9;8 9;9 9;10 9;11 9;12 9;13 9;14 9;15 9;...
        2 12;3 12;4 12;5 12;6 12;7 12;8 12;9 12;10 12;11 12;12 12;13 12;14 12;15 12;16 12;17 12];

A = importdata("electronegativity.txt","\n");
ln = length(A)
for i=1:ln
  b = strsplit(A{i});
  nAt = str2num(b{1});eNeg = str2num(b{4});h1 = eNeg;
  x1 = AtNo(nAt,1) + d;y1 = AtNo(nAt,2) - d;x2 = x1 + 1 - d;y2 = y1 -1 + d;
  v2 = [x1 y1 0;x2 y1 0;x2 y2 0;x1 y2 0;x1 y1 h1;x2 y1 h1;x2 y2 h1;x1 y2 h1];  
  f2 = [3 4 8 7;1 4 8 5;5 6 7 8];
  patch('Faces',f2,'Vertices',v2,'FaceColor',[ (h1*0.25) (1-h1*0.25) 0]);view(3);
endfor  

xlabel("groups",'FontSize',20);ylabel("periods",'FontSize',20);zlabel("Elecronegativity",'FontSize',20);
title("Electronegativity , Pauling scale",'FontSize',20);

for i=0.5:17.5
  text(i,7,num2str(i+0.5),'FontSize',20);   
endfor  
for i=18:-1:12
  text(-1.3,i-1,num2str(19-i),'FontSize',20);   
endfor  