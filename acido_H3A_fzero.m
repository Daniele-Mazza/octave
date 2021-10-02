%  H3PO4 <==> H+  +  H2PO4-      Ka1 = [H+]*[H2PO4-]/[H3PO4]
%  H2PO4- <==> H+  +  HPO4--     Ka2 = [H+]*[HPO4--]/[H2PO4-]
%  HPO4-- <==> H+  +  PO4---     Ka3 = [H+]*[PO4---]/[HPO4--]
clear;clc;
Atot = 0.4;    % total concentration of H3PO4  (i.e. total amoun of H3PO4 added)
Natot = 0.4;   % total concentration of Na3PO4 (i.e. total amoun of Na3PO4 added)
global Ka1 = 7.5e-3;  % 1st dissociation constant
global Ka2 = 6.2e-8;  % 2nd dissociation constant
global Ka3 = 4.4e-13; % 3rd dissociation constant
global Kw = 1e-14;    % ionic water product
global Na = Natot*3;
global Ptot = Natot + Atot;
global H X
function y = neut(pH)
    global Ptot Na Ka1 Ka2 Ka3 Kw H X
    H = 10^(-pH);
    OH = Kw/H;
    % here we solve the system of 4 equations, 4 unknowns
    A = [-Ka1 H 0 0 ;     %  -Ka1*H3PO4 + H*H2PO4 + 0*HPO4 + 0*PO4 = 0
          0 -Ka2 H 0;     %   0*H3PO4 - Ka2*H2PO4 + H*HPO4 + 0*PO4 = 0
          0 0 -Ka3 H;     %   0*H3PO4 + 0*H2PO4 - Ka3*HPO4 + H*PO4 = 0
          1  1  1  1];    %   1*H3PO4 + 1*H2PO4 +  1*HPO4 +  1*PO4 = Ptot
    b =[0; 0; 0; Ptot];     
    X = A\b;        %   X(1)=H3PO4   X(2)=H2PO4-   X(3)=HPO4--   X(4)=PO4---
    y = H + Na - Kw/H - X(2) - X(3)*2 - X(4)*3;
    endfunction
[pH,fval,info] = fzero(@neut,[0,14])
T1 = ['H3PO4 initial concentration   = ',num2str(Atot,5)];disp(T1);
T1 = ['Na3PO4 initial concentration  = ',num2str(Natot,5)];disp(T1);
disp('----------------------------------');

T1 = ['pH        = ',num2str(pH,6)];disp(T1);
T1 = ['[H+]      = ',num2str(H,6)];disp(T1);
T1 = ['[OH-]     = ',num2str(Kw/H,6)];disp(T1);
T1 = ['[H3PO4]   = ',num2str(X(1),6)];disp(T1);
T1 = ['[H2PO4-]  = ',num2str(X(2),6)];disp(T1);
T1 = ['[HPO4--]  = ',num2str(X(3),6)];disp(T1);
T1 = ['[PO4---]  = ',num2str(X(4),6)];disp(T1);
disp('----------------------------------')






