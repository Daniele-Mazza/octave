%  H3PO4 <==> H+  +  H2PO4-
%  Ka1 = [H+]*[H2PO4-]/[H3PO4]
%  H2PO4- <==> H+  +  HPO4--
%  Ka2 = [H+]*[HPO4--]/[H2PO4-]
%  HPO4-- <==> H+  +  PO4---
%  Ka3 = [H+]*[PO4---]/[HPO4--]
clc;clear;format compact;
Atot = 0.4;    % total concentration of H3PO4  (i.e. total amount of H3PO4 added)
CNaOH = 1.00;  % molar concentration of added sodium hydroxide
Ka1 = 7.5e-3;  % 1st dissociation constant
Ka2 = 6.2e-8;  % 2nd dissociation constant
Ka3 = 4.4e-13; % 3rd dissociation constant
Kw = 1e-14;    % ionic product of water
T1 = ['H3PO4 initial concentration   = ',num2str(Atot,5)];disp(T1);
disp('----------------------------------');
k = 0;
for mL_NaOH = 0:10:1000 % mLitre NaOH 1M added to 1 litre solution of H3PO4
A1 = Atot*1000/(1000 + mL_NaOH);
Na = mL_NaOH*CNaOH/(1000 + mL_NaOH);
pH1 = 0;
pH2 = 14;
pHstep = 1;
for j = (1:8)     % each step of the loop pHstep decreases 10x
  for pH = (pH1:pHstep:pH2)           % loop of pH
    H = 10^(-pH);
    % here we solve the system of 4 equations, 4 unknowns (speed up needed!)
    %  Ka1 = H*H2PO4/H3PO4
    %  Ka2 = H*HPO4/H2PO4
    %  Ka3 = H*PO4/HPO4
    %  A1 = H3PO4 + H2PO4 + HPO4 + PO4
    PO4 = Ka1*Ka2*Ka3*A1/(H*H*H + Ka1*H*H + Ka1*Ka2*H + Ka1*Ka2*Ka3);
    HPO4 = H*PO4/Ka3;
    H2PO4 = H*HPO4/Ka2;
    H3PO4 = H*H2PO4/Ka1;
    Neut = H + Na - Kw/H - H2PO4 - HPO4*2 - PO4*3;
    if Neut<0;break;end  
  end
  pH2 = pH;
  pH1 = pH2 - pHstep;
  pHstep = pHstep/10;
end
T1 = ['NaOH ml added  = ',num2str(mL_NaOH,4),'    pH =  ',num2str(pH,4)];disp(T1);
k = k + 1;
a(k) = mL_NaOH;b(k) = pH;
end
disp('----------------------------------')
plot(a,b,'b','LineWidth',2) ;grid on;
xlabel('mL NaOH added');
ylabel('pH')
title('Titration of triprotic acid (H3PO4) with NaOH / Buffer system')






