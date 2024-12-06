
function [tav_MassC] = MassConSolverFun(IterFolder,nIter,tav_MassC)

% MassConSolver


load('d_x3')
load('Kt')

cd(IterFolder{nIter})

load('wav_ws')



Ks          = 1e-6;     % The same than the molecular viscosity 
d_ws        = wav_ws;   % Dimensional settling velocity 
Tol         = 1e-6;     % Tolerance 

d_dx3       = gradient(d_x3);
IntMassC    = sum(tav_MassC.*d_dx3); % Integrating the Mass Concentration 

GuessMassC  = 0; 

MassCn      = ones(length(d_x3),1);
MassCn(1)   = GuessMassC;
IntMassCn   = sum(MassCn.*d_dx3);
Error       = abs(IntMassCn - IntMassC);

IterSol = 0;
while Error > Tol
    
    for i = 1:length(d_x3)-1
        MassCn(i+1) = MassCn(i) -( d_dx3(i) * d_ws(i) * MassCn(i) / (Kt(i) + Ks));
    end 
    
    if Error > Tol 
    IntMassCn       = sum(MassCn.*d_dx3);
    Error           = abs(IntMassCn - IntMassC);
    MassCn(1)       = MassCn(1) + 0.0001; 
    IterSol         = IterSol + 1;
%     disp(IterSol);
    end 
    
end 

% Checking the numerical solution using a constant settling velocity 

% close all
% plot(MassCn,d_x3*100,'-b','LineWidth',1.5)
% hold on 
% plot(tav_MassC,d_x3*100,'--r','LineWidth',1.5)

% Saving 

tav_MassC = MassCn;
save('tav_MassC','tav_MassC')
cd('../../');

end 


