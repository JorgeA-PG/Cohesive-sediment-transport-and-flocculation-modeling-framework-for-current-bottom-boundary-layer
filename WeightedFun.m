
% Calculate the settling velocity using a weigth averaged from the
% concentration in the equillebrium state

function WeightedFun(fd,alpha,beta,IterFolder,nIter,dGsed,NN,MaxFloc)

% load('dGsed','dGsed')
% load('MassC','MassC')

grav        = 9.81;     % Gravity 
rhoref      = 1035;     % Density 
rhosp       = 2650;     % Density for a solid grain of sediment
f_dp0       = 4e-6;     % Primary particle diameter
vis         = 1e-6;

% NN          = 75;                                                 % Number of classes  
% MaxFloc     = 2500;                                               % Maximum Size Floc in micros 
f_diam      = (logspace(log10(4),log10(MaxFloc),NN))'.*1e-6;        % Size Class in m
% fd          = 2.4;                                                % Fractal dimension 

f_rho = rhoref+(rhosp-rhoref)*(f_dp0./f_diam).^(3.0-fd);            % Floc density distribution 
f_ws = (grav/(18.*vis))*((f_rho-rhoref)/rhoref).*f_diam.^2.0;       % Settilng velocity distribution

% Limiting the settling velocity -----------------------------------
%Lws = 5;
%f_ws(f_ws>Lws/1000) = Lws/1000;

% Stuber Formulation ------------------------------------------------
%f_ws = (0.0002.*(f_diam*1e6).^1.54)./1000;                        % Divided by 1000 to have m/s;


cd(IterFolder{nIter})

for yy = 1:length(dGsed) 
for ii = 1:length(alpha)
for jj = 1:length(beta)   
    

    newFolder = ['GG_',num2str(dGsed(yy)),'alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),...
            '_fd_',num2str(fd)];

    cd(newFolder);
    
    load(['MassConFlocs_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'MassConFlocs')    
    load(['tiempo_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'tiempo')
    load(['d_50_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'d_50')
    load(['d_84_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'d_84')
    load(['d_16_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'d_16')
    
    MassConFlocsEq(:,yy)    = MassConFlocs(:,end);
    MassConFlocsTotal       = sum(MassConFlocs(:,end));
    sde                     = MassConFlocs(:,end);
    wei                     = sde./MassConFlocsTotal;
    
    wav_ws(yy)      = sum(wei.*f_ws);   % Weighted Average settling velocity
    wav_d(yy)       = sum(wei.*f_diam); % Weighted Average floc diamenter
    wav_fdens(yy)   = sum(wei.*f_rho);  % Weighted Average floc density
    
    %d16(:,yy) = d_16;
    %d50(:,yy) = d_50;
    %d84(:,yy) = d_84;
    %tfloc     = tiempo;        % No need to storage because it is equal for each G case 
    
    	d16{yy} 	= d_16;
    	d50{yy} 	= d_50;
	    d84{yy} 	= d_84;
    	tfloc{yy}   = tiempo;  

    cd ..
    
end
end
end

save('MassConFlocsEq','MassConFlocsEq')
save('wav_ws','wav_ws')
save('wav_d','wav_d')
save('wav_fdens','wav_fdens')
save('d16','d16')
save('d50','d50')
save('d84','d84')
save('tfloc','tfloc')

cd('../..')
end
