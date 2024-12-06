
clear; close all; clc; 
% Loading Concentration profile and Mass Concentration from DNS and SolConcentrationProfile 

load('dGsed')
load('d_x3')
load('Kt')

Tsim        = 10;           % Time Simulation in Minutes     
GTOl        = 1e-6;     	% Tolerance
alpha       = 0.1;      	% Alpha 
rati 		= 10;
beta 		= alpha/rati;  	% 10, 5, 2.5

%beta           = 0.01;     	% Beta
fd              = 2.2;      	% Fractal Dimension
NN              = 75;
MaxFloc         = 2500;

% Making a folder for fractal dimension 
Folder_fd = ['fd_',num2str(fd,'%.2f')];

if ~exist(Folder_fd, 'dir')
    mkdir(Folder_fd)
end 

% myCluster = parcluster('local');
% myCluster.NumWorkers = str2double(getenv('SLURM_NTASKS'));
% myCluster.JobStorageLocation = getenv('TMPDIR');
% myPool = parpool(myCluster, myCluster.NumWorkers);

NfolT = dir(Folder_fd);
Nfol = length(NfolT)-2;

if Nfol == 0 
    InFolder = 1;
    load('tav_MassC')  % The first condition 
else
    InFolder = Nfol+1;
    
    for iFol = 1:Nfol
        IterFolder{iFol} = [Folder_fd, '/Iter_' num2str(iFol)];
    end 
    
    load([Folder_fd, '/Iter_' num2str(Nfol), '/tav_MassC'])
    
end 


for ki = InFolder:InFolder+20
    
    coi = 0;
    nIter = ki; 
    disp("nIter = " + nIter)
    
    IterFolder{ki} = [Folder_fd, '/Iter_' num2str(nIter)];
    mkdir(IterFolder{ki})
        
    %Call subroutines
    
    parfor iG = 1:length(dGsed)
        FlocModDNSFun(fd,alpha,beta,IterFolder,nIter,Tsim,tav_MassC(iG),dGsed(iG),NN,MaxFloc)        
    end 
    
        
    WeightedFun(fd,alpha,beta,IterFolder,nIter,dGsed,NN,MaxFloc)
  
    
    [tav_MassC] = MassConSolverFun(IterFolder,nIter,tav_MassC);  

    if nIter > 1
    
        cd(IterFolder{nIter-1})
        load('wav_ws')
        wwsOld = wav_ws;
        clear wav_ws;
        cd('../..')

        cd(IterFolder{nIter})
        load('wav_ws')
        wwsNew = wav_ws;
        clear wav_ws;
        cd('../..')
        

        RMSE = sqrt(mean((wwsOld - wwsNew).^2));
        
        disp("RMSE = " + RMSE)
        
        if RMSE < GTOl
            return 
        end     
    end  
    
end 
    
    

    


