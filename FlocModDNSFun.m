

% fm_flocmod_main keyvani and strom

function FlocModDNSFun(fd,alpha,beta,IterFolder,nIter,Tsim,tav_MassC,idGsed,NN,MaxFloc)
format long

cd(IterFolder{nIter});

time_min            = Tsim;          % Time simulation in minutes     to save = 1000  real = 1800 for ps1
IniConVector_tem    = tav_MassC;

% FLoc Inputs 

bite        = 0;            % ====> 0 = binary, and 0.5 ternary
fra         = 0;
nk          = 2;            % ====> number of partilce I want to erode. By the default
fc          = 1;            % ====> fragment class. Which class do you want to eroded. Values go from class #1 to the highest class (14) by the default is 1
dt_tem      = 10;     

% Floc size distribution 

% NN          = 75;                                               % Number of classes  
% MaxFloc     = 2200;                                             % Maximum Size Floc in micros 
f_diam      = (logspace(log10(4),log10(MaxFloc),NN))'.*1e-6;    % Size Class in m

IniSize = 20*1e-6;                                          % Location for Initial Concentration
[~, LocalIniCon] = min(abs(f_diam-IniSize));                % Search the Location in f_diam


for yy = 1:length(idGsed) 
    for ii = 1:length(alpha)
      for jj = 1:length(beta)

          GG = idGsed(yy);
          IniCon = IniConVector_tem(yy); % Solo coje la primera concentracion
          
% Changing the path to work in the original 

% cd(IterFolder{nIter});

% =======================================================================
% =======================================================================
% =================== Verney's model ====================================
% =======================================================================
% =======================================================================

vis = 1.05e-6; % ====> new variable 
grav = 9.81;
rhoref = 1035; % ====> new variable 
mu = vis*rhoref;
rhosp = 2650;
f_dp0 = 4e-6;
f_nf = fd;


l_ADS=0;        % ======================> differential settling aggregation
l_ASH=1;        % ======================> shear aggregation
l_COLLFRAG=0;   % ======================> collision-induced fragmentation
f_dmax=0.001500;
f_nb_frag = 2;        % number of fragments by shear erosion. If binary/ternary : 2. 
f_alpha=alpha(ii);     % flocculation efficiency, ranging from 0 .to 1. =======  
f_beta=beta(jj);       % shear fragmentation rate =============================  
f_ater = bite;
% f_ero_frac=0.0;          % fraction of the shear fragmentation term transfered to shear erosion. Ranging from 0. (no erosion) to 1. (all erosion)
f_ero_frac= fra;
f_ero_nbfrag= nk;         % Number of fragments induced by shear erosion. EROSION
f_ero_iv=fc;           % fragment size class (could be changed to a particle size or a particle distribution
f_mneg_param=0.000;        % negative mass tolerated to avoid small sub time step (g/l)
f_collfragparam=0.00;      % Fragmentation rate for collision-induced breakup
f_test=1;
dfragmax=0.00003;
epsilon = 1e-8;
% min concentration below which flocculation processes are not calculated
f_clim=0.001;

% dt = 10.0;
dt = dt_tem;
tstart = 0.0;
tend = time_min*60.0; % =============================================================================================> time simulation (min)*(conver to sec)
t_print = 1;   % This value used to be 600 to print 

G_value = GG;
nv_mud = length(f_diam); %c total number of the clases
% ================Initial concentrations ===============================
% TODO - cv_wat should (nv_mud,t)

cv_wat = zeros(nv_mud,1);
cv_wat(LocalIniCon)=IniCon; % this is based on the experiment concentration % =================================> Initial concentration 
% cv_wat(7)=0.5;

% cv_wat(5)=0.093; % UNITS [kg/m3] % Sediments were directly introduced in the
% % test chamber without any treatment at a concentration of
% % 93 mg l  1 , and deflocculated by high turbulence mixing to reach
% % the initial condition of the tidal test.

f_vol = (pi/6.0)*f_diam.^3;
f_area=(pi/4.0)*f_diam.^2;
f_rho = rhoref+(rhosp-rhoref)*(f_dp0./f_diam).^(3.0-f_nf);
f_mass = zeros(nv_mud+1,1);
f_mass(1:nv_mud) = f_vol.*(f_rho-rhoref);
f_mass(nv_mud+1) = f_mass(nv_mud)*2.0+1.0;
if (f_diam(1) == f_dp0)
   f_mass(1)=f_vol(1)*rhosp;
end

f_ws = (grav/(18.*vis))*((f_rho-rhoref)/rhoref).*f_diam.^2.0;
 
% ======================================================================
%fm_print_init

% % fm_print_init - Print out initial values
%   fprintf(1,'\n');
%   fprintf(1,'***********************\n')
%   fprintf(1,'    FLOCMOD\n')
%   fprintf(1,'***********************\n')
%   fprintf(1,'class  diameter (um)  volume (m3)  density (kg/m3)  mass (kg) Ws (m/s)\n')
%   for iv=1:nv_mud
%      fprintf(1,'% 3d % 7.1f % 14g % 6.1f % 14g % 12f\n',...
%         iv,f_diam(iv)*1e6,f_vol(iv),f_rho(iv),f_mass(iv),f_ws(iv))
%   end
%   fprintf(1,'\n')
%   fprintf(1,' *** PARAMETERS ***\n')
%   fprintf(1,'Primary particle size (f_dp0)                                : %f\n',f_dp0)
%   fprintf(1,'Fractal dimension (f_nf)                                     : %f\n',f_nf)
%   fprintf(1,'Flocculation efficiency (f_alpha)                            : %f\n',f_alpha)
%   fprintf(1,'Floc break up parameter (f_beta)                             : %f\n',f_beta)
%   fprintf(1,'Nb of fragments (f_nb_frag)                                  : %f\n',f_nb_frag)
%   fprintf(1,'Ternary fragmentation (f_ater)                               : %f\n',f_ater)
%   fprintf(1,'Floc erosion (pct of mass) (f_ero_frac)                      : %f\n',f_ero_frac)
%   fprintf(1,'Nb of fragments by erosion (f_ero_nbfrag)                    : %f\n',f_ero_nbfrag)
%   fprintf(1,'fragment class (f_ero_iv)                                    : %f\n',f_ero_iv)
%   fprintf(1,'negative mass tolerated before redistribution (f_mneg_param) : %f\n',f_mneg_param)
%   fprintf(1,'Boolean for differential settling aggregation (L_ADS)        : %d\n',l_ADS)
%   fprintf(1,'Boolean for shear aggregation (L_ASH)                        : %d\n',l_ASH)
%   fprintf(1,'Boolean for collision fragmenation (L_COLLFRAG)              : %d\n',l_COLLFRAG)
%   fprintf(1,'Collision fragmentation parameter (f_collfragparam)          : %f\n',f_collfragparam)
%   fprintf(1,'\n')
%   fprintf(1,'*** END FLOCMOD INIT *** \n')    
% 
%   if ~(l_ADS+l_ASH)
%      fprintf(1,'CAUTION : incompatible flocculation kernel options : \n')
%      fprintf(1,'*****************************************************\n')
%      fprintf(1,'l_ADS=%d\n',l_ADS)
%      fprintf(1,'l_ASH=%d\n',l_ASH)
%      error('simulation stopped')
%   end
  

% ======================================================================

% dimension arrays
f_coll_prob_sh=zeros(nv_mud,nv_mud);   % shear agregation collision probability 
f_coll_prob_ds=zeros(nv_mud,nv_mud);   % differential settling collsion probability 
f_g1_sh = zeros(nv_mud,nv_mud,nv_mud); % aggregation due to shear stress
f_g1_ds = zeros(nv_mud,nv_mud,nv_mud); % aggregation due to differential settling velocity 
f_g3 = zeros(nv_mud,nv_mud);           % shear fragmentation  
f_l3 = zeros(nv_mud);                  % fragmentation loss term  
f_g4=zeros(nv_mud,nv_mud,nv_mud);      % Collision fragmentation gain term 
f_l4=zeros(nv_mud,nv_mud);             % Collision fragmentation losst term 

% floc kernals

% fm_flocmod_aggregation_statistics
% flocmod_agregation_statistics

for iv1=1:nv_mud
   for iv2=1:nv_mud
      % Shear (eqn 9)
      f_coll_prob_sh(iv1,iv2)=1.0/6.0*(f_diam(iv1)+f_diam(iv2))^3.0;
      % Differential settling
      f_coll_prob_ds(iv1,iv2)=0.250*pi*(f_diam(iv1)+f_diam(iv2))^2.0 ...
         *grav/mu*abs((f_rho(iv1)-rhoref)*f_diam(iv1)^2.0 ...
         -(f_rho(iv2)-rhoref)*f_diam(iv2)^2.0);      
   end
end

%=========================================================================
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

% fm_aggregation_gain.
% fm_aggregation_gain

%********************************************************************************
% agregation : GAIN : f_g1_sh and f_g1_ds
%********************************************************************************

% the Fortran indexing starting at 0 is hard to replicate in Matlab,
% so I have used diffmass and f_masslo to work around.
diffmass = diff([f_mass; 0]); % getting the difference between adjancent elements of f_mass
for iv1=1:nv_mud
   for iv2=1:nv_mud
      for iv3=iv2:nv_mud
          
         if(iv1==1)
            f_masslo = 0.0;
         else
            f_masslo = f_mass(iv1-1);
         end
         if((f_mass(iv2)+f_mass(iv3)) > f_masslo ...
               && ((f_mass(iv2)+f_mass(iv3)) <= f_mass(iv1)))
            
            %f_weight=(f_mass(iv2)+f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1));
            f_weight=(f_mass(iv2)+f_mass(iv3)-f_masslo)/(diffmass(iv1));
            
         elseif ((f_mass(iv2)+f_mass(iv3)) > f_mass(iv1) ...
               && ((f_mass(iv2)+f_mass(iv3)) < f_mass(iv1+1)))
            
            if (iv1 == nv_mud)
               f_weight=1.0;
            else
               %f_weight=1.0-(f_mass(iv2)+f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1));
               f_weight=1.0-(f_mass(iv2)+f_mass(iv3)-f_mass(iv1))/(diffmass(iv1+1));
            end
            
         else
            f_weight=0.0;
         end
         
         f_g1_sh(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_sh(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1);
         f_g1_ds(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_ds(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1);
         
% f_coll_pro_sh = shear agregation collision probability A(i,j) = (1/6)*G*(D_i+D_j)^3
% f_coll_prob_ds =  differential settling collision probability.  This is neglected
         
                  
      end
   end
end
clear f_weight


% fm_shear_frag_gain
% fm_shear_frag_gain
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Shear fragmentation : GAIN : f_g3
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

for iv1=1:nv_mud
   for iv2=iv1:nv_mud
      
      if(iv1==1)
         f_masslo = 0.0;
      else
         f_masslo = f_mass(iv1-1);
      end
      
      if (f_diam(iv2)>dfragmax)  %--sediment must be higher that 0.00003
         % binary fragmentation   
         if (f_mass(iv2)/f_nb_frag > f_masslo ...f_ero_nbfrag
               && f_mass(iv2)/f_nb_frag <= f_mass(iv1))
            
            if (iv1 == 1)
               f_weight=1.0;
                                       
            else
               f_weight=(f_mass(iv2)/f_nb_frag-f_masslo)/(f_mass(iv1)-f_masslo); % m_i/2
          
            end
         elseif (f_mass(iv2)/f_nb_frag > f_mass(iv1) ...
               && f_mass(iv2)/f_nb_frag < f_mass(iv1+1))
            f_weight=1.0-(f_mass(iv2)/f_nb_frag-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1));
            
         else            
            f_weight=0.0;
            
         end
      else
         f_weight=0.0;
           
      end
      
       f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.0-f_ero_frac)*(1.0-f_ater)*f_weight*f_beta ...
         *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)^(3.0-f_nf)           ...
         *f_mass(iv2)/f_mass(iv1);
      
      % ternary fragmentation
      if (f_diam(iv2)>dfragmax)
         if (f_mass(iv2)/(2.0*f_nb_frag) > f_masslo ...
               && f_mass(iv2)/(2.0*f_nb_frag) <= f_mass(iv1))
            if (iv1 == 1)
               f_weight=1.0;
               
    
            else
               f_weight=(f_mass(iv2)/(2.0*f_nb_frag)-f_mass(iv1-1))/(f_mass(iv1)-f_masslo); % m_i/4
                        

               
            end
         elseif (f_mass(iv2)/(2.0*f_nb_frag) > f_mass(iv1) ...
               && f_mass(iv2)/(2.0*f_nb_frag) < f_mass(iv1+1))
            f_weight=1.0-(f_mass(iv2)/(2.0*f_nb_frag)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1));   
            
         else
            f_weight=0.0;
            
           
         end
         % update for ternary fragments
         

         f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1.0-f_ero_frac)*(f_ater)*f_weight*f_beta ...
            *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)^(3.0-f_nf)           ...
            *f_mass(iv2)/f_mass(iv1);
         
         % Floc erosion
         if ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) > f_mass(f_ero_iv))
            if (((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) >f_masslo) ...
                  && (f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) <= f_mass(iv1))              
               if (iv1 == 1)
                  f_weight=1.0;
                  
                  
               else
                  f_weight=(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag-f_masslo)/(f_mass(iv1)-f_masslo);
                  

               end
            elseif ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) > f_mass(iv1) ...
                  && (f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) < f_mass(iv1+1))
               f_weight=1.0-(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1));
            else
               f_weight=0.0;
            end
            
            % update for eroded floc masses       
            
            
            f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_weight*f_beta                    ...
               *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)^(3.0-f_nf)           ...
               *(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag)/f_mass(iv1);
            
            if (iv1 == f_ero_iv)               
               f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_beta                           ...
                  *f_diam(iv2)*((f_diam(iv2)-f_dp0)/f_dp0)^(3.0-f_nf)           ...                                      ...
                  *f_ero_nbfrag*f_mass(f_ero_iv)/f_mass(iv1);
            end
         end
      end % condition on dfragmax
   end
end
clear f_weight

%===========================================================================================

% fm_aggregation_loss
% fm_aggregation_loss
%**************************************************************************
%   Shear agregation : LOSS : f_l1
%**************************************************************************
for iv1=1:nv_mud
   for iv2=1:nv_mud     
      if(iv2 == iv1)
         mult=2.0;
      else
         mult=1.0;
      end
      f_l1_sh(iv2,iv1)=mult*f_alpha*f_coll_prob_sh(iv2,iv1);
      f_l1_ds(iv2,iv1)=mult*f_alpha*f_coll_prob_ds(iv2,iv1);      
   end
end
clear mult

% ====================================================================

% fm_shear_frag_loss
%**************************************************************************
%  Shear fragmentation : LOSS : f_l2 (stet...f13)
%**************************************************************************
% TODO - this is easy to vectorize
for iv1=1:nv_mud
   if (f_diam(iv1)>dfragmax)
      % shear fragmentation
      f_l3(iv1)=f_l3(iv1)+(1.0-f_ero_frac)*f_beta*f_diam(iv1)*((f_diam(iv1)-f_dp0)/f_dp0)^(3.0-f_nf);
      % shear erosion
      if ((f_mass(iv1)-f_mass(f_ero_iv)*f_ero_nbfrag) > f_mass(f_ero_iv))
         f_l3(iv1)=f_l3(iv1)+f_ero_frac*f_beta*f_diam(iv1)*((f_diam(iv1)-f_dp0)/f_dp0)^(3.0-f_nf);
      end
   end
end

% ==============================================================================================

% fm_kernal_stats
% fm_kernal_stats
knames = {'f_coll_prob_sh',...
'f_coll_prob_ds',...
'f_g1_sh',...
'f_g1_ds',...
'f_g3',...
'f_l3'};

for i=1:length(knames)
% fprintf('%s ',knames{i})
%eval(['size( ',char(knames{i}),')'])
eval(['val(1)=min( ',char(knames{i}),'(:));'])
eval(['val(2)=max( ',char(knames{i}),'(:));'])
eval(['val(3)=sum( ',char(knames{i}),'(:));'])
eval(['val(4)=sum(abs( ',char(knames{i}),'(:)));'])
%fprintf(1,'Min: %g Max: %g Sum: %g SumAbs: %g\n',val)
% fprintf('%g\n',val(3))
end

% =======================================================================================

% I think we can delete fmass(nv_mud+1) now

f_mass = f_mass(1:nv_mud,1);

% fid = fopen(sprintf('data_r%d.dat',ii),'w');
% fid = fopen('fm.dat','w');    % ==================================================================================================> FILE NAME

t = tstart;
nt = 0;

% figure
f_dt=dt;



dt50 = 0;
while (t<tend)
   nt=nt+1;
   dtmin=dt;
      
   dttemp=0.0;
   
   cv_tmp=cv_wat; % concentration of all mud classes in one grid cell
   cvtotmud=sum(cv_tmp);
   
   
   % TODO - fix calculation of G
   % ==================================================================================
   % fm_Gval
   
   % fm_Gval - Returns Gval(t) for test case
% Gval=0.0;

% Gval = sqrt(d_e/vis);  % ======================================================> New Variable 
Gval = G_value;

  
% =====================================================================================
      
   f_gval = Gval;
   
   if( mod(t,t_print) == 0 ) % each 600 seconds is printing resutls
%       fprintf(1,'t, G, cvtotmud: %f %f %f\n',t,Gval,cvtotmud)
   end
   
  NNin=cv_tmp./f_mass; % =======================================================================> Initial concentration 
   
   if( any (NNin<0.0) )
%       fprintf(1,'***************************************\n')
%       fprintf(1,'CAUTION, negative mass at t = %f\n', t)
%       fprintf(1,'***************************************\n')
   end
   
   if (cvtotmud > f_clim)
      while (dttemp <= dt)
         %     print*, 'f_dt:',f_dt
         
         
% ==================================================================         
%          fm_comp_fsd % NNin -> NNout
                  
% fm_comp_fsd
% This processes NNin and returns NNout
NNout = zeros(size(NNin));
tmp_g1=0.0;
tmp_g3=0.0;
tmp_g4=0.0;
tmp_l1=0.0;
tmp_l3=0.0;
tmp_l4=0.0;
f_g1_tmp=zeros(nv_mud,nv_mud,nv_mud);
f_l1_tmp=zeros(nv_mud,nv_mud);

if (l_COLLFRAG)
%    fm_collfrag   % ====> this is another subrouting ====<> collision-induced fragmentation (Paper section 3.2.3)

        % fm_collfrag
f_fp=0.10;
f_fy=1e-10;
f_cfcst=3.0/16.0;

for iv1=1:nv_mud
   for iv2=1:nv_mud
      for iv3=iv2:nv_mud
         % fragmentation after collision probability based on Gval for particles iv2 and iv3
         % gcolfrag=(collision induced shear) / (floc strength)
         gcolfragmin=2.0*(Gval*(f_diam(iv2)+f_diam(iv3))).^2.0*f_mass(iv2)*f_mass(iv3)  ...
            /(pi*f_fy*f_fp*f_diam(iv3).^2.0*(f_mass(iv2)+f_mass(iv3))         ...
            *((f_rho(iv3)-rhoref)/rhoref).^(2.0/(3.0-f_nf)));
         
         gcolfragmax=2.0*(Gval*(f_diam(iv2)+f_diam(iv3))).^2.0*f_mass(iv2)*f_mass(iv3)  ...
            /(pi*f_fy*f_fp*f_diam(iv2).^2.0*(f_mass(iv2)+f_mass(iv3))         ...
            *((f_rho(iv2)-rhoref)/rhoref).^(2.0/(3.0-f_nf)));
         
         % first case : iv3 not eroded, iv2 eroded forming 2 particles : iv3+f_cfcst*iv2 / iv2-f_cfcst*iv2
         if (gcolfragmin<1.0 && gcolfragmax>=10)
            if (((f_mass(iv3)+f_cfcst*f_mass(iv2))>f_mass(iv1-1)) &&  ...
                  ((f_mass(iv3)+f_cfcst*f_mass(iv2))<=f_mass(iv1)))             
               f_weight=((f_mass(iv3)+f_cfcst*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)));
            elseif (f_mass(iv3)+f_cfcst*f_mass(iv2)>f_mass(iv1)  && ...
                  f_mass(iv3)+f_cfcst*f_mass(iv2)<f_mass(iv1+1));
               if (iv1.eq.nv_mud)
                  f_weight=1.0;
               else
                  f_weight=1.0-((f_mass(iv3)+f_cfcst*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)));
               end   
            else
               f_weight=0.0;
            end
            f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   ...
               *(f_mass(iv3)+f_cfcst*f_mass(iv2))/f_mass(iv1);
            
            if (f_mass(iv2)-f_cfcst*f_mass(iv2)>f_mass(iv1-1)   && ...
                  f_mass(iv2)-f_cfcst*f_mass(iv2)<=f_mass(iv1));
               f_weight=((f_mass(iv2)-f_cfcst*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)));
            elseif (f_mass(iv2)-f_cfcst*f_mass(iv2)>f_mass(iv1)  &&  ...
                  f_mass(iv2)-f_cfcst*f_mass(iv2)<f_mass(iv1+1));               
               if (iv1.eq.nv_mud)
                  f_weight=1.0;
               else
                  f_weight=1.0-((f_mass(iv2)-f_cfcst*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)));
               end
            else
               f_weight=0.0;
            end
            f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   ...
               *(f_mass(iv2)-f_cfcst*f_mass(iv2))/f_mass(iv1);   
         % second case : iv3 eroded and iv2 eroded forming 3 particles : iv3-f_cfcst*iv3 / iv2-f_cfcst*iv2 / f_cfcst*iv3+f_cfcst*iv2
         elseif (gcolfragmin>=1.0 && gcolfragmax>=10)   % iv2 and iv3 eroded forming new (third) particle
            if (f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)>f_mass(iv1-1) &&  ...
                  f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)<=f_mass(iv1))
               f_weight=((f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)));
            elseif (f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)>f_mass(iv1) &&  ...
                  f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)<f_mass(iv1+1));
               if (iv1==nv_mud)
                  f_weight=1.0;
               else
                  f_weight=1.0-((f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)));
               end
            else
               f_weight=0.0;
            end
            f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   ...
               *(f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3))/f_mass(iv1);
            if ((1.0-f_cfcst)*f_mass(iv2)>f_mass(iv1-1) &&  ...
                  (1.0-f_cfcst)*f_mass(iv2)<=f_mass(iv1))
               f_weight=((1.0-f_cfcst)*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)); 
            elseif ((1.0-f_cfcst)*f_mass(iv2)>f_mass(iv1) &&  ...
                  (1.0-f_cfcst)*f_mass(iv2)<f_mass(iv1+1))   
               if (iv1==nv_mud)
                  f_weight=1.0;
               else
                  f_weight=1.0-(((1.0-f_cfcst)*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)));
               end            
            else
               f_weight=0.0;
            end           
            f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   ...
               *((1.0-f_cfcst)*f_mass(iv2))/f_mass(iv1);
            if ((1.0-f_cfcst)*f_mass(iv3)>f_mass(iv1-1) &&  ...
                  (1.0-f_cfcst)*f_mass(iv3)<=f_mass(iv1))              
               f_weight=((1.0-f_cfcst)*f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1));
            elseif ((1.0-f_cfcst)*f_mass(iv3)>f_mass(iv1) &&  ...
                  (1.0-f_cfcst)*f_mass(iv3)<f_mass(iv1+1));
               if (iv1.eq.nv_mud)
                  f_weight=1.0;
               else
                  f_weight=1.0-(((1.0-f_cfcst)*f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)));
               end;               
            else
               f_weight=0.0;
            end
            f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   ...
               *((1.0-f_cfcst)*f_mass(iv3))/f_mass(iv1);
         end % end collision test case
      end
   end
end

for iv1=1:nv_mud
   for iv2=1:nv_mud
      
      gcolfragiv1=2.0*(Gval*(f_diam(iv1)+f_diam(iv2))).^2.0*f_mass(iv1)*f_mass(iv2)  ...
         /(pi*f_fy*f_fp*f_diam(iv1).^2.0*(f_mass(iv1)+f_mass(iv2))         ...
         *((f_rho(iv1)-rhoref)/rhoref).^(2.0/(3.0-f_nf)));
      
      gcolfragiv2=2.0*(Gval*(f_diam(iv1)+f_diam(iv2))).^2.0*f_mass(iv1)*f_mass(iv2)  ...
         /(pi*f_fy*f_fp*f_diam(iv2).^2.0*(f_mass(iv1)+f_mass(iv2))         ...
         *((f_rho(iv2)-rhoref)/rhoref).^(2.0/(3.0-f_nf)));
      
      mult=1.0;
      if (iv1.eq.iv2); mult=2.0; end
      if (iv1.eq.MAX(iv1,iv2) && gcolfragiv1>=1.0)
         f_l4(iv2,iv1)=f_l4(iv2,iv1)+mult*(f_coll_prob_sh(iv1,iv2));
      elseif (iv1.eq.MIN(iv1,iv2) && gcolfragiv2>=1.0)
         f_l4(iv2,iv1)=f_l4(iv2,iv1)+mult*(f_coll_prob_sh(iv1,iv2));
      end
   end
end

f_g4(1:nv_mud,1:nv_mud,1:nv_mud)=f_g4(1:nv_mud,1:nv_mud,1:nv_mud)*f_collfragparam;
f_l4(1:nv_mud,1:nv_mud)=f_l4(1:nv_mud,1:nv_mud)*f_collfragparam;

% =============================== END ========================================

end

for iv1=1:nv_mud
%    for iv2=1:nv_mud
%       for iv3=1:nv_mud
         %if (l_ASH)
            f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ASH*f_g1_sh(:,:,iv1)*Gval;
         %end
         %if (l_ADS)
            f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ADS*f_g1_ds(:,:,iv1); % l_ADS = 0 because diferential settling aggregation is neglected
         %end
              
         tmp_g1=tmp_g1+(NNin'*(f_g1_tmp(:,:,iv1))*NNin); % Temporal Aggregation Paper section (3.2.1)
         
         %if (l_COLLFRAG)
            tmp_g4=tmp_g4+l_COLLFRAG*(NNin'*(f_g4(:,:,iv1)*Gval)*NNin); % l_COLLFRAF = 0 because Collision-induced breakup or fragmentations is neglected
         %end
%       end
      
      tmp_g3=tmp_g3+f_g3(:,iv1)'*NNin*Gval^1.5;  % Temporal Shear fragmentation 
      
      %if (l_ASH)
         f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ASH*f_l1_sh(:,iv1)*Gval;
      %end
      %if (l_ADS)
         f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ADS*f_l1_ds(:,iv1)*Gval;
      %end  
      
      
      
      tmp_l1=tmp_l1+(f_l1_tmp(:,iv1))'*NNin;
      
      %if (l_COLLFRAG)
         tmp_l4=tmp_l4+l_COLLFRAG*(f_l4(:,iv1)*Gval)'*NNin; % l_COLLFRAF = 0 because Collision-induced breakup or fragmentations is neglected
      %end
%    end
   
% f_l4 = collision fragmentaion loss term 

   tmp_l1=tmp_l1*NNin(iv1);
   tmp_l4=tmp_l4*NNin(iv1);
   
   tmp_l3=f_l3(iv1)*Gval^1.50*NNin(iv1);
   
   NNout(iv1)=NNin(iv1)+f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+tmp_l3+tmp_l4));
   
   tmp_g1=0.0;
   tmp_g3=0.0;
   tmp_g4=0.0;
   tmp_l1=0.0;
   tmp_l3=0.0;
   tmp_l4=0.0;
end
                  
% ==================================================================         
         
         % fm_mass_control         
         ineg = find(NNout<0.0);
         mneg = sum( -NNout(ineg).*f_mass(ineg) );
         
         %     fprintf(1, 'mneg',mneg
         if (mneg > f_mneg_param)
            while (mneg > f_mneg_param)
               f_dt=min(f_dt/2.0,dt-dttemp);
               
               
%fm_comp_fsd % NNin -> NNout  =======> My Modification 
               
               NNout = zeros(size(NNin));
tmp_g1=0.0;
tmp_g3=0.0;
tmp_g4=0.0;
tmp_l1=0.0;
tmp_l3=0.0;
tmp_l4=0.0;
f_g1_tmp=zeros(nv_mud,nv_mud,nv_mud);
f_l1_tmp=zeros(nv_mud,nv_mud);

if (l_COLLFRAG)
   fm_collfrag
end

for iv1=1:nv_mud
%    for iv2=1:nv_mud
%       for iv3=1:nv_mud
         %if (l_ASH)
            f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ASH*f_g1_sh(:,:,iv1)*Gval;
         %end
         %if (l_ADS)
            f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ADS*f_g1_ds(:,:,iv1);
         %end
         
     
         tmp_g1=tmp_g1+(NNin'*(f_g1_tmp(:,:,iv1))*NNin);
         
         %if (l_COLLFRAG)
            tmp_g4=tmp_g4+l_COLLFRAG*(NNin'*(f_g4(:,:,iv1)*Gval)*NNin);
         %end
%       end
      
      tmp_g3=tmp_g3+f_g3(:,iv1)'*NNin*Gval^1.5;
      
      %if (l_ASH)
         f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ASH*f_l1_sh(:,iv1)*Gval;
      %end
      %if (l_ADS)
         f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ADS*f_l1_ds(:,iv1)*Gval;
      %end
      
      
      tmp_l1=tmp_l1+(f_l1_tmp(:,iv1))'*NNin;
      
      %if (l_COLLFRAG)
         tmp_l4=tmp_l4+l_COLLFRAG*(f_l4(:,iv1)*Gval)'*NNin;
      %end
%    end
   
   tmp_l1=tmp_l1*NNin(iv1);
   tmp_l4=tmp_l4*NNin(iv1);
   
   tmp_l3=f_l3(iv1)*Gval^1.50*NNin(iv1);  

   
   NNout(iv1)=NNin(iv1)+f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+tmp_l3+tmp_l4));
   
   tmp_g1=0.0;
   tmp_g3=0.0;
   tmp_g4=0.0;
   tmp_l1=0.0;
   tmp_l3=0.0;
   tmp_l4=0.0;
end


% =========================================================================
               
               
               
               % fm_mass_control
               ineg = find(NNout<0.0);
               mneg = sum( -NNout(ineg).*f_mass(ineg) );
            end
            %         if (f_dt<1.0)
            %           fprintf(1, 'apres : Gval,f_dt',Gval, f_dt,dttemp
            %  end
         else
            
            if (f_dt<dt)
               while (mneg <=f_mneg_param)
                  
                  if (dttemp+f_dt == dt)
 %fm_comp_fsd % NNin -> NNout  ==========> My Modification 
 
 % fm_comp_fsd
% This processes NNin and returns NNout
NNout = zeros(size(NNin));
tmp_g1=0.0;
tmp_g3=0.0;
tmp_g4=0.0;
tmp_l1=0.0;
tmp_l3=0.0;
tmp_l4=0.0;
f_g1_tmp=zeros(nv_mud,nv_mud,nv_mud);
f_l1_tmp=zeros(nv_mud,nv_mud);

if (l_COLLFRAG)
   fm_collfrag
end

for iv1=1:nv_mud
%    for iv2=1:nv_mud
%       for iv3=1:nv_mud
         %if (l_ASH)
            f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ASH*f_g1_sh(:,:,iv1)*Gval;
         %end
         %if (l_ADS)
            f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ADS*f_g1_ds(:,:,iv1);
         %end
         
     
         tmp_g1=tmp_g1+(NNin'*(f_g1_tmp(:,:,iv1))*NNin);
         
         %if (l_COLLFRAG)
            tmp_g4=tmp_g4+l_COLLFRAG*(NNin'*(f_g4(:,:,iv1)*Gval)*NNin);
         %end
%       end
      
      tmp_g3=tmp_g3+f_g3(:,iv1)'*NNin*Gval^1.5;
      
      %if (l_ASH)
         f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ASH*f_l1_sh(:,iv1)*Gval;
      %end
      %if (l_ADS)
         f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ADS*f_l1_ds(:,iv1)*Gval;
      %end
      
      
      tmp_l1=tmp_l1+(f_l1_tmp(:,iv1))'*NNin;
      
      %if (l_COLLFRAG)
         tmp_l4=tmp_l4+l_COLLFRAG*(f_l4(:,iv1)*Gval)'*NNin;
      %end
%    end
   
   tmp_l1=tmp_l1*NNin(iv1);
   tmp_l4=tmp_l4*NNin(iv1);
   
   tmp_l3=f_l3(iv1)*Gval^1.50*NNin(iv1);  

   
   NNout(iv1)=NNin(iv1)+f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+tmp_l3+tmp_l4));
   
   tmp_g1=0.0;
   tmp_g3=0.0;
   tmp_g4=0.0;
   tmp_l1=0.0;
   tmp_l3=0.0;
   tmp_l4=0.0;
end


% ===================================================================================
 
 
 
                     break
                  else
                     dt1=f_dt;
                     f_dt=min(2.0*f_dt,dt-dttemp);
% fm_comp_fsd % NNin -> NNout  ===================> My modification 

% fm_comp_fsd
% This processes NNin and returns NNout
NNout = zeros(size(NNin));
tmp_g1=0.0;
tmp_g3=0.0;
tmp_g4=0.0;
tmp_l1=0.0;
tmp_l3=0.0;
tmp_l4=0.0;
f_g1_tmp=zeros(nv_mud,nv_mud,nv_mud);
f_l1_tmp=zeros(nv_mud,nv_mud);

if (l_COLLFRAG)
   fm_collfrag
end

for iv1=1:nv_mud
%    for iv2=1:nv_mud
%       for iv3=1:nv_mud
         %if (l_ASH)
            f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ASH*f_g1_sh(:,:,iv1)*Gval;
         %end
         %if (l_ADS)
            f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ADS*f_g1_ds(:,:,iv1);
         %end
         
     
         tmp_g1=tmp_g1+(NNin'*(f_g1_tmp(:,:,iv1))*NNin);
         
         %if (l_COLLFRAG)
            tmp_g4=tmp_g4+l_COLLFRAG*(NNin'*(f_g4(:,:,iv1)*Gval)*NNin);
         %end
%       end
      
      tmp_g3=tmp_g3+f_g3(:,iv1)'*NNin*Gval^1.5;
      
      %if (l_ASH)
         f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ASH*f_l1_sh(:,iv1)*Gval;
      %end
      %if (l_ADS)
         f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ADS*f_l1_ds(:,iv1)*Gval;
      %end
      
      
      tmp_l1=tmp_l1+(f_l1_tmp(:,iv1))'*NNin;
      
      %if (l_COLLFRAG)
         tmp_l4=tmp_l4+l_COLLFRAG*(f_l4(:,iv1)*Gval)'*NNin;
      %end
%    end
   
   tmp_l1=tmp_l1*NNin(iv1);
   tmp_l4=tmp_l4*NNin(iv1);
   
   tmp_l3=f_l3(iv1)*Gval^1.50*NNin(iv1);  

   
   NNout(iv1)=NNin(iv1)+f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+tmp_l3+tmp_l4));
   
   tmp_g1=0.0;
   tmp_g3=0.0;
   tmp_g4=0.0;
   tmp_l1=0.0;
   tmp_l3=0.0;
   tmp_l4=0.0;
end
                    

% ========================================================                     
                     
                     
                     % fm_mass_control
                     ineg = find(NNout<0.0);
                     mneg = sum( -NNout(ineg).*f_mass(ineg) );
                     if (mneg > f_mneg_param)
                        f_dt=dt1;
%  fm_comp_fsd % NNin -> NNout ========> My modification 

% fm_comp_fsd
% This processes NNin and returns NNout
NNout = zeros(size(NNin));
tmp_g1=0.0;
tmp_g3=0.0;
tmp_g4=0.0;
tmp_l1=0.0;
tmp_l3=0.0;
tmp_l4=0.0;
f_g1_tmp=zeros(nv_mud,nv_mud,nv_mud);
f_l1_tmp=zeros(nv_mud,nv_mud);

if (l_COLLFRAG)
   fm_collfrag
end

for iv1=1:nv_mud
%    for iv2=1:nv_mud
%       for iv3=1:nv_mud
         %if (l_ASH)
            f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ASH*f_g1_sh(:,:,iv1)*Gval;
         %end
         %if (l_ADS)
            f_g1_tmp(:,:,iv1)=f_g1_tmp(:,:,iv1)+l_ADS*f_g1_ds(:,:,iv1);
         %end
         
     
         tmp_g1=tmp_g1+(NNin'*(f_g1_tmp(:,:,iv1))*NNin);
         
         %if (l_COLLFRAG)
            tmp_g4=tmp_g4+l_COLLFRAG*(NNin'*(f_g4(:,:,iv1)*Gval)*NNin);
         %end
%       end
      
      tmp_g3=tmp_g3+f_g3(:,iv1)'*NNin*Gval^1.5;
      
      %if (l_ASH)
         f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ASH*f_l1_sh(:,iv1)*Gval;
      %end
      %if (l_ADS)
         f_l1_tmp(:,iv1)=f_l1_tmp(:,iv1)+l_ADS*f_l1_ds(:,iv1)*Gval;
      %end
      
      
      tmp_l1=tmp_l1+(f_l1_tmp(:,iv1))'*NNin;
      
      %if (l_COLLFRAG)
         tmp_l4=tmp_l4+l_COLLFRAG*(f_l4(:,iv1)*Gval)'*NNin;
      %end
%    end
   
   tmp_l1=tmp_l1*NNin(iv1);
   tmp_l4=tmp_l4*NNin(iv1);
   
   tmp_l3=f_l3(iv1)*Gval^1.50*NNin(iv1);  

   
   NNout(iv1)=NNin(iv1)+f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+tmp_l3+tmp_l4));
   
   tmp_g1=0.0;
   tmp_g3=0.0;
   tmp_g4=0.0;
   tmp_l1=0.0;
   tmp_l3=0.0;
   tmp_l4=0.0;
end

% ============================================================================================



                        break
                     end
                  end
               end
            end
         end
         dtmin = min(dtmin,f_dt);
         dttemp = dttemp+f_dt;
         NNin = NNout; % update new Floc size distribution
         
		 
% ======================================================================
		 
         %fm_mass_distribute % redistribute negative masses in NNin (if any) over positive classes,
         % depends on f_mneg_param
		 
		 % fm_mass_redistribute -  based on a tolerated negative mass parameter, negative masses
% are redistributed linearly towards remaining postive masses

% TODO - Vectorize this
mneg=0.0;
npos=0.0;
NNtmp=NNin;

for iv=1:nv_mud
   if (NNin(iv) < 0.0) 
      mneg=mneg-NNin(iv)*f_mass(iv);
      NNtmp(iv)=0.0;
   else
      npos=npos+1.0;
   end
end

if (mneg > 0.0) 
   if (npos == 0.0) 
      error('All floc sizes have negative mass!')
   else
      for iv=1:nv_mud
         if (NNin(iv) > 0.00) 
            NNin(iv)=NNin(iv)-mneg/sum(NNtmp)*NNin(iv)/f_mass(iv);
         else
            NNin(iv)=0.0;
         end         
      end
   end
end % and negative masses are set to 0
		 
		 
		 % ======================================================================
		 
         
         if (abs(sum(NNin.*f_mass)-cvtotmud)>epsilon*100.0)
%             fprintf(1, 'CAUTION flocculation routine not conservative!\n')
%             fprintf(1, 'time = %f\n',t);
%             fprintf(1, 'f_dt= %f\n',f_dt);
%             fprintf(1, 'before : cvtotmud= %f\n',cvtotmud);
%             fprintf(1, 'after  : cvtotmud= %f\n',sum(NNin.*f_mass));
%             fprintf(1, 'absolute difference  : cvtotmud= %f\n',abs(cvtotmud-sum(NNin.*f_mass)));
%             fprintf(1, 'absolute difference reference  : espilon= %f\n',epsilon);
%             fprintf(1, 'before redistribution %f\n', sum(NNout.*f_mass));
%             fprintf(1, 'after redistribution %f\n', sum(NNin.*f_mass));
            error('Simultation stopped')
         end
         
         if (dttemp == dt); break; end
      end % loop on full dt
   end % only if cvtotmud > f_clim
   
   if (abs( sum( NNin.*f_mass )-cvtotmud) > epsilon*10.0)
%       fprintf(1, 'CAUTION flocculation routine not conservative!\n');
%       fprintf(1, 'time = %g\n',t);
%       fprintf(1, 'before : cvtotmud= %f\n',cvtotmud)
%       fprintf(1, 'after  : cvtotmud= %f\n',sum( NNin.*f_mass ))
%       fprintf(1, 'absolute difference  : cvtotmud= %f\n',...
%          abs(cvtotmud-sum( NNin.*f_mass )))
%       fprintf(1, 'absolute difference reference  : espilon= %f\n',epsilon);
      error('Simultation stopped')
   end
   
   % update mass concentration for all mud classes
   cv_wat  = NNin.*f_mass; % maybe this is the full distribuition concentration
   % Units of NNin is 1/m3 = m-3
      
  
   % compute floc distribution statistics before output
   f_csum=0.0;
   f_ld50=1;
   f_ld16=1;
   f_ld84=1;
   
   f_davg = sum(NNin.*f_mass.*f_diam)./(sum(NNin.*f_mass)+eps);   
   f_dtmin = dtmin;     
   
      
  for iv1=1:nv_mud
      f_csum= f_csum + NNin(iv1)*f_mass(iv1)/((sum(NNin.*f_mass))+eps);
      matrix_csum(iv1,nt) = f_csum;   % Temporal evolution of accumulative sum for each class       
      if (f_csum > 0.16 && f_ld16)
         
         f_d16 =f_diam(iv1);
         f_ld16 = 0;
      end
      
      if (f_csum > 0.5 && f_ld50)
         f_d50 = f_diam(iv1);         
         NNin_50 = NNin(iv1);
         f_mass_50 = f_mass(iv1);
         f_area_50 = f_area(iv1); 
         f_ld50=0;
      end
      
      if (f_csum > 0.84 && f_ld84)
         f_d84=f_diam(iv1);         
         f_ld84=0;      
      end
  end          
      
   
%    fprintf(fid,'%f %f %f %f %f\n',...
%       t, f_dt, Gval, f_d50*1e6, f_d_area_weighted*1e6);

   f_d_area_weighted = (NNin.*f_area)'*f_diam/(sum(NNin.*f_area));
   

   
   Nfs(:,nt) = NNin;    % Number of particles for each band m^-3
   MassConFlocs(:,nt) = cv_wat;  % Mass concentration 
      
   d_area(nt) = f_d_area_weighted*1e6;
   
   d_16(nt) = f_d16*1e6;
   d_50(nt) = f_d50*1e6;
   d_84(nt) = f_d84*1e6; 
   
   dt50(nt+1) = d_50(nt);  
     
       if t > 5*60              
         
            dif = d_50(nt) -  dt50(nt); disp("dif = " + dif)                               
            
            if dif<1e-6 && t>(5+3)*60          
 
               break
               
            end 
            
       end 
   
       tiempo(nt) = t; 
       t = t+dt ;
end

newFolder = ['GG_',num2str(idGsed(yy)),'alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),...
            '_fd_',num2str(fd)];

mkdir(newFolder)
cd(newFolder);


save(['Nfs_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'Nfs')
save(['MassConFlocs_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'MassConFlocs')
save(['d_area_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'d_area')
save(['tiempo_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'tiempo')

save(['d_50_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'d_50')
save(['d_84_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'d_84')
save(['d_16_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'d_16')

save(['matrix_csum_alpha_',num2str(alpha(ii)),'_beta_',num2str(beta(jj)),'_fd_',num2str(fd),'.mat'],'matrix_csum')

save('f_diam.mat','f_diam')     % Discrete diameter 
save('dt.mat','dt')             % time step 
save('tstart.mat','tstart')     % Initial time
save('tend.mat','tend')         % Final time 

cd('..')
      end 
    end 
end 

cd('../../');

end 


