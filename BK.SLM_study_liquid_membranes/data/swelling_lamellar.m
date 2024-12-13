% Kilian Frank 2018-09-10

% changes:
% - replaced eval statements by function calls
% - fixed some typos
% - added a legend to the parameter plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        swelling_lamellar.m                              %
%                Sebastian Aeffner, March 23rd, 2010                      % 
%                                                                         %
%             tool to obtain a reasonable phase setting                   %
%                on the basis of the swelling method                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               user input                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear;
%close all;

% sample name and phase, filenames for data input and output
sample = 'DOPC';
phase  = 'L';
file_input  = '09-Jun-2022_DOPC_reflectivity_L.mat' %'27-Apr-2017_DOPC_reflectivity_L.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              main program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
file_output = strcat([date,'_',sample,'_swellingdata_reflectivity_L']);

load(file_input);
disp(['Reading ' file_input ' ... done']);

% display RH values corresponding to the available form factor sets and 
% define which ones will be used to apply the swelling method
disp('Available d-spacings in the dataset:');  
disp(' ');
strarray = cell(1,4);            
strarray{1} = 'index';
strarray{2} = 'RH';
strarray{3} = 'd';
strarray{4} = 'std(d)';
disp(strarray);
disp([linspace(1,length(d_obs(:,1)),length(d_obs(:,1)))',RH, d_obs]);
firstindex = input('Enter index of first form factor set: ');
lastindex =  input('Enter index of last  form factor set: ');

% use only the data we have just defined 
F_obs = F_obs(firstindex:1:lastindex,:);
d_errors = d_obs(firstindex:1:lastindex,2);
d_obs = d_obs(firstindex:1:lastindex,1);
RH = RH(firstindex:1:lastindex);
                             
n_max  = length(F_obs(1,:));       % number of recorded diffraction orders
n_obs  = 1:1:n_max;                % observed Bragg orders 
num_RH = length(F_obs(:,1));       % number of RH values
d_mean = mean(d_obs);              % mean d-spacing of the dataset
F_mean = sum(F_obs)/num_RH;        % mean form factors for each order n

% generate all 2^(n_max) possible phase combinations
% (the phase for F(-)(0) will be determined from the other phases!)
nu = [1,-1]';
for i = 1:n_max-1
    nu = [ [nu ones(length(nu(:,1)),1)] ; -[nu ones(length(nu(:,1)),1)] ];
end
num_phases = length(nu(:,1));                                  

% check if no phase combination occurs twice
for i = 1:num_phases
    for j = (i+1):num_phases
        if nu(i,:) == nu(j,:)
            disp('Error: ');
        end
    end
end

% generate all 2^{n_max} phased form factors including F(-)(0)
% For L phase: Use formula of Nagle, F^(-)(0) = 2*sum_n (-1)^(n+1)*...
F_mean_phased = zeros(num_phases,n_max+1);
for i = 1:num_phases
    F_mean_phased(i,1) = 2*sum( cos(pi*(1:1:n_max))*(-1).*nu(i,:).*F_mean);
    F_mean_phased(i,(2:1:n_max+1)) = nu(i,:).*F_mean;
end
   
% due to centrosymmetry: incorporate "negative" orders of diffraction
n_obs = [-fliplr(n_obs) 0 n_obs];  
F_mean_phased = [fliplr(F_mean_phased(:,2:1:(n_max+1))) F_mean_phased];
F_obs   = [fliplr(F_obs) zeros(num_RH,1) F_obs];   
q_z_obs = (1./d_obs)*2*pi*n_obs;


% data required for interpolation of continuous form factor 

% generate q_z axis for form factor plot
delta_q_z = 50;                    
n_q = 2*(n_max+1)*delta_q_z;
q_z = linspace(-(n_max+1)*2*pi/d_mean,+(n_max+1)*2*pi/d_mean,n_q);

% sinc functions centered around q_z = n*2*pi/d_mean, n=-n_max...n_max
sinc_mean = zeros((2*n_max)+1,n_q); 
for i = 1:length(n_obs)
    argument = d_mean/2*q_z - n_obs(i)*pi;
    sinc_mean(i,:) = sin(argument)./(argument);  
end

% generate all possible mean continuous form factors for d_mean, use en-
% tries in F_mean as coefficients. 
F_mean_cont = F_mean_phased * sinc_mean;

% determine sum of squared errors 
% vector containing sum of squared residuals for each phase setting

errors = zeros(num_phases,1);

for i = 1:num_phases   % for all possible phase combinations    
    sum_squares = 0;  
    F_obs_phased = (ones(num_RH,1)*[fliplr(nu(i,:)) 0 nu(i,:)]).*F_obs;
    
    for j = 1:num_RH           % for all RH_num humidity value   
        for k = 2:n_max+1      % for all n_max measured datapoints
        % calculate function value of the mean, continuous form
        % factor at the observed q_z value and subtract the
        % observed F(q_z) value (-> NO F(-)(0))
        argument = 0.5*d_mean*q_z_obs(j,n_max+k)-n_obs*pi;
        delta = F_mean_phased(i,:)*(sin(argument)./(argument))' - ...
                                                   F_obs_phased(j,n_max+k);
        sum_squares = sum_squares + delta^2;  
        end
    end
    errors(i) = sum_squares;  
end
    
% sort phase combinations based on their fit quality 
[errors_sorted, ind] = sort(errors,1);    
nu_sorted = [];                              
for i = 1:length(nu(:,1))
    nu_sorted = [nu_sorted; nu(ind(i),:)];
end

% indices, phases and phased observed form factors for the two best phase
% combinations related by total phase factor (-1) 
i_opt = find(errors==min(errors));
nu_opt = nu(i_opt,:);
%if strcmp(phase,'L') == 1;
    if nu_opt(1,1) == -1
        i_opt  = i_opt(1);
        nu_opt = nu_opt(1,:);
    else i_opt  = i_opt(2);
        nu_opt = nu_opt(2,:);
    end
%end;
F_phased_opt = (ones(num_RH,1)*[fliplr(nu_opt) 1 nu_opt]).*F_obs;


% number of swelling plots which will be shown
num_plots = input('Enter number of swelling plots you would like to see: ');

% Plot swelling graph for the num_plots best independent phase combinations
for h = 1:2:(2*num_plots-1)
    
    % out of 2 equivalent phase settings differing by global factor (-1),
    % use the one which leads to EDP with bilayer interior centered around
    % origin 
    if nu_sorted(h,1) ~= -1
        h = h + 1;
    end
    
    
    F_obs_phased  = ...
         (ones(num_RH,1)*[fliplr(nu_sorted(h,:)) 1 nu_sorted(h,:)]).*F_obs;
     
    % create 'linear' versions of observed data required for plot
    q_z_obs_lin = [];                
    F_obs_phased_lin = [];        
    for i = 1:num_RH
        q_z_obs_lin = [q_z_obs_lin q_z_obs(i,:)];
        F_obs_phased_lin = [F_obs_phased_lin F_obs_phased(i,:)];
    end

    % plot swelling graph
    F_continuous = F_mean_cont(ind(h),:);
    figure;
    hold on;
    plot(q_z_obs_lin,F_obs_phased_lin,'o','color','black');
    plot(q_z,F_continuous,'-','color','black');
    plot(q_z,zeros(length(q_z)),'--','color','black');
    plot(0,F_mean_phased(ind(h),n_max+1),'s','MarkerEdgeColor','k',...
                                     'MarkerFaceColor','k','MarkerSize',8);
    xlim([0 max(q_z)]);
    xlabel('$q_z$ (\AA$^{-1}$)','interpreter','latex','fontsize',16);
%     ylim([min(F_continuous)-0.1 max(F_continuous)+0.2]);
    ylim([-2 0.3]);
    ylabel('$F^{(-)}(q_z)$ (a.u.)','interpreter','latex','fontsize',16);
    title(['Swelling plot for ' sample],'fontsize',16);
    text(0.3,0.5,['phase ranking: ' num2str(ceil(h/2)) ...
                                 '/' num2str(num_phases/2)],'FontSize',14);

    % add lines showing the residuals and sum um their squares
    %F_mean_phased = [ nu_opt 1 fliplr(nu_opt) ].*F_mean;
    error_min = 0;
    for i = 1:num_RH
        for j = (n_max+2):(2*n_max+1)
            argument = d_mean/2*q_z_obs(i,j)-n_obs*pi;
            xpos_1 = q_z_obs(i,j);
            xpos_2 = q_z_obs(i,j);
            ypos_1 = F_obs_phased(i,j);
            ypos_2 = F_mean_phased(ind(h),:)*(sin(argument)./argument)';
            line([xpos_1 xpos_2],[ypos_1 ypos_2],'linestyle','-');                            
            error_min = error_min + ( ypos_1 - ypos_2 ).^2;
        end
    end
    
    % show error message if sum of resuduals as calculated here and first entry 
    % of errors_sorted are not the same
    if error_min ~= errors_sorted(h)
        disp('Mismatch of errors: Check program!');
    end
    
end

% plot all num_RH resulting electron density profiles for the best phase
% setting, determine approx. values for d_hh and d_w for reach RH
figure;
hold on;
for i = 1:num_RH
    z(i,:)    = linspace(-0.5*d_obs(i),0.5*d_obs(i),1000);
    rho(i,:)  = sum((F_phased_opt(i,:)'*ones(1,1000))...
                                       .*cos(n_obs'*2*pi/d_obs(i)*z(i,:)));
    d_hh(i,:) = 2*mean(abs(z(i,rho(i,:)==max(rho(i,:)))));
  
    plot (z(i,:),rho(i,:),'-');
    d_w(i,:)  = d_obs(i) - d_hh(i);
    % add offset 
    rho(i,:) = rho(i,:)%- 0.5*i; % offset der Kurven in der Darstellung
    plot (z(i,:),rho(i,:),'-');
end                    
title('Most reasonable EDPs for all hydrations','fontsize',16);
xlabel('$z$ (\AA)','interpreter','latex','fontsize',16);
ylabel('$\Delta\rho(z,RH)$ (a.u.), offset','interpreter','latex','fontsize',16);

% plot bilayer parameters
figure;
hold on;
errorbar(RH,d_obs,d_errors,'-o');
plot(RH,d_hh,'-o');
plot(RH,d_w,'-o');
legend('d','d_{HH}','d_w')
xlabel('$RH$ (\%)','interpreter','latex','fontsize',16);
ylabel('$d_w,d_{HH},d$ (\AA)','interpreter','latex','fontsize',16);
title(['bilayer parameters of ' sample],'fontsize',16);


% data output
eval(['save ' file_output]);
disp(['Data saved into ' file_output '.mat']);
disp('Finished!!!');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
