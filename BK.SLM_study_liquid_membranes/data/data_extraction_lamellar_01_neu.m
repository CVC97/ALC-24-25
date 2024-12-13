

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        one-click-version extraction.m for Sebastian and Yihui           %
%                       complete scans only                               %
%              use with UNSPEC, options: -w 1,6 -1 -B -#                  %
%                change filename to .txt after unspec                     %
%          pay attention that _ shall not exist in the file name          %           
%                Sebastian Aeffner, January 11, 2011                      % 
%                 revised by Yihui, March 25, 2013                        %
%                                                                         %
%    tool to obtain form factors and d-spacings from reflectivity scans   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% KF 23.09.2018
% - use the actual values of th from the data
% - deal with slightly inaccurate th beetween datasets
% - avoid running out of the range of the dataset when finding peaks
% - fixed some typos
% - removed the traces of some incomplete functions


%new.m version
%changed: 1. merge unspeced datafiles into one; 2. using [0:0.01:10]
%instead of real values; 3. more notes 

%one-click version
%changed: automatically search for the rest peaks after the first one is clicked out manually 
%Attention: for ugly plots, its better to use the old version

close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               user input                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sample name and phase, filenames for data input and output
sample = 'DOPC';
phase  = 'L';
file_input  = 'messung.txt';
% RH = [85:-1:65 40 30]';
%RH = [90 80 70 60 52 48 46 44]
%RH = [89:-1:59]';
RH = [90:-2:60]';

% file_input  = 'reflectivities_2017_04_26.txt';
% RH = [90:-1:65]';

lambda      = 1.54055;      % used x-ray wavelength (Angstrom)
beamwidth   = 0.5;           % beamwidth: Use value of s1hg (mm)
waferlength = 25;            % sample length in beam direction (mm)
d_film      = 0.000007;      % approx. thickness of lipid film: 7 microns
abs_coeff   = 783;           % absorption coefficient of lipid film (m^-1)
peakwidth   = 40;            % half width of interval for peak search (points)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           main program                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scrsz = get(0,'ScreenSize'); % screen size
file_output = strcat([date,'_',sample,'_reflectivity_',phase]);

% refl = load(file_input);
A = importdata(file_input);
% raw_Detector = A(:,32);
% raw_Icorr = A(:,6);
% raw_th = A(:,1);
%raw_humidity = A.data(:,19);
%plot(1:length(raw_humidity), raw_humidity)

refl = A(:,[1 4]);

disp(['Reading experimental data: ' file_input ' ... done']);
num_RH = length(RH);
num_pts = uint16(length(refl)/num_RH); % we assume that all scans have the same length

th = reshape(refl(:,1),num_pts,num_RH)'; % theta angles
Icorr = reshape(refl(:,2),num_pts,num_RH)'; % reflectivities (one RH is one column)
tth = 2*th; % assume specular reflectivity
q_z = 4*pi/lambda*sind(th);
clear refl;

clear A raw_Detector raw_th raw_Icorr;

% remove entries for q_z = 0 
% (They become "Inf" after corrections and would cause problems)
th(:,1)    = [];
tth(:,1)   = [];
q_z(:,1)   = [];
Icorr(:,1) = [];

num_pts = num_pts-1;
%I_0(find(I_0 == 1)) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     apply correction factors                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% apply illumination correction
q_z_c   = 4*pi/lambda*sin(beamwidth/waferlength); % critical angle    
for i = 1:num_RH
    j = 1;
    while q_z(i,j) < q_z_c
        if Icorr(i,j) > 1
            Icorr(i,j) = Icorr(i,j)*q_z_c/q_z(i,j);  
        end
        j = j + 1;
    end
end
disp('Applying illumination correction ... done');

% apply polarization correction (sealed tube -> unpolarized beam):
C_pol = 0.5*(1+(cosd(tth)).^2);
for i = 1:num_RH
    for j = 1:num_pts
        if Icorr(i,j) > 1
            Icorr(i,j) = Icorr(i,j)./C_pol(j);
        end
    end
end
disp('Applying polarization correction ... done');

% apply absorption correction
C_abs = (1-exp(-1*2*d_film*abs_coeff./sind(0.5*tth)))./(2*d_film*abs_coeff./sind(0.5*tth));
for i = 1:num_RH
    for j = 1:num_pts
        if Icorr(i,j) > 1
            Icorr(i,j) = Icorr(i,j)./C_abs(j);
        end
    end
end
disp('Applying absorption correction ... done');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  plot corrected reflectivities                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Plotting reflectivity curves ... done');
for i = 1:num_RH
    semilogy(q_z(i,:),Icorr(i,:),'-k');
    hold on;
end
% xlim([0 max(q_z)]);
xlabel('$q_z (Angstrom^{-1})$','interpreter','latex','fontsize',16);
ylabel('corrected intensity (cps)','interpreter','latex','fontsize',16);
title([sample ', RH=' num2str(min(RH)) '-' num2str(max(RH)) '%'],'fontsize',16);

dummy = input('Check reflectivity curves and press return to continue');
num_peaks = input('Enter number of observable diffraction orders (all you see, the first and last will be removed later anyway): ');
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       determine integrated peak intensities and d-spacings              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Integrating peak intensities...');
q_z_sampled  = zeros(num_RH,num_peaks);
peakareas    = zeros(num_RH,num_peaks);

for i=1:num_RH    
    % allow for nonequispaced q   
    disp (['RH = ' num2str(RH(i)) '%...']);
    
    % plot reflectivity curve corresponding to RH(i)
    figure('Position',[1  1 0.9*scrsz(3) 0.9*scrsz(4)]);
    semilogy(1:1:length(Icorr(i,:)),Icorr(i,:),'-','color','blue');
    xlim([0 length(Icorr(i,:))]);
    xlabel('data points','interpreter','latex','fontsize',16);
    ylabel('corrected intensity (cps)','interpreter','latex','fontsize',16);
   
    title(['RH = ' num2str(RH(i)) '%. Click on the 1st Bragg peak'],'fontsize',16);
    %x_pos_1_pts = 110; semilogy(1:1:length(Icorr(i,:)),Icorr(i,:),'-','color','blue');
    xlim([0 length(Icorr(i,:))]);
    xlabel('data points','interpreter','latex','fontsize',16);
    ylabel('corrected intensity (cps)','interpreter','latex','fontsize',16);
    [x_pos_1_pts, y_pos] = ginput(1);
       
    for j = 1:num_peaks
        x_pos_j_pts = round(j*x_pos_1_pts); % search start position
      
        v = max(0,(x_pos_j_pts - peakwidth)):1:min((x_pos_j_pts + peakwidth),length(Icorr(i,:))); % search interval
        [peakheight, peakpos] = max(Icorr(i,v));
        peakpos = peakpos + x_pos_j_pts - peakwidth; % this offset is really important
        v = max(0,(peakpos-peakwidth)):1:min((peakpos+peakwidth),length(Icorr(i,:)));
        line([peakpos peakpos],[1 peakheight],'color','black'); % mark the peak position
        
        % extract corresponding momentum transfer 
        q_z_sampled(i,j) = q_z(i,peakpos);
        
        % create linear background in integration interval v and 
        % and subtract it from array 'y', integrate remaining peak area
        y = Icorr(i,v); 
        y_lborder = mean(y([1 2 3]));
        y_rborder = mean(y([length(y)-2 length(y)-1 length(y)]));
        backgr = linspace(y_lborder,y_rborder,length(y));
        line([min(v) max(v)],[y_lborder y_rborder],'color','black')
        y = y - backgr; 
        peakareas(i,j) = trapz(y); % trapezoidal integration
        %dummy = input('Press return to continue');
    end
    dummy = input('Press return to continue');
    close all
end

% determine d-spacings and q_z-intevals at each RH value
v = (q_z_sampled./(ones(num_RH,1)*(1:1:num_peaks)))'; % q/j " project back on first peak"
v([1 num_peaks],:) = [];  % removes first and last
% Do not take into account positions of 1st and last Bragg peak:
% The former show systematic shift towards higher q_z (presumably because 
% slope of reflectivity curve is still quite strong), the latter is often
% very weak
d_obs(:,1) = (2*pi./mean(v))';                      % d-spacing
d_obs(:,2) = (2*pi./(mean(v)).*(std(v)./mean(v)))'; % standard deviation
disp('Determining d-spacings ... done');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Lorentz correction and normalization                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lorentz correction: multiply peak area with peak order or squared peak
% order
disp('Which type of Lorentz correction would you like to apply?');
disp('F=sqrt(I*n)   (oriented samples) (1)');
disp('F=sqrt(I*n^2) (powder samples)   (2)');
choice_Lorentz = input('My choice: ');
peakorders = linspace(1,num_peaks,num_peaks);
for i = 1:num_RH
    peakareas_Lorentz(i,:) = peakareas(i,:).*(peakorders.^choice_Lorentz);
end
disp('Applying Lorentz correction ... done');

% normalization procedure:
% make sum of peak intensities proportional to stacking distance
% For motivation: see thesis P69
for i = 1:num_RH
    c_norm(i) = (d_obs(i,1)/d_obs(1,1))*sum(peakareas_Lorentz(1,:))./sum(peakareas_Lorentz(i,:));
    peakareas_normalized(i,:) = c_norm(i)*peakareas_Lorentz(i,:);
end
disp('Normalization of the different datasets ... done');

% convert corrected and normalized intensities into form factors
% remove imaginary values (do sometimes appear for weak and noisy peaks)
F_obs = real(sqrt(peakareas_normalized/max(max(peakareas_normalized))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           data output                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data output: form factors , RH values, d-spacings, choice of Lorentz
% factor
save (file_output,'F_obs','d_obs','RH','choice_Lorentz');
disp(['Data saved into ' file_output '.mat']);

save (file_output,'F_obs','d_obs','RH','choice_Lorentz');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
