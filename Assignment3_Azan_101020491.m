%%
%Name: Mohamed Azan
%ID: 101020491

clc
clear
%% Defining Constants
C.m0 = 9.109E-31; % Elctron Rest Mass
C.k = 1.381E-23; % Boltzman Constant 
C.m = C.m0 * 0.26; % Electron Effictive Mass
C.q = 1.602E-19; % Charge of an Electron
%% Defining Variables 
x_max = 200E-9; % Length of area 
y_max = 100E-9; % Width of Area 
Temp = 300; % Room Temp. 300K
tau = 0.2E-12; % Mean Time Between Collisions = 0.2ps

No_particles = 1000;
iterations = 3000;

appliedVoltage_x = 0.2;
appliedVoltage_y = 0;

vx = zeros(1,No_particles);
vy = zeros(1,No_particles);

figure(1)
subplot(2,1,1)
axis([0 200*10^-9 0 100*10^-9])
xlabel('')
ylabel('')

subplot(2,1,2)
title('Drift Current Density vs Time')
xlabel('Time (4ps)')
ylabel('Current Desnity (A/m)')
hold on
% Initial Conditions 
vT = sqrt((3*C.k*Temp)/C.m);

% Since tau = lambda/velocity
lambda = tau * vT; % Mean Free Path

x = x_max * rand(1,No_particles);
y = y_max * rand(1,No_particles); 

dt = 1e-15; % Time step is set to 1fs.


random_angles = randi([0 360],1,No_particles);
vx = vT*cosd(random_angles);
vy = vT*sind(random_angles);

%% Question 1a: What is the electric field on the electrons? 

Ex = appliedVoltage_x/x_max; 
Ey = appliedVoltage_y/y_max; 

fprintf('The electric field on the electrons is %i',Ex);
fprintf('\n');

%% Question 1b: what is the force on the electrons? 

Fx = Ex * C.q;
Fy = Ey * C.q;

fprintf('The force on the electrons is %i',Fx);
fprintf('\n');


%% Question 1c: calculate the accelration of the electrons. 

ax = Fx/C.m; 
ay = Fy/C.m;

fprintf('The accelration of the electrons is %i',ax);
fprintf('\n');

Jn_trace = zeros(1,iterations);
color = rand(1,20);
for i = 1:iterations 
    
    vx = vx + (ax*dt); 
    vy = vy + (ay*dt);
    x = x + vx*dt;
    y = y + vy*dt;
    
    x_out_right = x > x_max;
    x(x_out_right) = x(x_out_right) - x_max;
    
    x_out_left = x < 0;
    x(x_out_left) = x(x_out_left) + x_max;
    
    y_out_top = y > y_max;
    y(y_out_top) = y_max;
    vy(y_out_top) = -vy(y_out_top);
    
    y_out_bottom = y < 0;
    y(y_out_bottom) = 0;
    vy(y_out_bottom) = -vy(y_out_bottom);
    
    p_scat = 1 - exp(-dt/tau);
    
    random = rand(1,No_particles);
    
    scat = random < p_scat;
 
    new_theta = randi([0 360], No_particles);
    
     vx(scat) = vT * cosd(new_theta(scat));
     vy(scat) = vT * sind(new_theta(scat));
    
%% Question 1d: 

% Average carrier velocity 
% Since we are investigating the drift current
%in the x-direction, only the x-component of the
% velocity matters 

v_avg = mean(sqrt((vx.^2)+(vy.^2)));

v_avgx = mean(vx); 

% Drift Current Desnity Calculation 
% Jn = -enVnE  
% where e = carrier charge
%       n = carrier concentration 
%      Vn = aerage carrier velocity
%  and  E = Electric Field 

Jn = C.q * (10E15) * v_avgx * Ex;

mfp = v_avg * tau;
        for j = 1:20
            x_trace(j,i) = x(j);
            y_trace(j,i) = y(j);
        end
           
        Jn_trace(i) = Jn;
         
   
  
end
 figure(1)
 subplot(2,1,1)
 title(['Mean Free Path = ',num2str(mfp)])
 hold on
 for i=1:20
 scatter(x_trace(i,:),y_trace(i,:),2);
 end

           
 figure(1)
 subplot(2,1,2)
 xtime = linspace(1,iterations,iterations);
 plot(xtime,Jn_trace(1,:),'-');

 %% Question 1e
 
 %Desnity Plot 
 figure(2)
 histogram2(x,y,20);
 title('Density Plot')
 xlabel('X')
 ylabel('Y')
 
 %Temp PLot 

 %% Part 2
clc
clear
 
nx = 100;
ny = 200;
G = sparse(nx*ny);
B = zeros(1, nx*ny);
delta = 1; 

V0 = 1; 

V = zeros(nx,ny);
Ex = zeros(nx,ny);
Ey = zeros(nx,ny);
E = zeros(nx,ny);


conductance = ones(nx,ny); % create matrix

% Define Bottle nick 

top_box = [40 50 1 80];
bottom_box = [40 50 120 200];

% Define conductivty of each region

for i=1:nx
    for j=1:ny
        if (i > top_box(1) && i < top_box(2) && (j < top_box(4) || j > bottom_box(3)))
        conductance(i,j) = 0.01;
        end
    end
end


for i=1:nx
    for j=1:ny
        
        n = j + (i-1)*ny;
        nxm = j + (i-2)*ny;
        nxp = j + (i)*ny;
        nym = (j-1) + (i-1)*ny;
        nyp = (j+1) + (i-1)*ny;
        
        if i == 1
            
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = V0;
            
        elseif i == nx
            
            G(n,:) = 0;            
            G(n,n) = 1;
            B(n) = 0;
            
        elseif j == 1
            
            G(n,nyp) = (conductance(i,j+1)+conductance(i,j))/2;
            G(n,nxp) = (conductance(i+1,j)+conductance(i,j))/2;
            G(n,nxm) = (conductance(i-1,j)+conductance(i,j))/2;
            G(n,n) = -G(n,nxp)-G(n,nxm)-G(n,nyp);
            
        elseif j == ny
            
            G(n,nym) = (conductance(i,j-1)+conductance(i,j))/2;
            G(n,nxp) = (conductance(i+1,j)+conductance(i,j))/2;
            G(n,nxm) = (conductance(i-1,j)+conductance(i,j))/2;
            G(n,n) = -G(n,nxp)-G(n,nxm)-G(n,nym);
        else 
           
            G(n,nym) = (conductance(i,j-1)+conductance(i,j))/2;
            G(n,nyp) = (conductance(i,j+1)+conductance(i,j))/2;
            G(n,nxp) = (conductance(i+1,j)+conductance(i,j))/2;
            G(n,nxm) = (conductance(i-1,j)+conductance(i,j))/2;
            G(n,n) = -G(n,nxp)-G(n,nxm)-G(n,nyp)-G(n,nym);

        end
    end
end

V = G\B';

mapped = zeros(ny, nx);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        mapped(i,j) = V(n);
    end
end

%V(x,y) Surface Plot
figure(4)
surf(mapped)
title('Voltage Map')

% Question b: Calculate Electric Field
% The Electric field can be caluclated from the graidiant
% of the voltage disterbution 

[X,Y] = gradient(mapped);
figure(5)
quiver(X,Y);
title('Electric Field')
xlim([-10 200]);
ylim([-10 100])

%Question c
C.m0 = 9.109E-31; % Elctron Rest Mass
C.k = 1.381E-23; % Boltzman Constant 
C.m = C.m0 * 0.26; % Electron Effictive Mass
C.q = 1.602E-19; % Charge of an Electron

% Defining Variables 
x_max = 200E-9; % Length of area 
y_max = 100E-9; % Width of Area 
Temp = 300; % Room Temp. 300K
tau = 0.2E-12; % Mean Time Between Collisions = 0.2ps

No_particles = 1000;
iterations = 3000;

vx = zeros(1,No_particles);
vy = zeros(1,No_particles);

vT = sqrt((3*C.k*Temp)/C.m);

% Since tau = lambda/velocity
lambda = tau * vT; % Mean Free Path

x = x_max * rand(1,No_particles);
y = y_max * rand(1,No_particles); 

dt = 1e-15; % Time step is set to 1fs.


random_angles = randi([0 360],1,No_particles);
vx = vT*cosd(random_angles);
vy = vT*sind(random_angles);

%  What is the electric field on the electrons? 

Ex = X;
Ey = Y;

fprintf('The electric field in the X on the electrons is %i',Ex);
fprintf('\n');
fprintf('The electric field in the Y on the electrons is %i',Ey);
fprintf('\n');

% calculate the accelration of the electrons. 
Fx = Ex * C.q;
Fy = Ey * C.q;

ax = Fx/C.m; 
ay = Fy/C.m;

fprintf('The accelration of the electrons is %i',ax);
fprintf('\n');

Jn_trace = zeros(1,iterations);
color = rand(1,20);

% if electron apwns in box, relocate to middle of region 
for i=1:No_particles
    if (x(i) > top_box(1) && x(i) < top_box(2))
        if (y(i) < top_box(4) || y(i) > bottom_box(3))
            x(i) = x_max/2;
            y(i) = y_max/2;
        end
    end
end

for i = 1:iterations 
    
    for j = 1:No_particles 
        if (y(j) < 1e-6)
            index_y(j) = 1;
        elseif (y(j) == y_max)
            y(j) = y_max * 10^9;
        else
            index_y(j) = ceil(y(j)*10^9); 
        end
        
        if (x(j) < 1e-6)
            index_x(j) = 1;
        elseif (x(j) == x_max)
            x(1) = x_max * 10^9;
        else
            index_x(j) = ceil(x(j)*10^9); 
        end
            
    end
    
     
        
    vx = vx + ax(index_x,index_y).*dt; 
    vy = vy + ay(index_x,index_y).*dt;
    
    x = x + vx*dt;
    y = y + vy*dt;
    
 % If an electron hits the side of the box, the electron 
 % will elastically rebound. 
 
   for ii=1:No_particles
    if (x(ii) == top_box(1) || x(ii) == top_box(2))
        if (y(ii) < top_box(4) || y(ii) > bottom_box(3))
            vx(ii) = vx(ii) * -1;
        end
    end
   
    
    if(y(ii) == top_box (4) || y(ii) == bottom_box(3))
        if(x(ii) > top_box(1) && x(ii) < top_box(2))
            y(ii) = y(ii) * -1;
        end
    end
   end
    
    x_out_right = x > x_max;
    x(x_out_right) = x(x_out_right) - x_max;
    
    x_out_left = x < 0;
    x(x_out_left) = x(x_out_left) + x_max;
    
    y_out_top = y > y_max;
    y(y_out_top) = y_max;
    vy(y_out_top) = -vy(y_out_top);
    
    y_out_bottom = y < 0;
    y(y_out_bottom) = 0;
    vy(y_out_bottom) = -vy(y_out_bottom);
    
    p_scat = 1 - exp(-dt/tau);
    
    random = rand(1,No_particles);
    
    scat = random < p_scat;
 
    new_theta = randi([0 360], No_particles);
    
     vx(scat) = vT * cosd(new_theta(scat));
     vy(scat) = vT * sind(new_theta(scat));
   

% Average carrier velocity 
% Since we are investigating the drift current
%in the x-direction, only the x-component of the
% velocity matters 

v_avg = mean(sqrt((vx.^2)+(vy.^2)));

v_avgx = mean(vx); 

% Drift Current Desnity Calculation 
% Jn = -enVnE  
% where e = carrier charge
%       n = carrier concentration 
%      Vn = aerage carrier velocity
%  and  E = Electric Field 


mfp = v_avg * tau;
        for j = 1:20
            x_trace(j,i) = x(j);
            y_trace(j,i) = y(j);
        end
end     
     %   Jn_trace(i) = Jn;

 figure(7)
 figure('Renderer', 'painters', 'Position',  [10 10 1100 600])
 hold on
 
line([40 80],[50 80])
line([40 1],[40 80])
line([50 1], [50 80])

   
line([40 120],[50 120])
line([40 200],[40 120])
line([50 200], [50 120])

 title(['Mean Free Path = ',num2str(mfp)])
 for i=1:20
 scatter(x_trace(i,:),y_trace(i,:),2);
 end