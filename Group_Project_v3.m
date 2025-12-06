clear
%----Given Values-------
d = 84 * 10^-6;     %Diameter of droplet (m)
r = 42 * 10^-6;     %Radius of droplet (m)
D = 3 *10^-3;       %Distance to paper (m)
vx = 20;            %Velocity in x direction (m/s)
L1 = 0.5 * 10^-3;   %Length of Capacitor (m)
L2 = 1.25 * 10^-3;  %Distance between paper and cap (m)
q = -1.9*10^-10;    %Charge of electron (C)
pd = 1000;          %Density of droplet (kg/m^3)
W = 1 *10^-3;       %Diameter between Cap plates (m)
Volt = 0;

%------Equations------
m = pd * (4/3)*pi*(r)^3;  %Mass of droplet

%----Condition values----
dt = .000001;  %Time step
x(1) = 0;
z(1) = 0;
o = 1;

%The letter I at font 11 and 300 dpi will be 3.88mm tall

%Points is where the center of the droplet has to land on the paper
%Found by 3.88mm / 2 = 1.94mm --> 1.94mm - .042mm = 1.898mm
points(1) = 1.898*10^-3;
%Example: points(2) = points(1) - radius of droplet
while points(o) > -1.94*10^-3
    o = o + 1;
    points(o) = points(o-1) - d;
end

numDistPoints = floor((D/dt)/20);

%--------- Part 2: Varried Voltages ---------------------------------------

%Calculates the voltage for every location the droplet has to land on the
%paper
for it = 1:length(points)
    V(it) = (m*points(it)*W*(vx)^2)/(q*L1*L2);
end

j = zeros(length(points),numDistPoints);   %Initialize the array for plotting the z-axis

for test = 1:length(points)
    Vol = V(test);
    E = Vol / W;
    % ay = q*E/m;

    %Reset indexes
    x(1) = 0;
    z(1) = 0;
    vy = 0;

    for i = 2:numDistPoints      %iterations in mm
        % Update index

        beforeCap = x(i-1) < (D-L1-L2);
        afterCap = x(i-1) > (D-L2);

        % Calculate new velocities (only y changes)
        if ~beforeCap && ~afterCap
            ay = q*E/m;
            vy = vy + ay * dt;
        else
            ay = 0;
        end

        % Calculate new positions

        x(i) = x(i-1) + vx * dt;
        z(i) = z(i-1) + vy * dt + ay * dt^2;
        j(test,i) = z(i);
    
    end
end
y = zeros(size(x)); % Initialize z positions for 3D plot


%-------------Capacitor Display-------------------

InkMotionGraph(x,y,j,points);

%--------------- Part 3: Voltage graph -----------------------------

%%
%--------------- Part 4: Adjusted Parameters -----------------------
%--------- Part 4a -------------

%----Given Values-------
d = 84 * 10^-6;     %Diameter of droplet (m)
r = 42 * 10^-6;     %Radius of droplet (m)
vx = 20;            %Velocity in x direction (m/s)
L1 = 0.5 * 10^-3;   %Length of Capacitor (m)
L2 = 3.75 * 10^-3;  %Distance between paper and cap (m)
D = 1.25 * 10^-3 + L1 + L2;  %Distance to paper (m)
q = -1.9*10^-10;    %Charge of electron (C)
pd = 1000;          %Density of droplet (kg/m^3)
W = 1 *10^-3;       %Diameter between Cap plates (m)

%------Equations------
T = L1 / vx;        %Droplet time in capacitor
m = pd * (4/3)*pi*(r)^3;  %Mass of droplet

%----Condition values----
dt = .000001;  %Time step
o = 1;

% calculated max Voltage = 2613.389
CapTime = L1 / vx; % s
TimeToPaper = L2 / vx;

maxV = floor((W^2 * m) / (CapTime^2 * q)); % V

E = maxV / W;
ay = q*E/m;
vy = ay * CapTime;

peak = (ay * CapTime^2) + vy * TimeToPaper;

clear points
points(1) = peak;
while points(o) > -peak
    o = o + 1;
    points(o) = points(o-1) - d;
end

numDistPoints = floor((D/dt)/20);

clear V
for it = 1:length(points)
    V(it) = (m*points(it)*W*(vx^2))/(q*L1*L2);
end

clear x
clear z
clear y
clear j
clear test

for test = 1:length(points)
    Vol = V(test);
    E = Vol / W;
    % ay = q*E/m;

    %Reset indexes
    x(1) = 0;
    z(1) = 0;
    vy = 0;

    for i = 2:numDistPoints      %iterations in mm
        % Update index

        beforeCap = x(i-1) < (D-L1-L2);
        afterCap = x(i-1) > (D-L2);

        % Calculate new velocities (only y changes)
        if ~beforeCap && ~afterCap
            ay = q*E/m;
            vy = vy + ay * dt;
        else
            ay = 0;
        end

        % Calculate new positions

        x(i) = x(i-1) + vx * dt;
        z(i) = z(i-1) + vy * dt + ay * dt^2;
        j(test,i) = z(i);
    
    end
end
y = zeros(size(x)); % Initialize z positions for 3D plot

%-------------Capacitor Display-------------------

InkMotionGraph(x,y,j,points);
%%
%--------- Part 4b -------------

%----Given Values-------
d = 84 * 10^-6;     %Diameter of droplet (m)
r = 42 * 10^-6;     %Radius of droplet (m)
vx = 20;            %Velocity in x direction (m/s)
L1 = 1 * 10^-3;     %Length of Capacitor (m)
L2 = 1.25 * 10^-3;  %Distance between paper and cap (m)
D = 1.25 * 10^-3 + L1 + L2;  %Distance to paper (m)
q = -1.9*10^-10;    %Charge of electron (C)
pd = 1000;          %Density of droplet (kg/m^3)
W = 1 *10^-3;       %Diameter between Cap plates (m)

%------Equations------
m = pd * (4/3)*pi*(r)^3;  %Mass of droplet

%----Condition values----
dt = .000001;  %Time step

o = 1;

% Calculate max Voltage
CapTime = L1 / vx; % s
TimeToPaper = L2 / vx;

maxV = floor((W^2 * m) / (CapTime^2 * q)); % V

E = maxV / W;
ay = q*E/m;
vy = ay * CapTime;

peak = (ay * CapTime^2) + vy * TimeToPaper;


clear points
points(1) = peak;
while points(o) > -peak
    o = o + 1;
    points(o) = points(o-1) - d;
end

numDistPoints = floor((D/dt)/20);

clear V
for it = 1:length(points)
    V(it) = (m*points(it)*W*(vx)^2)/(q*L1*L2);
end

clear x
clear z
clear y
clear j
clear test

for test = 1:length(points)
    Vol = V(test);
    E = Vol / W;
    % ay = q*E/m;

    %Reset indexes
    x(1) = 0;
    z(1) = 0;
    vy = 0;

    for i = 2:numDistPoints      %iterations in mm
        % Update index

        beforeCap = x(i-1) < (D-L1-L2);
        afterCap = x(i-1) > (D-L2);

        % Calculate new velocities (only y changes)
        if ~beforeCap && ~afterCap
            ay = q*E/m;
            vy = vy + ay * dt;
        else
            ay = 0;
        end

        % Calculate new positions

        x(i) = x(i-1) + vx * dt;
        z(i) = z(i-1) + vy * dt + ay * dt^2;
        j(test,i) = z(i);
    
    end
end
y = zeros(size(x)); % Initialize z positions for 3D plot

%-------------Capacitor Display-------------------

InkMotionGraph(x,y,j,points);
%%
%--------- Part 4c -------------

%----Given Values-------
d = 840 * 10^-6;     %Diameter of droplet (m)
r = 42 * 10^-6;     %Radius of droplet (m)
vx = 20;            %Velocity in x direction (m/s)
L1 = 0.5 * 10^-3;     %Length of Capacitor (m)
L2 = 1.25 * 10^-3;  %Distance between paper and cap (m)
D = 1.25 * 10^-3 + L1 + L2;  %Distance to paper (m)
q = -1.9*10^-10;    %Charge of electron (C)
pd = 1000;          %Density of droplet (kg/m^3)
W = 1 *10^-3;       %Diameter between Cap plates (m)

%------Equations------
m = pd * (4/3)*pi*(r)^3;  %Mass of droplet

%----Condition values----
dt = .000001;  %Time step
o = 1;

% Calculate max Voltage
CapTime = L1 / vx; % s
TimeToPaper = L2 / vx;

maxV = floor((W^2 * m) / (CapTime^2 * q)); % V

E = maxV / W;
ay = q*E/m;
vy = ay * CapTime;

peak = (ay * CapTime^2) + vy * TimeToPaper;


clear points
points(1) = peak;
while points(o) > -peak
    o = o + 1;
    points(o) = points(o-1) - d;
end

numDistPoints = floor((D/dt)/20);

clear V
for it = 1:length(points)
    V(it) = (m*points(it)*W*(vx)^2)/(q*L1*L2);
end

clear x
clear z
clear y
clear j
clear test

for test = 1:length(points)
    Vol = V(test);
    E = Vol / W;
    % ay = q*E/m;

    %Reset indexes
    x(1) = 0;
    z(1) = 0;
    vy = 0;

    for i = 2:numDistPoints      %iterations in mm
        % Update index

        beforeCap = x(i-1) < (D-L1-L2);
        afterCap = x(i-1) > (D-L2);

        % Calculate new velocities (only y changes)
        if ~beforeCap && ~afterCap
            ay = q*E/m;
            vy = vy + ay * dt;
        else
            ay = 0;
        end

        % Calculate new positions

        x(i) = x(i-1) + vx * dt;
        z(i) = z(i-1) + vy * dt + ay * dt^2;
        j(test,i) = z(i);
    
    end
end
y = zeros(size(x)); % Initialize z positions for 3D plot

%-------------Capacitor Display-------------------

InkMotionGraph(x,y,j,points);
%%
%--------- Part 4d -------------

%----Given Values-------
d = 84 * 10^-6;     %Diameter of droplet (m)
r = 42 * 10^-6;     %Radius of droplet (m)
vx = 20 * 2;            %Velocity in x direction (m/s)
L1 = 0.5 * 10^-3;     %Length of Capacitor (m)
L2 = 1.25 * 10^-3;  %Distance between paper and cap (m)
D = 1.25 * 10^-3 + L1 + L2;  %Distance to paper (m)
q = -1.9*10^-10;    %Charge of electron (C)
pd = 1000;          %Density of droplet (kg/m^3)
W = 1 *10^-3;       %Diameter between Cap plates (m)

%------Equations------
m = pd * (4/3)*pi*(r)^3;  %Mass of droplet

%----Condition values----
dt = .000001;  %Time step
o = 1;

% Calculate max Voltage
CapTime = L1 / vx; % s
TimeToPaper = L2 / vx;

maxV = floor((W^2 * m) / (CapTime^2 * q)); % V

E = maxV / W;
ay = q*E/m;
vy = ay * CapTime;

peak = (ay * CapTime^2) + vy * TimeToPaper;


clear points
points(1) = peak;
while points(o) > -peak
    o = o + 1;
    points(o) = points(o-1) - d;
end

numDistPoints = floor((D/dt)/20);

clear V
for it = 1:length(points)
    V(it) = (m*points(it)*W*(vx)^2)/(q*L1*L2);
end

clear x
clear z
clear y
clear j
clear test

for test = 1:length(points)
    Vol = V(test);
    E = Vol / W;
    % ay = q*E/m;

    %Reset indexes
    x(1) = 0;
    z(1) = 0;
    vy = 0;

    for i = 2:numDistPoints      %iterations in mm
        % Update index

        beforeCap = x(i-1) < (D-L1-L2);
        afterCap = x(i-1) > (D-L2);

        % Calculate new velocities (only y changes)
        if ~beforeCap && ~afterCap
            ay = q*E/m;
            vy = vy + ay * dt;
        else
            ay = 0;
        end

        % Calculate new positions

        x(i) = x(i-1) + vx * dt;
        z(i) = z(i-1) + vy * dt + ay * dt^2;
        j(test,i) = z(i);
    
    end
end
y = zeros(size(x)); % Initialize z positions for 3D plot

%-------------Capacitor Display-------------------

InkMotionGraph(x,y,j,points);
%%
%--------- Part 4e -------------

%----Given Values-------
d = 84 * 10^-6;     %Diameter of droplet (m)
r = 42 * 10^-6;     %Radius of droplet (m)
vx = 20;            %Velocity in x direction (m/s)
L1 = 0.5 * 10^-3;     %Length of Capacitor (m)
L2 = 1.25 * 10^-3;  %Distance between paper and cap (m)
D = 1.25 * 10^-3 + L1 + L2;  %Distance to paper (m)
q = -1.9*10^-10 * 5;    %Charge of electron (C)
pd = 1000;          %Density of droplet (kg/m^3)
W = 1 *10^-3;       %Diameter between Cap plates (m)

%------Equations------
m = pd * (4/3)*pi*(r)^3;  %Mass of droplet

%----Condition values----
dt = .000001;  %Time step
i = 1;
u = 1;
o = 1;

% Calculate max Voltage
CapTime = L1 / vx; % s
TimeToPaper = L2 / vx;

maxV = floor((W^2 * m) / (CapTime^2 * q)); % V

E = maxV / W;
ay = q*E/m;
vy = ay * CapTime;

peak = (ay * CapTime^2) + vy * TimeToPaper;


clear points
points(1) = peak;
while points(o) > -peak
    o = o + 1;
    points(o) = points(o-1) - d;
end

numDistPoints = floor((D/dt)/20);

clear V
for it = 1:length(points)
    V(it) = (m*points(it)*W*(vx)^2)/(q*L1*L2);
end

clear x
clear z
clear y
clear j
clear test

for test = 1:length(points)
    Vol = V(test);
    E = Vol / W;
    % ay = q*E/m;

    %Reset indexes
    x(1) = 0;
    z(1) = 0;
    vy = 0;

    for i = 2:numDistPoints      %iterations in mm
        % Update index

        beforeCap = x(i-1) < (D-L1-L2);
        afterCap = x(i-1) > (D-L2);

        % Calculate new velocities (only y changes)
        if ~beforeCap && ~afterCap
            ay = q*E/m;
            vy = vy + ay * dt;
        else
            ay = 0;
        end

        % Calculate new positions

        x(i) = x(i-1) + vx * dt;
        z(i) = z(i-1) + vy * dt + ay * dt^2;
        j(test,i) = z(i);
    
    end
end
y = zeros(size(x)); % Initialize z positions for 3D plot

%-------------Capacitor Display-------------------

InkMotionGraph(x,y,j,points);
%%
%--------------- Part 5: Drawing 'H' -----------------------

%----Given Values-------
d = 84 * 10^-6;     %Diameter of droplet (m)
r = 42 * 10^-6;     %Radius of droplet (m)
vx = 20;            %Velocity in x direction (m/s)
L1 = 0.5 * 10^-3;     %Length of Capacitor (m)
L2 = 1.25 * 10^-3;  %Distance between paper and cap (m)
D = 1.25 * 10^-3 + L1 + L2;  %Distance to paper (m)
q = -1.9*10^-10;    %Charge of electron (C)
pd = 1000;          %Density of droplet (kg/m^3)
W = 1 *10^-3;       %Diameter between Cap plates (m)

%------Equations------
m = pd * (4/3)*pi*(r)^3;  %Mass of droplet

%----Condition values----
dt = .000001;  %Time step
i = 1;
u = 1;
o = 1;

% Calculate max Voltage
CapTime = L1 / vx; % s
TimeToPaper = L2 / vx;

maxV = floor((W^2 * m) / (CapTime^2 * q)); % V

E = maxV / W;
ay = q*E/m;
vy = ay * CapTime;

peak = (ay * CapTime^2) + vy * TimeToPaper;


clear points
points(1) = peak;
while points(o) > -peak
    o = o + 1;
    points(o) = points(o-1) - d;
end

numDistPoints = floor((D/dt)/20);

clear V
for it = 1:length(points)
    V(it) = (m*points(it)*W*(vx)^2)/(q*L1*L2);
end

clear x
clear z
clear y
clear j
clear k
clear test

state1 = false;
state2 = false;
state3 = false;

for h = 1:3
    
    if (h == 1) % First Leg of H
        state1 = true;
        state2 = false;
        state3 = false;
    elseif h == 2
        state1 = false;
        state2 = true;
        state3 = false;
    elseif h == 3
        state1 = false;
        state2 = false;
        state3 = true;
    end
    
    for test = 1:length(points)
        
        if state1
            VolY = V(test);
            VolZ = V(1);
        elseif state2
            VolY = 0;
            VolZ = V(test);
        elseif state3
            VolY = V(test);
            VolZ = V(end);
        end
        EY = VolY / W;
        EZ = VolZ / W;
    
        %Reset indexes
        x(1) = 0;
        y(1) = 0;
        z(1) = 0;
    
        vy = 0;
        vz = 0;
    
        for i = 2:numDistPoints      %iterations in mm
            % Update index
    
            beforeCap = x(i-1) < (D-L1-L2);
            afterCap = x(i-1) > (D-L2);
    
            % Calculate new velocities
            if ~beforeCap && ~afterCap
                ay = q*EY/m;
                az = q*EZ/m;
                vy = vy + ay * dt;
                vz = vz + az * dt;
            else
                ay = 0;
                az = 0;
            end
    
            % Calculate new positions
    
            x(i) = x(i-1) + vx * dt;
            y(i) = y(i-1) + vy * dt + ay * dt^2;
            z(i) = z(i-1) + vz * dt + az * dt^2;
            j(test,i,h) = y(i);
            k(test,i,h) = z(i);
        
        end
    end
end

%-------------Capacitor Display-------------------
width = .2*10^-3;
height = .05*10^-3;
Fa1 = 1.25*10^-3;     %First face
point1 = 0.5*10^-3;
point2 = 1.75*10^-3;

% Define the vertices of the cuboid
%% ADD THE OTHER 2 SHAPES %%
P = [Fa1, width/2, point1;               %V1
     point2, width/2, point1;
     point2, -width/2, point1;
     Fa1, -width/2, point1;           %V4
     Fa1, width/2, height+point1;          %V5
     point2, width/2, height+point1;
     point2, -width/2, height+point1;
     Fa1, -width/2, height+point1];     %V8

Q = [Fa1, width/2, -point1;               %V1
     point2, width/2, -point1;
     point2, -width/2, -point1;
     Fa1, -width/2, -point1;           %V4
     Fa1, width/2, -height-point1;          %V5
     point2, width/2, -height-point1;
     point2, -width/2, -height-point1;
     Fa1, -width/2, -height-point1];     %V8

% Create the alphaShape object
alpha_radius = 3*10^-3;
Tcap = alphaShape(P, alpha_radius);
Bcap = alphaShape(Q, alpha_radius);

%----------------------------------------------

figure; %[output:6c0cd270]
%plot3(x,k(:,:,1),j(:,:,1))
hold on; %[output:6c0cd270]
% plot3(x,k(:,:,2),j(:,:,2))
% plot3(x,k(:,:,3),j(:,:,3))
view(3); %[output:6c0cd270]
box on; %[output:6c0cd270]
grid on; %[output:6c0cd270]
title('Ink Motion'); %[output:6c0cd270]
xlabel('Horizontal Distance (mm)'); %[output:6c0cd270]
zlabel('Vertical Distance (mm)'); %[output:6c0cd270]

plot(Tcap); %[output:6c0cd270]
plot(Bcap); %[output:6c0cd270]
xlim([0 (x(end) + x(end)/4)]); % Makes graph look nicer %[output:6c0cd270]

% Animate the projectile
for p = 1:3
    for l = 1:length(points)
        h = plot3(x(1), k(l,1), j(l,1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); %[output:6c0cd270]
        for u = 1:length(x)
    
            h.XData = x(u);
            h.YData = k(l,u, p);
            h.ZData = j(l,u, p);
            %drawnow; 
        end
    end
end
hold off; %[output:6c0cd270]

%%
%------------------- Functions -------------------------
function c = InkMotionGraph(x_, y_, j_, points_)
    width = .2*10^-3;
    height = .05*10^-3;
    Fa1 = 1.25*10^-3;     %First face
    point1 = 0.5*10^-3;
    point2 = 1.75*10^-3;
    
    % Define the vertices of the cuboid
    P = [Fa1, width/2, point1;               %V1
         point2, width/2, point1;
         point2, -width/2, point1;
         Fa1, -width/2, point1;           %V4
         Fa1, width/2, height+point1;          %V5
         point2, width/2, height+point1;
         point2, -width/2, height+point1;
         Fa1, -width/2, height+point1];     %V8
    
    Q = [Fa1, width/2, -point1;               %V1
         point2, width/2, -point1;
         point2, -width/2, -point1;
         Fa1, -width/2, -point1;           %V4
         Fa1, width/2, -height-point1;          %V5
         point2, width/2, -height-point1;
         point2, -width/2, -height-point1;
         Fa1, -width/2, -height-point1];     %V8
    
    % Create the alphaShape object
    alpha_radius = 3*10^-3;
    Tcap = alphaShape(P, alpha_radius);
    Bcap = alphaShape(Q, alpha_radius);
    
    %----------------------------------------------
    
    figure;
    plot3(x_,y_,j_)
    view(3);
    box on;
    grid on;
    title('Ink Motion');
    xlabel('Horizontal Distance (mm)');
    zlabel('Vertical Distance (mm)');
    hold on;
    
    plot(Tcap);
    plot(Bcap);
    xlim([0 (x_(end) + x_(end)/4)]); % Makes graph look nicer
    
    % Animate the projectile
    for l = 1:length(points_)
        h = plot3(x_(1), y_(1), j_(l,1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
        for k = 1:length(x_)
    
            h.XData = x_(k);
            h.ZData = j_(l,k);
            %drawnow;
        
        end
    end
    hold off;
    c = true;
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":39}
%---
%[output:6c0cd270]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7tvQuwVtWZJvzmcJGrwTNEwi2GFgL80fIQWyXdhmFIWixMYyQWVZDUL5NAjymK1rIQiEM0sTApbCaVaX5+LYEUqbLJNL9F2wzQPaHpHKmUEQIi7WWkhMEAolwCR0EOImP\/PKvP+madzb6svfbl25dnV53iHPbea73redd31nPe66daW1v\/VXgRASJABIgAESACRKBECHyKBKZE2qKoRIAIEAEiQASIgEKABIYbgQgQASJABIgAESgdAiQwpVMZBSYCRIAIEAEiQARIYLgHiAARIAJEgAgQgdIhQAJTOpVRYCJABIgAESACRIAEhnuACBABIkAEiAARKB0CJDClUxkFJgJEgAgQASJABEhguAeIQA0RWLdunYwdO1Z27NghixcvtkJg+fLlMmnSJOt39PMY\/Ny5c7J06VLZuXOnmmvBggUye\/Zs9f3Fixdl1apVsmHDBis5vA\/NnDlT5s+fL6dOnZKFCxfKoUOHnMbhS0SACJQLARKYcumL0hKBVBDIm8BA6PXr18vKlSuV\/Ca5iUtgRo0aJT\/+8Y\/l7\/7u7xTpIYFJZUtwECJQOgRIYEqnMgpMBJIjkCeBOX36tLS2tjYsNyAgK1askGHDhsW2wNx2222ybNky6d27dyKrTXIEOQIRIALNRoAEptka4PxEoAkIeAmM6R4aMmSIci\/hOnbsWMMt43UhaTdQkAVFP\/\/2228rsqJdPLfccoty+cCthGvAgAHdyIiWTcOi3VyavOB5fe3fv1+2bt3q60IKGgfvmlabXbt2yTe+8Y3YZKoJauOURIAIGAiQwHA7EIEaIhBEYPyg0ATCJDB79uxRpAFXUPyKfv53v\/udDB8+XK6++moVB3Prrbeq+Bf9\/4MHD1Zj4GfTMmPKAqLy1FNPKetLFIHBe2HjzJkzp0FgYMnxXiZpq+HW4JKJQGkQIIEpjaooKBFID4EgAmMG2+pnQB5w6HstKjj8zbgWr3R+Vh08f\/PNNysLz\/PPP6\/IjCYwsPyA2JgWHW0p0UTp97\/\/\/RUuJG8MzNe\/\/vXIcTAeCJi5Bm1R8gYcp4c6RyICRCBNBEhg0kSTYxGBkiAQRGA0WcEy9IHuJTB6iVHBtyaBOXr0aDerC6wxsJL8xV\/8RYPAgNj4ZTmZsm7cuDGSwNx\/\/\/2R42gLEtags6O0iwrrMzOmSqJSikkEaocACUztVM4FEwGRsBgYnVYdRWCAo0l4wiwwmjBolw3cNE8++aQsWrSoqQTGTL0mgeEngwiUCwESmHLpi9ISgVQQSEJgEBPz9NNPN+JMgtxIpgXGfB4LMMeI40LScTL6Hb806jguJBKYVLYTByECTUGABKYpsHNSItBcBJISGFhpouqveLOWzKwgkJ7NmzcrEmQbxIs4HG8Ktl8WEpC1DeIlgWnuPuTsRCAJAiQwSdDju0SgpAikQWCw9LB6MkFp1zpI9sSJE90IjK7EG5b+jDnNInhwRf3qV79S8TXeSry2adS6ei9dSCXdzBS7tgiQwNRW9Vw4ESACRIAIEIHyIkACU17dUfIYCJhF0OL0\/4kxBR8lAkSACBCBHBEggckRbE7VPATgdkAqL+IuHnvsMVUUTTcWbJ5UnJkIEAEiQARcESCBcUWO7+WOgA7gbG9vbzQFhBDe4FDdMNBPQFhivve978mPfvQjdi3OXYOckAgQASKQHgIkMOlhyZEyRMDMPjHTdlGrZPr06Y0S9fp7P+uKzpp56aWXRNc6yVBkDk0EiAARIAIZIkACkyG4HDodBMx0XVRw3bRpU8MCA+sLLjPFFhaa48ePX1EqXkuj3Ulhlpp0JC\/GKJ988okSpKWlpRgCUQoiQASIQAoIkMCkACKHSIYACMqUKVMEJeBx4ef77rtPHn\/8cRWnMnXqVPX\/HR0dqoy8JjDaKnPgwAFlUfH+bEoFooNqsCAt5vfJJC\/22x9\/\/LF80tkpd1+6JGNE5C0Reb6lRXr07y+9evUKFH7o0KEybdo0+dKXvtR4RjduLPaKKR0RIAJ1QoAEpk7aLvBaQSpAQPbt2yfjx4\/37UWjM4m8BMaMicE4sL54XURmFlJY+fsCQxRLtM7OTlnU2Sk\/9HkL\/7eif3+56qqrGndBWkBYNHF59913BV+axOD7LVu2yLPPPisfffRRLFn4MBEgAkQgCwRIYLJAlWM6IQDXzujRo0UXFvMOkoTAOAlU0pdgeXno7Flf8qKXBBLz04ED5XOf+5wiLXPnzlW3QFTWrl2ryArIy6pVq5Qr7stf\/rJ8+9vflpdffllZwfAcLyJABIhAMxEggWkm+py7gUASC4yNC6lOUH\/0wQdy7tKlyCX\/+4kTZfUvfqHICAjL1q1buxETTWBmzJgh58+fl7vvvlu+8Y1vCKw1a9asUUSHFxEgAkSgWQiQwDQLec7bQCAqBkY\/6LXA4P9Nl1FQmnWdoEbA7rSODvkbi0V\/S0ROfO1ryqpiXnAtDRkyRCZMmCBLlixR7ji49m644QbZvXu3ssTAYgPiA+sMrTEWYPMRIkAEUkeABCZ1SDlgVgj4EZg4adRZyVWkcUFgHujoCHUfaXnhRvqvgwap7CSQlk9\/+tNy7bXXqn8R5zJ8+HBFYJ544gk5e\/astLa2KisN7sEKA\/cSrTFF0j5lIQL1QoAEpl76LvVq\/QiMtsKMHTtWrc2sEVPqxToKH9cC8+uhQ5W1ZeTIkWpGkJPDhw8LGi1qF9KsWbOkX79+MmnSJNm1a5e8\/vrrjUDe7373uw1rDGJjvNYcx2XwNSJABIhAJAIkMJEQ8QEiUHwEdBYRLCj\/\/b\/8F9lvITIo38ivflV69OihCAuyt8wMI28MDAgkLDGw2Bw5ckQRHVyYe+nSpYrwMDbGAng+QgSIQCoIkMCkAiMHIQLNQUDXbNFZRO+88448MWWK3H5ZHL8Uai3lz0Rks4j8r1Gj5P333\/cV3ktgdAyMttiA7Lz66qu0xjRH9ZyVCNQeARKY2m8BAlA2BPxqtiCLaO\/evSrI9tqODvn+5erEbweQGJCXDhH5JYJ4u2Jg\/DAIIjAgLoiZufHGG9W\/tMaUbQdRXiJQDQRIYKqhR66i4gho0oLMoLvuuktl\/iDeBEG1ZtwJYmBu7uhQFpjJIvITERks0qjEi3yjxV0EBhaYva2tgciFERj9EurIIH4GpOatt95qWHMYG1PxDcnlEYECIEACUwAlUAQiEISA10Wka7Z4a7DoLCK0CPjoslvHa4GBNebzXZOkYYEx5cXcY8aMUdlLQdYYWIgQ5MuLCBABIpAWAiQwaSHJcYhASgiEuYiCarboLCIUnJPt23OzwJhLDrLGwGKEIF+QL2YqpbRJOAwRIAJCAsNNQAQKgICfi+iNN96QkydPytNPP90tO8ivZguyiBCMe+bMGd8YmCwtMHGtMSBhbA5ZgE1HEYhAyREggSm5Ail+uRHQLiJYKfC96SKCS0Zn\/mCV+FkXm0PMCQiLJi4ahbxjYILQD7LG6LiaIFdYubVJ6YkAEcgTARKYPNHmXESgq26Kt\/OzziIyXUQgK3ju6NGjqkIuLhAXkBZdg8ULKAiMXxZSXhaYIGsMZEaQLy4zrofNIfmRIAJEwBUBEhhX5PgeEYiBAA5t09qis4iQ+gzy4j34QV7Q2wlZR+hDhPousLgE1WwpmgXGXA\/IF4J8vZlKIGeIjWE7ghgbiY8SASLQQIAEhpuBCGSIQJiLyI+0mC4iWFOGDRsm27Zt6xYDEyZukSwwQdYYkLDXXnutcdtMuWZzyAw3I4cmAhVDgASmYgrlcpqPgF8WkV\/NFkiqOz\/DSoHvTReRGQNjlviPIjB51oGJi7Zek9cVxuaQcZHk80SACJDAcA8QgZQQ0HEtZqG5MBeR2flZB+OaLiJXAlOUGJggWDVpQ+o31ovYGE3QWAAvpc3IYYhADRAggamBkuu8xJkzZ8q8efNk9erVsmHDhtShCHIRoUIu4lz0hUO7T58+KhgXXzqLSGcS+QnmSmCKbIEx14n1ITaG7QhS35YckAjUAgESmFqouZ6LRBDsihUrZPDgwbJq1arUCExaLqIorbgSmKJbYLzrNlOudXNIrP2hhx6SO+64gwXwojYK7xOBmiJAAlNTxddh2QsWLJDp06erpaZhgQFxgYtDu4hgYYGlJSiLyHQR+dVsidKBK4EpiwXGXL+3OSTwQg0cZF89\/PDDKp18zZo14m2hEIUh7xMBIlBdBEhgqqvbWq\/stttuk0WLFkl7e7siMa4ExtZFBLB1BpGtiyhKQa4EpmwWGBMHbY3p2bOnci298MILyt3G2Jio3cL7RKB+CJDA1E\/ntVjx8uXL1Tr37NkTOwYmDRfR8ePHrVOfgxTiSmDKaIHxWmPa2tpUDRykkB88eFDdhl5QNwbWGDaHrMXHmIskAqEIkMBwg1QOAQTu3nPPPfLII4\/ILbfcYk1gTBcRQNGpz1m5iKKAdyUwZbbAaEyw9ilTpkhHR4ecPXtWZSrpDC02h4zaObxPBOqBAAlMPfRcq1WuW7dOWV5WrlwptllIug4JgAJh8WYR4f\/h0oCLQ7uILly4oMr64yuLy5XAlN0CAyz12lHwDunW+PnIkSONFgqmNYbNIbPYfRyTCBQfARKY4uuIEsZAALEvy5YtkwEDBlzx1vr16xWpMS\/8NQ\/CopsMzpgx44r05yFDhijSYhaaS8NFFLUsVwJTFQuMbmSJGBg2h4zaLbxPBOqHAAlM\/XReqxVHWWB++9vfNjpAz507V0BgTp8+3QjIBYkI6vycNZCuBKZKFpjdu3c3YolAIFE3BriwOWTWu4\/jE4HiI0ACU3wdUcIECEQRGNMVgWkeeOAB6dWrl5pRpz5n5SKKWpYrgamiBcbEis0ho3YO7xOBeiBAAlMPPXOVEQgguwXupEcffVT2798vebiIopTiSmCqaoEx8TKtMWwOGbWTeJ8IVBMBEphq6pWriolAUAxMzGFSfdyVwFTdAmOCzOaQqW45DkYESoUACUyp1EVhs0KgSgSmDhYYrzUGgdZsDpnVp4PjEoFiIkACU0y9UKqcEagSgamTBcZrjWFzyJw\/OJyOCDQRARKYJoLPqYuDQJUITN0sMN5d5NccEs+wHUFxPm+UhAikgQAJTBoocozSI1AlAlNXC4zXrXTjjTeq2j1BBfDYHLL0H1suoOYIkMDUfANw+f+GQJUITN0tMOaeDiqAR2sMP\/lEoPwIkMCUX4dcQQoIVInA0ALTfUOYKddB1hg2h0zhQ8QhiEDOCJDA5Aw4pysmAlUiMLTA+O+xIGsMm0MW8zNJqYhAFAIkMFEI8X5mCCxfvlwmTZrkO\/6OHTtk8eLFmc3tHbhKBIYWmOBtE2SNGTdunKCY4fXXX6+6kM+fPz+3vceJiAARcEOABMYNN77liIDZbDGMpGhyc\/HiRVm1apVs2LDBcUa716pEYGiBida51xqDN9A88pNPPpGf\/vSnjf5Ya9eujR6MTxABItAUBEhgmgJ7PScFeVm0aJE8+eSTsnPnTisQ9Dvf\/OY3rZ53fahKBIYWGLtdYFpjOjs7pW\/fvoLmka2trTJt2jRBc09YY9Dd\/N1337UblE8RASKQGwIkMLlBzYmKjECVCAwtMPF2GppDtrW1yYgRI2TTpk1y6tQpNQD2BNxKaPjJlOt4mPJpIpAHAiQweaDMOXwRQKdoxBr07t37ivvnzp1Th4etpSYpxFUiMLTAxN8NIDFTp06Vo0ePysmTJ+W1115rDGKmXGO\/0hoTH1++QQSyQIAEJgtUOWYkAqNGjZIVK1bI2bNnZc6cOZHPZ\/1AlQhMUgsMXCvoLTRhwgRZsmSJCqbet2+fihGBi+Wjjz7KWh3i0sgyiVB6vrfeekvQjgBrPHHihBw+fFgNCysMYrFojUmCMt8lAukiQAKTLp4czRIBHcwLk\/3KlSst38rusSoRGBcLDEgLDnFYIvAvDvDhw4crAvPEE08ooonYkK1bt1aawICg4WJzyOw+axyZCKSFAAlMWkhynFgIaAtMe3s7CUwAci5WCGTRxLHAvDF0aOOwhhggLrA6wPqgSd2sWbOkX79+KuV9165d8vrrr2dOYlzWHmsDeh72mw\/\/x+aQSVDlu0QgWwRIYLLFl6OHILBgwQKZPn16rrEuQeLU1QIz8qtflR49eijCcvz48W7ExMTk\/PnzAqsZLDEtLS3d+gtlscmLQGD0utgcMgsNc0wikBwBEpjkGHIERwQYxBsOnMshHtcC879GjZL333\/fVxAvgdExMNq9AmvNq6++mok1xmXtjttQvRY1H1xsbA6ZBGG+SwTSR4AEJn1MOaIFAtqFdODAgVwr7tICI\/IzEekQkV+KyIlBg5RFxe8KIjAgLkEHuoXqrR6JIhRWg8R4yHY+NoeMASofJQIZI0ACkzHAHN4fAQbxRu8M20PVHCmuBWZva2ugIGEERr8UdKBHry5961OSOeNgzeaQSZDmu0QgPQRIYNLDkiPFQIBBvNFgxT1U8XyvXr3ko8tune+LyNuXp\/hh1zT4\/vNd36dhgTGlDzrQo1cY\/ESctSeZR7\/rMh+bQ6aBPMcgAu4IkMC4Y8c3EyKAGJj77rtPHn\/88dwK1gWJXNYgXl2zZeTIkWppCLaV7dvldhGZLCI\/EZHBIjJGRN4SkZdFBC0y4ULaLCJJLTAmnmlaY1wIRZLt6DqfjTWGzSGTaIbvEoFgBEhguDuagoDZ1NFPAFbiDQ4s9avZgiwiBOOeOXNGXOrA+OnAxoXkfS8ta4wroXDdzEnnCyJvGkNU792yZYuwOaSrhvgeEbgSARIY7goi0NX3BpVWZ8yYUZhS8eahCiXhZ11sDoG0ICyauJhKnHD6dFMsMGlaY5ISiribOo35TPIGvaCqLy5U72VzyLga4fNEIBoBEphojPhEDRAoqgsJcqE\/Dyrk4vKWuPeqBkG8Azs65HGLGJifi8iHCYN4w7ZG0IFus53SIBQ28+hn0pwPutLtCEBidJo6m0PG0QifJQLRCJDARGPEJzJCgHVg\/IHVLiIEOqMfEfoQvfPOO+ogDKrZokcCgbm7o0MF7EbFwLwiIn\/vmEYdZ0sEHehhY6RJKGxkTXs+k7xBZ2wOaaMFPkME4iFAAhMPLz6dEgJs5tgdSE1aTBcRyMiwYcNk27Zt1sXi8M60jg75MwsLDBwc\/5ihBcZcYdiB7rel0iYUUds2q\/n0uGwOGaUB3icC8REggYmPGd9IAQHWgfk3EHUWEawU+N486FwOVRAYBPHuv2yBaQ\/JQkL\/77EJCtm5boGgA907nsvaXWXCe1nOZ2aKwRoDtxL0jDkfeughueOOO1Tc1bJlywQZS7yIABGwQ4AExg4nPpUyAnWuAxOWRWS6iFwOVbOQna4BA9WZdWDwM2rBpJ1GbbtFgg50832XtdvO3yyLD9ZkNoeErtGeAe7Bhx9+WDXPXLNmDTOVkiiS79YKARKYWqm7WIvNqpmjN0X72LFjsnDhQjl06FAgAFkH8eLQ7tOnjwrGxZfOItKZRGkdqtoC4y1kZ46fdiE7113lPdDRBVtfVSQwem065bpnz56qjcOLL76o9sN3v\/tdmTt3Lq0xrhuK79UOARKY2qm8OAvOIojX22PJNtYmKwIT5iKK0oTLIZ53K4GoNdjc9+v27LJ2m7mCnsl7PuyLiRMnyvjx41WM08GDB5VoSLleunQprTFJlMl3a4MACUxtVF2shebZzHH58uUyevToUCtMmgTGz0UUVLMlTCsuh2qQBSbLVgJp7Cxvc0jtXtm9e7d1AHMSOVywTjIf3sWcU6ZMkY6ODjl79qyKjdEuRFpjkqLL9+uAAAlMHbRcwDXmGcSbF4HRGUS2LqIotbgcqmW0wJg4mO4VkJoXXnih0gQGMTBIsUYrCOj7yJEjol1ppjUGVXwR5MuLCBCB\/4MACQx3Q1MQyCuIVxMlZHcsXowuQP6XqwUmyEV0\/PjxxAevK4FpZiuBNDYTMG1ra1M1cEz3ShpjB43hgnVSebxzsjlkUkT5ft0QIIGpm8YLtN6smzlqkoQlpxnEm5aLKEoVLodq2S0wGpMw90oUbi73XbB2mcd8x29ONodMiirfrxMCJDB10naB1pp1M8c45AWw2FhgcLjgr2TtIrpw4YLqRYSvLC6XQ7WsMTBe\/PTag9wraePtgnVSGcLmZHPIpOjy\/TogQAJTBy3XbI22mUcmLEEEJksXUZRaXA7VKllgEB+ig3iDDvQoDG3vu2BtO3bQc1FzsjlkUoT5ftURIIGpuoZruL5169bJwIEDI91GCJJEB2oESO7du1d9j27Up0+fbnR9xiET1vk5S3ijDji\/uatmgTGzkJI0h4zSkwvWUWNG3bedk80ho5Dk\/boiQAJTV803Yd1wGy1atEiefPJJ2blzp5UE+p1vfvOb1s8jW2PAgAHdnj937pyqr+Gd10xXBaF54IEHpFevXupdnfqclYsoakG2B5w5TlUtMOYaXZpDZoF11JhR9+Pol80ho9Dk\/ToiQAJTR603cc1m7Mv69etl5cqVvtIg9XnSpEly8eJFZRnZsGFDZlKb6aqPPvqo7N+\/X9LIIkoqcJwDTs9VZQuMiWfc5pBRunDBOmrMqPsuc+p32BwyCl3erwMCJDB10HJB1whXz9ixaCl45bVjx47QtOe0l2QTxJv2nFHjuRxwdbDAmLgFHehR2Hrvu2Add4605gzrJcUCeEm1wvfLhAAJTJm0RVkzQ6BKBKbsdWCg5DiEwqY5ZNTGiTNf1Fi295POiffN5pB+BfDYHNJWG3yujAiQwJRRa5Q5dQSqRGBu7uiQ20Vksoj8REQGi8gYEXlLRF4WEZTz62hiN2ob5bkc7kEHelbz2Ywb9ozLGv3G8+slhedojUmqIb5fdARIYIquIcqXCwJVIjB1s8B4N0jQgZ4HmYizWdMiMJjT20uK1pg4muCzZUWABKasmqPcqSJQJQJTVwuMuSGCDvSgTZMmmbDdmFnMGVQvh9YYW63wuTIhQAJTJm1R1swQqBKBqbsFxtwktgXwsiATUZs1qzmD2hGMGzdOlRK4\/vrrVe0jNoeM0hDvFx0BEpiia4jy5YJAlQgMLTDdt0zQgW4+lRWZyMNthcwzXC0tLd2m85I33ER14yFDhsjDDz8s7777riIxaHTKiwiUEQESmDJqrUIy63ovWBLqwowYMUJGjx4dWUU3bQiqRGBogfHfHWHWmDISmI8\/\/lg+6eyUuy9dagRpP9\/SIj36928UYzTJ25kzZ+Saa65R7RlaW1uVNQY1kObPn6\/IDC8iUDYESGDKprEKyatL\/j\/zzDOKsGzatEk2b94sK1askLNnz8qcOXNyW22VCAwtMMHbJsgaUzYC09nZKYs6O+WHPkvF\/63o318F9uoL5A0uJPyBgM\/ZqVOnrBqY5vYB5EREwAEBEhgH0PhKcgR0RV78Mt21a5cyZeN7VOZdsGCBTJ8+3bf0f\/KZ\/UeoEoGhBSZ6lwS5V8zeS9GjJHvClTTB8vLQ2bO+5EVLBBLz04EDG5YY\/D9aMEydOlWOHj2qvtAvTPf\/ogUmmS75dnMQIIFpDu61n5UEJnoLuBxwdavEG42inTUGFo2+ffs2ul8nGdf2XRf9YuyPPvhAzl26FDnNwKuukt79+zee0\/MdOXJERo4cqawxDz74oGpgSgITCScfKCACJDAFVEpdREL8C+JdTBeStsYgsHDxYpRcy+eiBeZKnE1Mzp8\/rwJA87JQuB7uLrsFlom2trZu7hWXceK+47JGENRpHR3yNxaTfUtE\/rG19QoCAx3iuvPOO2XJkiUkMBZY8pFiIkACU0y91EYquItmz57dbb1590HC5FUiMIyBif\/xMd0rJ0+elNdeey3+IDHfcCUwD3R0hLqPtBhwI\/3XQYMa2Une+Yq452NCyMdrjgAJTM03QDOXP2rUKBWwe+DAAWVt8f6cp2xF\/GXuesAljYHRvYUmTJig\/kKHbvbt21dZCwz2mcb6rbfeUv2FvN2es9iLrvpNwwKD9RVxz2eBM8esLgIkMNXVbeFXprOQkIF06NAhJa+OjaELKV5DQ61s1xgYkBYcqLBE4F8ccMOHD1cE5oknnlBZYUi93bp1q7qX9eVyuCeRyZwP46BWCuJE3n\/\/fQGpyWLNLmuEfkFQ91ssFn3eT9ACY4EUHykrAiQwZdVcyeU2g3iReWRecCtNnjw511owRfxrNMkB930RefsyqDrNFt9\/vgvkB0Xkv4sIKn\/0+3f\/Th3U+MKFgxp9dE6cONH4C33WrFnSr18\/mTRpksoYe\/311zM50M094LL2JB8Jv\/mSNIe0kcVljSZB9Uuh1vP+zKdZJ11INlrhM2VCgASmTNqqkKxRBIZp1OlbYH4jIi+KyD1d3alfuRzI+ZKItI4fr2qGHD9+vBsx8QbxQmewxKDiKzJZdMPALLaly+GeRI6w+VyaQ9rI4rJGbYHxElRzPpAXdBv\/JS0wNmrgMyVGgASmxMors+hh8S46O8l0LWW91qpaYCaLSLuI\/F5EruuyyODnn4jIYBH5TJelZpungivwDspC0u4VWGteffXVTKwxLod7kj0SNV\/c5pA2skTN6TeGq4sQY9ECY6MVPlMmBEhgyqStiskKV9G9996rimlt2LBBrW7mzJmqtPlzzz2nitrldVWJwJhZSNNEZGsXiYHLQbsdQGw0iRkjItoac86o4BqWRp3FgW7q2uVwT7JXbOezbQ5pI4vtnOZYQRYY00VIC4wN+nymCgiQwFRBiyVeg18aNXoi5UlevNaGohT1cj3gBnZ0yHe63AgnReSUiPyPLhID64u+QGbCrDFwGelKrUF1YNI80MtAYCCjTXNIm4+kq35d0+RpgbHRCp8pEwIkMGXSFmXNDIGyW2B0FlGvXr3k9Kuvytwu0nKbiIDEoPDZ1MtBuoidAGmxscZ84StfkTVr1qhCZ2GF7NI60MtCYLScScmbK4FxTZMngcns1wcHbhICJDBNAp7TFguBshIYXbPheOYgAAAgAElEQVRFZxGBaPx++3b5ooj8TxEZLyJtl+Nc0BbzPztYY14bPlxOfPCBDB48OLIOTNIDvWwEJqk1xpXA0AJTrN8dlKZ5CJDANA\/72s+s41169+59BRbnzp1jM8dPf9qXNPjVbEHaM2qWnDlzRnDADRCRQZetLr+9nE77xyLyxOWfERsRxxqjA3x3ikjLsGFy8803R7YSSMsa43K4J\/lAJZ3Phby5zMkYmCRa5rtVQ4AEpmoaLcl6dBYS0nLnzIF9oLlX0S0wQAcHni42hwwgEBZNXDR65gEHV9FRETkiogqffcPSGhMU4Du0re2KVOsgrbkc6GW0wJgym+QNekEBvLDLlcDQAtPc3xWcvTgIkMAURxe1kiSsDkwzgCgqgYFcR48eVRVycUWVuPdLs\/2gK\/4F1hikTodZYzT2QQG+v2ppkZa+fVUga9QV90AvO4HR8kNXuh0BSAyIpt\/lSmAYAxO183i\/LgiQwNRF0wVbp7bAtLe3555x5AdFkQiMdhEBI\/QjQh+id955Rx2EQYehnwXGrMS7rit4FyQmyBoTFuD7v7ssOXBHhdWO8cPW9kCvCoHBOkzyBp35NYd0JTC0wBTslxnFaRoCJDBNg54TI4U674q7Qag3m8Bo0mK6iGBNGTZsmGzbts26WFxYoTNU4kVBu291xcZ4rTFB6dYaM5vaMUH42hzoVSIwei2apPhZzlwJDC0w\/N1JBP4NARIY7oSmIKBdSAMGINz0yqsuQbw6iwhWCnxvHnRZHXA64si0xqBPkl+Ar67kqzUUVTsGadxhV9iBXkUCo60xfs0hXfVLC0xTfmVx0gIiQAJTQKVQpPwRyNMCE5ZFZLqIsjzgUHn3\/xOReV3WmL4BAb42tWPMbKWzRiXfMGtMVLdnl7Un2TV5zIc5EBsD\/aOXFHR9ww03RGZ2metiFlISLfPdqiFAAlM1jXI9TghkTWBwaPXp00cF4+JLZxHpTCI\/oV0O1bgH3JG+fWVmZ6dvujWCfcNqxyRxKWG93gPdbA7psnYnxXe9lOd8OkOrZ8+eqjHmiy++mIqLEDlPL4vI4q4qzJtFZG9rawMWFrJLskP4bhERIIEpolZqIlOWdWDQEHLSpEkKyR07dsjixfi1HnxlRWDCXERRanY5VF2a\/X388cfycWenfO7SpSsCfINqx2jZk7qUMI5ft2eXtUfhGXY\/7\/mwLyZOnCjjx49XMU4HDx60Ej8uQT0xaJAiSZowmhafrPa81UL4EBFIAQESmBRA5BDxETC7UT\/99NOyYsUKQUbSrl27ZNmyZbJp0ybn7CQQo3nz5snq1auVYPp73TDST9o0f5n7uYj8arZEoeZyqCY54GAVWvjhh8oagwBftB\/wqx0Tpx2B2RwybL3e5pAu7pUoPItEYDShmDJlinR0dAjqIYWlXGvZXQiqfpcWmCQ7hO8WEQESmCJqpQYyeevArFu3ThVJg6UE2UmTJ0+WhQsXyqFDh2KjAevL6NGjG++bYwcNlgaB0RlEti6iqIW5EhjXIE\/IAxJzbxeJQbo1YmBss5XStMbAvQJS88ILL1i7V6LwLCKBgUUEKdZoBQF9IzbGdKV5ZU5CUElgkuwQvltEBEhgiqiVGsjkJTAgHQjsRFXepOnVICy4dIVf789pWmCCXEQgYyADSS5XAuOaZqtl1S6lxZcuCeJc\/pNROyYoW8mmOWQca0xbW5uqgRPHvZI31knm0xYY06VjU72YFpikqPP9KiFAAlMlbZZsLSaxMEkLuh+bFpS4y\/JaXExylIYFJi0XUdS6XAlMEguMKVNra6t868ABAWlBUTy4loKyldC2AF\/6SmqNwdrjulei8CyqBWb37t0NshvVS4oWmCRa5rtVQ4AEpmoaLdF6zDgYuI5APMaOHSsXL16UVatWSVjMStgysyIwOFzwV7J2EV24cEH1IsJXFpcrgUlqgdFrgVvtvvvuk\/v\/43+UGy9dkrEpNoeMSrfWa4\/jXkmiAxesk8znZ4ExxwuyxtACkxR1vl8lBEhgqqRNrkUhkKYLKUsXUZS6XA7VJAecVx4zLuj8+fMydOhQeXnXLt9spbTTrb1rt3GvROFZBguMKaNfLylaYJJome9WDQESmKpplOsRr8soLIh36dKlsnXrVoUarD5wX50+fbrR9RkHaVDn56yhdiUwaVpgNCYgMDpe44MPPuiWrZRFc0i\/tSdpDhmlKxeso8aMum87p9lL6s0335Qxhw\/L7Zdddt46PKwDE4U471cNARKYqmm0JOsJ60adNIg3Tho1DmhYGl5++WX17wMPPCC6HL5Ofc7KRRSlKtsDzhwnSwuMGXDqzVaKm24d1RwybO0uzSGzwDpqzKj7cfSryRv25pHt2wWp7GazTnyPWCVcqN3TISK\/FBHWgYnSAu+XGQESmDJrr8SyZ0lgAEucQnZ33XWXwBKD69FHH5X9+\/erlO6kWURJ1RPngNNzJXExeOX1upC8Ze\/NbCVdO8a0xiRpDgmSElZmP25zyChduGAdNWbUfZc50TuMFpgoZHm\/LgiQwNRF0wVZp0kswkQCidBp0HmIrkkMXEjvvvtuHlNGzuFywOVlgTGF97PGpNEcsm3iRFUHKIxIaoz8uj1HAmw84IJ1nPH9nnWZMwlBZR2YpBrj+0VDgASmaBqpiTxhFphmQJBGIbu05c77gItrgTGf91pjitYcMko3LlhHjRl132XOJASVBCZKI7xfNgRIYMqmsQrJa5IY3UIAJvJz584pl87OnTtzW22VCExadWCiXEh+yjGtMX8sIqZLqZnNIaM2kguZiBoz6r7LnLTARKHK+3VCgASmTtou2FrNdGe4lnBggrjAjaOr8uYlcpUITNZZSFGxQUVsDhm1j1zIRNSYUfdd5qQFJgpV3q8TAiQwddJ2gdZqWl82b96smjmioV0arQRcllklAtNMC4yJfZGaQ4b1F4LMLmTCZZ+Z77jMSQtMUtT5fpUQIIGpkjZLtBY\/9xFSmXUzx+nTp+fqRqoSgWm2BcZLYorQHBJkKqzbswuZSPpxc5kzrgXmD21tjeaQjIFJqjG+XzQESGCKppGayGO2ETh69Kjce++9jfYBNs0X04apSgSmKBYYraMiNIccM2ZMaLdnFzKRdA+6zBnXAjPyq1+VHj16KPKGy0xNL+KeT4op368XAiQw9dJ3oVaLgnPz58+X3r17q9orcB+BvAwcOFAWLlyoUmjzuor4yzyPA84sdObF2iWIN0xf2qXUjOaQkCusHYEL1kn3psuccS0wbwwdKpq8nTlzRq655hrRzSOLuOeTYsr364UACUy99M3VBiBQxF\/meRxwe1tbA\/dE2gQGE8Ea80lnZ+7NIfUig7o9u2Cd9MPkMmdcC4wmqCBv48aNkxEjRsimTZvk1KlTKmhet4ooSu2jpJjy\/XohQAJTL31ztTUgMEWKgQnbcGHZSmk3h\/TK4bXGeN0reXxQXAmMq4sQ1Y2nTp0qcNniC5ZOEpg8NM05skKABCYrZDnuFQgE1X3xgyrvWjBF\/GvU5YADlhNOn3Zq9pe1CynoI+HNVkrSHHJbS4v06N+\/0c8q7GNoWmM6Ozulb9++DfdKHh9fF\/26WmCwHj3fkSNHZOTIkcoa8+CDD6qyBbTA5KFxzpE2AiQwaSPK8UqJQFUIDA64gR0d8rhFs7+fi8iHObuQwkiMma0UtznkYBEZIyKviMhLIvJB376KkNhcsEy0tbV1c6\/YvJf0GVcC42qBMeeD7HfeeacsWbKEBCapIvl+0xAggWka9Jy4SAhUicDc3dGhOhNPFpGfiIg+3JGH8rKILO7qVozD\/u8HDZKWlhZfVWQRAxOmc9fmkD+8vNb2rqq\/+0XkKwaZ+VVLi\/S0sMiY7pWTJ0\/Ka6+9lvn2dCUwri5CplFnrlJOkDMCJDA5A87p\/g2BBQsWyOzZs6+AY\/369bJy5crcYaoSgZnW0SF\/ZmGBAaH5Rx8LDFwrqIQ8YcIE9Rc6avPs27cvtDt0mgqL0xwS5AVfL4jIv+\/6HmRGE7fPdOEQ5VrShzvSjZG1k7Q5pA0eLgQG47q6CElgbLTCZ8qEAAlMmbRVAVl1\/Zerr776ikJ1Okbmgw8+YBq1Y3VYHSMBS4R5kMO9Ylpg5ojIWBEx06hxwMESgS8c4MOHD1cE5oknnlBVkltbW2Xr1q2h3aGDtujQoUNV1su0adPUv7jWrFkja9eu9X0FtUtaPvlE7v\/DH+S3IhLUHBJrXHf5CxYnTWYwoGl90iQGnbXO9u8vIGjey+teAYFDnMj777+vaqhEtU9w+Wi6EJgkLkISGBct8Z0iI0ACU2TtVFA29DwaPXp0IEExC9zhL\/+8ripZYHSMBA50fb3ddcjrn38mIptFBHVCcFiDtOBgNy0PGpNZs2ZJv379ZNKkSYKmm6+\/\/rr1ga5Jy1133dVNlRc+OSl9Wj6jgkdRC8gbRKoP2xdffFG+3kVivM0h\/0ZEpnaN+j+6CBsIjb68LjQdH3POh8T4kQn8H6wxwAWBr1HtCOLuVVcC4+oiJIGJqyE+X3QESGCKrqEKyafJSXt7e6ibKIrkZAFJlQiMX4yEiRnIS8dlC8wvRURXaoWl4cSJE8rioC9vDAwsZLDEIGYm7ECHtQWWFpAWfB91bdmyRZYtW9btMfOwRUbax52d8rlLlwSWJTPA90ER6X855kWTme93ERk9mI6PiXIphZEJM+X61VdftSZvUet2JTCuLkISmCiN8H7ZECCBKZvGSiyvmUYdFueC+Bj2QnJrMBi3Uus748Yp4uJ3BQXxavcKrDXeA12\/Y7tNoywwumosxgtKtx50Oc4FLjFdO0aTmCCXkpmtpK0xUWQCVpgbb7wxVWtM1Jx+GCZ1EbKVgO3O5HNlQIAEpgxaqoiMJDDxFJnkgMMhDreRdiOZLiTTAuPaSiDoQN+4cWPD6qLdRH6rBnGB5QUxNX41SILW7g3w\/VbX4GbtmFMi4udSCrPGtA4bpoKWTcLkJ3dYO4J42k1OUG1chGalZVpg4mqIzxcdARKYomuoQvKRwMRTpiuBca0T4pXOJo3ae6AjVgYkJoi0oOM4SAv+DbvC1u5Nt8Y4ICxwLSGY96SPS8nGGjO0rU2OHz8e6SIKakcQT7vuBCaOi9AbpE0LTFwt8fkiI0ACU2TtVEw2Eph4Co1LYHCw9urVS\/ocPixZW2DMlXgPdBCMpUuXNrKNYGFBthEsLraXzdq1NQZjgrhEZSuZAb62sTFh8ia1xtis0Tt\/XBchLTC2O47PlREBEpgyaq2kMmsCM2DAgMgVsJWA\/V\/oIBA4TJFJdP78eZHt25vSSsDPGoMg3ihri99msD3ctTVm6qVLglgYkJigbKWo2BibdGuvrEmsMbZrNOdMo5UAu1FH\/vrhAyVBgASmJIqimNkiULYsJF1sDvV0cBDq9Oe3335bXCu1ehG2cSGleaCbY8U93EFkLnV2yjWXLgnosTdbKe\/mkGY2V9DOjbtGjEMLTLa\/Bzh6uRAggSmXvihtRgiUhcB4i81505+THHBpEBg9RjPcK3ruIjSHRGYXCuCFXa4ExpWgMog3o18eHLZpCJDANA16TlwkBIpMYNCXRxMXb7E5L4ZJXAxpEhiMZbpXbA70JBYYr+x+7QjaEqRbx20OqdsRgMQEWWNcCYxrkDYJTJF+41CWNBAggUkDRY5RegSKRmBw+I8YMUL+5E\/+RI4ePaoKyPkVm\/MjMK4HXNoERo+H2BybAz1NAoOxmtkc0iRv0Jtfc0hXAkMLTOl\/3XABKSFAApMSkBym3AgUhcDo2BYc+ug9BBLzT\/\/0T9Zl7ItkgTF3hM2BnjaB0eM1ozmknluTFL\/mkK4ExpWg0gJT7t9RlP5KBEhguCuIgIhK+V21apXMmDHDt7BaliDhcNcuIjMgF3OioWBUcTVTtqLEwAThFXagZ0Vg\/KwxWTeH9JI3v+aQrgSGFpgsP40cu0wIkMCUSVsll5Vp1N0VqK0tICn4C\/3ChQuqrL8u7Z\/3AZeVC8k7rrnuoG7PLmu3+XiY1pgsm0P6yYI1mc0hsXazsJyN\/EkIKi0wNgjzmTIhQAJTJm1R1swQyMsCE2Rt8asA63KIJzng8iIwpnslqNuzy9ptN4eOjcmyOWSYLDpDq2fPnqoxJjpug1jZXElchCQwNgjzmTIhQAJTJm1RVisEvJaeY8eOycKFC+XQoUOB72dNYMxiczisbAJyXQ7xJAdc3gRGz+fX7dll7Vabw3goy+aQUbJgP0ycOFHGjx8v27Ztk4MHD0a9ou4nIagkMFYQ86ESIUACUyJlVU3UmTNnyvz586V3795XLM21Eu+oUaNkxYoVcuDAAVm8eLHon5HFM2cOehb7X1kQmKBic4cPH7ZSpcshnuSAaxaBwbze5pAu7hUrUD0PZdUcskf\/\/qqtQ9gF\/U6ZMkU6OjpUlllYyrUeJwlBJYFx2SF8p8gIkMAUWTsVls0kGk8\/\/bQiHe3t7bJr1y5ZtmyZbNq0SVauXJkKAsuXL5fRo0eHWmHSJDBRxeZsF+VKYFyDPJtJYPTcpnsFpOaFF16wdq\/Y4up9Di6lln\/9V3nw3DnVigBXGs0hz\/Xvr4hZ0KX1ixRrxEHh5yNHjoRmnCUhqCQwrjuE7xUVARKYomqm4nJ5GzuuW7dOdQKG1WTBggUyefLkSLePLUR5EBgz\/Tmq2Jyt3K4ExjXNtggERltj2traZMKECbHcK7a4+j0HrPF10yuvqNt5NIf06temejEtMEm0zHerhgAJTNU0WpL1eAkMSAZSTeHmAYGZPn266mi8c+fORCvS86ChIMhR0OVigfELyLWJbbFdkCuBKbMFRmPj4l6xxTWIwCAj6KWXXpKzZ85IHs0h\/fQb1RySFpgkWua7VUOABKZqGi3RemB1weUlLajFEuXysVmmdlPh2TSDeP2sLYhr0enPNrLZPONKYIpqgQFJnDZtmtx1112qQzW+1q5d6wuFi3vFBtOgZ7xY59EcMky\/QdYYWmCSaJnvVg0BEpiqabRE6\/EG3ILQjB07Vi5evKiKym3YsCFyNd5A4PXr16vYmTjkBZNEWWCC0p9hcbHpPBy5EJ8HXAlMkSwwQ4cObZAWfI\/rwicnpU\/LZ9T3IDEI5PZeLu4VF4z1O2FYZ9UcMkq\/fr2kaIFJomW+WzUESGCqplGuxzrzCFDBGrBly5ZAAhNVbC5LuKMOOL+5kxxw3vFMUnf+\/HnromsgKtragn+9l0lgcA9B29CBeUW5V+I2h4zSUxTWWTSH\/OxnP2uFqdlL6s0335Q+hw\/L9y9bLt++vKgfdi0M3yNuB9fPRKRDRH4pIicGDVK1ZnAxiDdqF\/B+2RAggSmbxiokb1DKs06Bdl0qLDkDBw6MdBvhoN24caNqHQBLAMgM3FenT5\/2Le3vV2zOVUab96IO1SAC00wLjCY9NuvTz8ACA\/yjCIy+79IcMkoeG6yzaA75Qr9+8r979IhMudbWGJWavX273C4iky+nn\/9ERAaLyBgReQsWLRFBpBcIzGYR2dva2lg6CUzULuD9siFAAlM2jVVIXj+iYRt0GwRDULuCoLoy2sUxd+5cNeQDDzygDhPbYnNZqsPmUPXO32wLDAih6Sq68K8nZVCP\/+sKmDRp3Lp16xXkxc9a4B0gbnPIKD3FwbqZzSEHDBhAC0yUMnm\/NgiQwNRG1cVaqDcLyZQu7TRqm5XD+oKsp0cffVT2799v3f3ZZmzXZ+IcqnqOJEGeXjldXEjaqoWxvK4i\/B+ICwJ3vS4j79y2a9fP+XV7joO77Xx6TK81Jq\/mkEkIKi0wcXYEny0DAiQwZdBSBWWMIjBppVHbQhcVxGs7TprPxT1UMXeSAy4NAoMxvG4kkBYQlr179\/paW\/wwi7N2m+aQUXqJM585Vt7NIZMQVBKYqF3A+2VDgASmbBqriLze+BdzWTaF59KGoUoEppkxMFov2o0EXKOsLUkJjH4fB3RQc8io\/eJKYDBuns0hkxBUEpioXcD7ZUOABKZsGquQvHAV3Xvvvd1SpnVa9HPPPZdaKwEbyKpEYIpaB8ZGDyYZQWG53bt3x24l4NccMmruJARGj51Hc0haYKI0yft1QoAEpk7aLuBaQWJmz57dTTJdyyVPcatEYIpggUmqu6SEwtscMqqBZtL5TBJz74cfqp5K+y9nBn2r68YTXenNJ0XkVFevpfbLmUT40hdSovGzzixCpRykR29raRHdHJIWmKQ7i+9XCQESmCppk2txRqBKBKbuFhhzE9j0F8LzaREYjOUN8MX\/uTaH1CQGDTXOdnW4jkNQ5QtfkFOnQJlYB8b5lwNfLCwCJDCFVQ0FyxOBKhGYOAecWejMi7dLFlJaOkuTUET1F0qbwHitMfg5bnNIb40XtJh8SUQ+6NtX\/qSz07oOTL8\/\/3PlgnvrLVSJkW6F84q459PaPxynHgiQwNRDz4VYpZl5tGvXLlWBFXUt\/K6gui1ZLaSIv8xdDvEkLoaqEhi9rjBrjAvWNntRW2PiNIfU4wa5lP6niHUl3veHDFHtObC+M2fOyDXXXNOIKyrinrfBlM8QAY0ACQz3AhGw6IXUDJBcDtUkQZ5VJzBYX5A1xgXrOHvCtjkkWgQgDgbkRbcJMK0xcCm9JmJtgdGVeEHexo0bJyNGjJBNmzYptxIJTBwN8tkiIkACU0St1EAm1oGJVrLLoUoLTDSueMJrjcH\/uWY92c34f54Kaw4ZFeD7YxH5fQwLjOkiRAuGqVOnytGjR9UX2m2gaSraZ6BWDy8iUDYESGDKprGKyEsCE61IVwLDGJhobL3WmM7OTunbt69T2rbdbN2fcmlHgBFgjVnmYIHBu3o\/HTlyREaOHKmsMQ8++CAJjIsC+U4hECCBKYQa6iMEitRNmjQpcsEo5z9nzpzI59J6oIjmdFcCwyykeLsClom2trZu7pV4I7g93aNHD2n55BO5\/w9\/UGnXYe0IdLo1mjT+g6MFxtxPkPjOO++UJUuWkMC4qY9vFQABEpgCKKGOIoRZYJqBR5UIDC0w8XeQ6V45efKkvPYaIk2yvTShePHFF+XrXSTmj0XErBnzNyIytYuwgMSA0qPtKLtRZ6sbjl4OBEhgyqGnykmpWwm0t7fnWnE3CMgqERhaYOJ\/XDSZQLox2hEkbQ5pI4FpEUHW3cednfK5S5dUAbxviEhbF2H5z57id\/8pBQsM1lfEPW+DG58hAhoBEhjuhaYgQAtMNOyuLqSkFhjdHHHChAnKxbB48WLZt29fbkGuQMZl7dGIBj\/hda8MGTJExYm8\/\/77qoYKDvy0L781hgX4whpzBxp20gKTtio4XkkRIIEpqeKqIHYzmjbSAiOyWEQQS7FZRHSaLUgLDlS4UvAvDtLhw4crAvPEE0\/I2bNnpbW1VbZu3ZrJYe7VSzMJjCYrSZpD2nw+g9boF+CrrTEPikicOjBmFhKbOdpohc+UCQESmDJpq0KyagsMC9nZWQVsLQBx68Cg0Bk6R8PagAvzoG\/QiRMnGi6GWbNmSb9+\/VTwNQoQvv7665mTmCIQGK0Zl+aQNh\/VsDV62xEMNmJj4tSB2TNokLS0tChxSGBstMJnyoQACUyZtEVZM0OgiPEALod43DowutQ8CMvx48e7ERNvKwGQTlhicCAiFTeqQWISZbmsPcv54jaHtJHFZo1+1piXReTxrkaPutgdmj6iXQGun3VZ2H6O\/kkkMDaq4DMlRYAEpqSKo9jpIlAlAhMnBubC5z4nCCD1u4J6Ien4EByur776aibWGJvDPc0dYDufbXNIG9ls5zStMS+KyBtdmUjefknodgRyo12EKy5\/39na2hCFFhgbrfCZMiFAAlMmbVVQVrMuzPr161UtjtGjR8vChQvl0KFDua24SgQmaRYSXErTpk1TLiR8vfzyy\/LII490C+LNwiJhKtv2cE9rg8SZz6Y5pI1ccebEeCCMX\/vwQ3ldRL7YlaUUZIF5oavlwClaYGxUwWdKigAJTEkVVwWx161bp8qZP\/PMM4qwoEfL5s2bZcWKFcpVwUJ2n46d+RM3BkYHeYK0gKxo4oLS8iAue\/fulaVLl6pS87C2\/PjHP+5mcUnTIlEWAqPlTLr2uARGk5iJH34oaIE6SETQAADxMWNERFtgbhOR9+E+MoK08S4tMFX4rck1mAiQwHA\/NAUBv87UIDArV66UBQsWyPTp09XBuXPnzlzkK7sFRmcR9erVS6599VXrQmfvffGLirTMnYvyaKKIytq1a2XLli3qZ4x70003yV133SV33HGHIjXoIm72zknLIlE2AqPxQd0YkIO4cUEuBMYkqChspwkMmjye7HIhgcAgHuaXIsIspFx+fXCSJiFAAtMk4Os+LQlM9A6wOeB0zRadRXT+\/Hk5sn27oKsxAjvDgjxxwD39z\/+sgnJBWJAiDWLil1KNeiiDBw+Whx9+WGUtrVmzRhEd80pqkSgjgUlijbHRr3eXBAVpmwTGL02eFpjozxufKB8CJDDl01llJNZ1YEwXEtJ08Rc+\/tJHAbW8rjJZYPwIBrKIQDLOnDkjcWJgPvW1rymstTVBB+jiZzOl2tTDd7\/7XWWxAdmZP39+JtYYl8M9yV5JYz7TEgV9oABe2OUyp6uLkAQmye7gu0VFgASmqJqpiVxwF82ePbvbanfs2JErecHkRScw+gDSxeZALkBYNHHRAMY94HQdmKuvvrpRwM4vpdq7HWGFWbVqVWbWGJfDPclHJs35UAxQtyMAiYGe\/C6XOeOmyetChSQwSXYH3y0qAiQwRdUM5coVgaISGMh19OhRVSFXW0VAMMJqsEw4fdo6BkbXgfEjQzYKMK0x2nKm34trkTDnczncbeQNeibt+cy1A1u\/5pAuc8YlqIyBSbIr+G7RESCBKbqGKiSfWX13\/\/79uWYZRcFYJAKjXURoeIl+ROhD9M4776i\/5IP+mtdZRHj+7\/7bf5Pv7NoVGQPz\/wwYIJ8ZMUJZcZJcmBsB18DQLzbG1jKvjVMAACAASURBVCJRJQKj16JJil9zSFcCE8dFSAtMkp3Nd4uOAAlM0TVUQfmQPj127Fi1sosXLypXxIYNG5q60mYTGE1aTBcR\/toeNmyYbNu2LbBYnK7ZYmYRPfbYYzL5hRckrNDZv7S0yJr+\/QVZS2ldttaYIItEFQkM1mQGWpvNIV0JTJxChbTApLW7OU4RESCBKaJWaiKTtx9SM60yzSIw+nCDlQLfm3+pBx1wfjVbkEWEmi0IyMU4ffr0kf\/7ciFAvyykh1pa5OmrrpK+ffumvtOirDFhFomqEhjTGoPYGOgHKdcgMzfccIPs3r3bupoxY2BS37IcsMQIkMCUWHlVEn3mzJkqo6V3795NscrkSWDCsohMF5FJYNAJGjLCRYSaLLrQHFKfg7KI4Hba\/y\/\/IndfutQodPZ8S4t8qm9fdYhmeUVZY3S2k2mRqDqB0evT6eY9e\/ZUKewvvvhiLAJDC0yWO5djlwkBEpgyaasmsqZZyA7EaN68ebJ69epQN1XWBEZbRWBpwZfOItLBs36qBYH5yle+Itddd518+9vfVo+AuMDaomuw2JIh\/OWuuxLntY1MawxkRpCvl6SYFgkzMNnFvZJkXXnPB71NnDhRxo8fr1yEBw8etBKfFhgrmPhQTRAggamJoou+TNMCA1nRFwlVeZNcCIJFWwIUYIuKs8mKwIS5iILW5nURnT59WjZu3NjNRYQDNyqlOgl2ab4Li5FuR+DNVMI8ZgE83Rwyb0KR93xYN+acMmWKdHR0qNYZYSnXWh\/MQkpzZ3KssiNAAlN2DZZYfk0wEKiK69ixY6k2cdSWHIydpwUmqJKtt2aLqTozi8h0Eb3xxhty\/PhxZbFBkTp0JtZVd\/0yW4q6HUxrDFxecBeal7c5pEt8SJK1N4vAIAYGKdbQqU07AlpgkmiZ71YNARKYqmm0BOsxs5DSsrZ4l40A4UWLFkl7e7vqq5QHgdEWEVsXEWTWWUQgLfjez0WEeJFx48bJH\/3RHykSAwIQRoaKvAW0pcu7Ti2zGR8CUvPCCy9Yx4ckWXczCYwO4rVpxUALTBIt892qIUACUzWNFng9ZtZR2tYW77LRpgDXnj17Mo2BCXIRaauJnzpssoj8XEQ4vK655ppGFktYMbsCb4MGaUPqd1BzyLa2NhWwHCc+JMmai0BgIH9UY0xaYJJome9WDQESmKpptMDrAYFBQGrWNV8QT3PPPffII488IrfcckvqBMbVRWRaW3QWEVKfzc7PsLb4pVSbavWLGSmw2gNFgzUGsTF+zSFd4kOSYFAUAuO1RMFNaMbG0AKTRMt8t2oIkMBUTaM1W49f8O\/NN9+sLC8IAk4zCwnEBeRBu4guXLigXDlhlWxtXETa9YR\/dVxLWNVd\/Zf6jTfeWHprDNbi1xxSE4o48SFJtn7RCIzXGqObQ9ICk0TLfLdqCJDAVE2jNV+PtzieCYdfZhOsAHBjBGUhpeUiwhy6ZotLSnWQWm3iJsqwJbzNIZ977rluRd6yXmcRCYzWm9mK4c0335Q+hw\/L90UiW0X8UkRYibcMu58yuiJAAuOKHN8rBQJRFpjf\/va3isCAXMCdMWPGDEHastcqYtPsECRo2rRp3QrNubiI4gIbFTcRd7xmPm9aY55\/\/nn527\/920YQb5LmkFFrKjKBMa0xqvXD9u3WzTrZCylK87xfZgRIYMqsPcoeiUAUgTHjMDDYAw880OgPpEmLi4sIhAhxLi7xMpGLCnggayuFq1xx3gNeN910k\/zlX\/6lXH\/99ak1h4ySoegERss\/YMAAWmCilMn7tUGABKY2quZCwxCA9QWpzI8++qigJ1PcLCLTRYR5XFKq09BQWa0xfnhNmjRJkKkEIugtgGeu06Y5ZBS2ZSEwWMeE06dpgYlSKO\/XAgESmFqomYuMQsCmEi\/iNODi0IXmcLDC0hKVRRRGhqLkcr1fBmuMTZXitJpDRuFYFgKDIN6BHR3yuEUMzM9F5MPW1sbSvWu02fNRuPE+EWgmAiQwzUSfcxcGgaBf5kFZRM1wEcUFK8uYkbiy6OddXWpJm0NGyVsmAnN3R4d8XkQmX64b8xMRGSzSaNb5sogsFpEOEXnlsrXm7wcNavTAIoGJ2gW8XzYESGDKpjHKmwkCJoHBBDogF\/8GdX6Om1KdieAWg5pZLDb9diyGjP2In7UFhfjC4ou8kyRpDhklcJkIzLSODvkzCwvMWyLyj7TARKme90uMAAlMiZVH0dNDQBMYuIPgIsKl41qK6CKKu\/K0Y0Zs5rftlG0zlvmMS3PIqDnKRGCu7eiQ\/ZctMO0hFpg5IjKWadRRauf9kiNAAlNyBVL8dBDw1iEpg4vIZeX6oM6yEaS2tqBBIeaxKfgXdy1xm0NGtV0oE4G5uaNDBfH+0ADtbRHlVtLXz0Rks4j8oa1N9NrpQoq7y\/h80REggSm6hihfbgh4Gw3u2LFDVd3FZZNSnZugCScyCQbWBbcSiEaSK8jaknUAs21zSG9Jfu9ay0RgYIHxFrIz1wPyghgYFLIb+dWvSo8ePZSOcaH7tW4eySDeJDue7xYBARKYImiBMhQGAbMuzMGDB+Wv\/\/qvZd++fYkP+MIs0BAEh\/aYMWMStSPwtlewKfiXNhY60DqsOSTWifUeOXKkYZEw5SgTgdEWmKggXlhg3hg6VOkY60MnczQDJYFJewdyvGYhQALTLOQ5b6EQAHFBposZtIvv\/RoNFkrwFISJ2xxSW3CuvvpqdTBm6Y6Ks7yw5pAYJyy1vEwExs8CY7qQTAuMbiWAtY8bN05GjBghmzZtklOnTgW2z4iDOZ8lAs1EgASmmehz7sIgoLOOdL8iLZhfo8HCCJ2iICAlUc0h\/YrNIYsIVpciXWE6Cyr0VyYCE8cCY7YSgDt06tSpcvToUfU1cOBAWbVqlWqfgUw7XkSgbAiQwJRNY5Q3dwS8Ab5r167NXYa8JvRaKRCAO2TIEBULhMO\/KNaWKDyidOZdJ8Yz40Oixk\/jvgtpQiE7FwsM5NXzwY2GAGtYYx588EESmDSUyTGaggAJTFNg56RlRCCsmFoZ1xMkMw66trY25XI5d+6cihtpRmxLGphGFcDT8SGdnZ3St2\/fRnxIGnNHjeFKYFwtMOZ8kO3OO++UJUuWkMBEKYr3C4sACUxhVUPBiohAVGn7IspsK5O32BwOvJaWFnnvvfdUFkvRXEW264rSGaxLIGxmfIjt2EmecyUwSS0wDOJNojW+WyQESGCKpA3KkhkCUV2p405cFWtMVLG5sjaH9NNnmM7M+JCTJ0\/Ka6+9FndLxH7elcCkYYGBK5Bp1LFVxhcKhgAJTMEUQnHSR2DUqFGyYsUKGTx4sApa3LBhQyqTRP1ln8okGQ0St9hcGZpD2kAVpDNNJmBpglspj1gfVwJDC4yNpvlMHRAggamDlmu+xgULFsj06dMVCqtXr06NwGhYy2KNSVpsrojNIV23tldnqPmjg3gxJgKXEeiaVqE\/PzldCQwtMK5a53tVQ4AEpmoa5Xq6IXDbbbfJokWLpL29XZGYLAgMJoxqNNhMtaRdbK4IzSHTwNPU2fbt2+Uf\/uEfugXxplHoL0xOVwJDC0wa2ucYVUCABKYKWuQaAhFYvny5urdnzx6ZN29eZgRGCxDVaDAvVWVdbK4ZzSGzwk7r7PTp0\/KDH\/xANfE0r7iF\/mzldCUwtMDYIsznqo4ACUzVNVzj9SFw95577pFHHnlEbrnlllwIjNcag8Nw\/vz5uWkh72JzeTSHzBI8TfRQxG\/WrFly\/fXXKwLj1ZlNob+4croSGFpg4iLN56uKAAlMVTXLdcm6deuU5WXlypWSdhaSDbxRjQZtxrB5xpv+nEcAqilXFs0hbdbt+oyOBdJkD3jpOjcgMAj0RmXaLVu2iLdoYZrBzK4EhhYYV83zvaohQAJTNY1yPQoBxL4sW7ZMBgwYcAUi69evV6Qmjyuq0aCrDH4Buc0uNpd1zIgrVvo9k2jh\/0BcDh8+LGiHYF5ROksrtdyVwNACk3Qn8P2qIEACUxVNch2hCDTDAmMKFNVo0FZ9ftYWv0PYdrwsnssqZsRF1qg6N2FjRuksqTXGlcDQAuOyE\/hOFREggamiVrmmKxBoNoHRArk0h0xyCDdrK2QRMxJnLX6xQNpCFWccPOvSHNJmDlcCQwuMDbp8pg4IkMDUQctcY6EQiGo0qIWNW2yuUIvsEiaplSLOmrKMBYrSmcs6XQkMLTBxdgWfrTICJDBV1i7XVmgE\/ArgJS02V8QFpxUz4re2vGOBbJtDIq4GVX3DLlcCQwtMEXc5ZWoGAiQwzUCdcxKBLgTMYmrPPvusvPLKKyq4tNkBuVkoyMVKESSHNyBX4+UNyM1iHVEtJGwL\/bkSGFpgstAqxywjAiQwZdQaZa4EAjgI8Rc9CqmhiFpra6v616+YWiUWLCJJrDFFs07ZWmNArvyaQ7oSGFpgqvJp4DqSIkACkxRBvk8EHBDw1ojZunWrGmXp0qWqS\/CaNWuuqEHiME1hX7G1xoTVbAExaPYVZY0JK\/TnSmBogWm21jl\/URAggSmKJihH7RAAUfGWrQcIZWkOmVRhYc0hswzITSq33\/tR1hi\/5pCuBIYWmCw0yDHLiAAJTBm1Rpkrj0CRm0OmDb6OGfnkk09UUbm+ffsKDveyxQJF6cxb6A8WJN0BG2u1uYARLTA2SPGZOiBAAlMHLXONpUWgKM0hswRQW2ImTJggvXr1UiRm165dV1TIzVKGNMeO0pl2n\/Xs2VNaWlrkxRdfVGTN5gKBoQXGBik+UwcESGDqoGWusdQImH\/Z590cMivggrplX7hwQcaMGaMOdBAZVBku4xWlM6x\/4sSJMn78eNm2bZscPHjQapm0wFjBxIdqggAJTE0UzWVmj4C3\/9KxY8dk4cKFcujQoVQmz6s5ZCrC+gxiG5BbtuaQYXiF6QwupSlTpkhHR4ecPXtW1Y2JCkymBSar3clxy4gACUwZtUaZC4fAqFGjZMWKFXLgwAFZvHix6J9xMM2ZMyc1eaMaDaY2UYoD2TZR9E5Z9OaQthAF6UwH8SLFeuTIkSru58iRI6FWJ1pgbFHnc3VAgASmDlrmGpuCwPLly2X06NGpWmH0QqIaDTZlwcakafZvKlJzyCS4enX23HPPdQvitUktpwUmiQb4btUQIIGpmka5nsIgkCWB0Yt0aQ6ZJUBpNlE05Wx2c8g0MTN19vOf\/1zFwOgg3qhCf7TApKkJjlV2BEhgyq5Byl9IBHQ8DIJu4VLK8opqNJjl3Bg7z5otNlaKrNebdHzgddNNN8ljjz2mqi\/7FS0MWictMEnR5\/tVQoAEpkra5FoKgYCOf4EwaQbxRi0uzwJ4eTdR9FpjkKlkEzMShVme9\/2sU5MmTZK5c+fKu+++K8uWLetW2NCv0B8tMHlqjHMVHQESmKJriPKVCoFmkRcNUlRp+6RgNrOJolf2MlhjbKxTUTozm0O++eab0ufwYfm+iLx9GZAfdoGC7z\/f9f3PRKRDRH4pIicGDVK1ZnB5K\/\/qDKkZM2YoAsWLCJQNARKYsmmM8hYWgawyj1wWnKY1pmhNFItujXG1TkW1I4DVCYX+ZPt2uV1EJl923\/1ERAaLyBgReUtEXhYROCxBYDaLyN7W1gZcJDAunyS+U2QESGCKrB3KVioE1q1bJwMHDszVbRQGUNRf9mHv2tZsKYqCimCN8bO2oBAfCvLZXlE6GzBgAC0wtmDyucojQAJTeRVzgXkg4C1ip+c8d+6c6jC9c+fOPMTwnSOONcbG5dG0hURMHNYcMiuZ00wXN2UM0hljYLLSJMctIwIkMGXUGmUmAjERCGs06OryiClCbo+bMSM21W1dBDNjgZACjRYIsLTEsbZEzeunM2YhRaHG+3VCgASmTtrmWmuPgNlo8K\/+6q\/k1KlTggNfd35GKfs0D+FmAW5aY7AmVLtNejUrFsjU2Y9+9CNp+fWvGQOTVJl8vxIIkMBUQo1cBBGwQwB\/1cM9gUMRF5oIPvLII3L8+HHrjsh2MxXjKR24mqQ5JIgLYmxMogeSF9W3KE0EtDVmyJAhcv+UKcxCShNcjlVaBEhgSqs6Ck4E4iGAWBwQF6TMosDe3r17VXwOft6yZYusXbs23oAledqlOWRQt+xmd8dGATxaYEqy8Shm5giQwGQOMScgAsVAAJYXXCZRKWNzSFc0bZpD+hWby9vaEra+uDEw8oUvKDchLqZRu+4cvldUBEhgiqoZypU6Amam0I4dOzIv8Z\/6AjIcsOjNIdNcurc5JMaGawYuIlhekrib0pTTb6y4WUj9\/vzP1XoQzIzrhhtukN27d6v\/YyG7rLXF8bNGgAQma4Q5fmEQQHPFo0ePyubNm1Ufmqeeeqqp6c2FAcYQpGjNIbPCCNYIEFoQF6S6HzlyRMW0FMnaEkRgru3osI6BeX\/IEBk7dqyyvpw5c0auueYaEpisNhXHzR0BEpjcIeeESRDQ1W7b29tl5cqVjaFQRA6\/qHGtX7++2z3vfDi4vve97wkyOg4dOpREnEq+2+zmkFmC6q1zg4Mdpfbfe+89ZaXIMzDXZZ1xLTC6Ei+sTuPGjZMRI0bIpk2blFuJFhgXDfCdIiFAAlMkbVCWUAQ0eRk2bFg3krJgwQKZPn26Cki99dZbG9\/7FY+bOXOmzJ8\/X1566aVKuJBgVUJDQFxpu8XiFMAr8taNKjZnplzDEtPsQN0wLOPGwJi9kOAimzp1qrJC4gtVo1etWiXshVTk3UvZwhAggeH+KAUCmnjgL8err75a\/RWpLTCwvuCaM2eOmBYapAaDrPTu3fsKq4x2J5lWnFIAYQgJTObNmyerV69W\/6u\/37BhQ2pLiSptn9pEGQwUt9hcEdoRRMEQ1wLzpYcfbgRt6yBekLSRI0cqa8yDDz5IAhMFOu8XFgESmMKqpl6C4TCeMmWK3H\/\/\/Wrh+Pm+++6Txx9\/XMWp4C9HXB0dHbJs2bIGgdGE5cCBA8qi4v3ZRBFEZ8+ePYr4mN+XFWmQsNGjRzd6L2FNIG3AIe2rLNaYpMXmim6NiWuBefqf\/1m5yPCZQc0fHcSL\/XHnnXfKkiVLSGDS\/rBwvNwQIIHJDWpOFIUADmAQkH379sn48eN9ewjpTCJtgfGLiQk6yM0spP379yuLTZkv0\/KEdXh\/TnttRbbGpF1srqjWmLgWmPe++EX1OUK8y7PPPiuvvPIKg3jT\/mBwvKYhQALTNOg5sR8CXquC95kkBKZqiHuJGrBDVk3WxKwo1pisi801ozlk1B6Na4HRMTBaZ6dPn5Yf\/OAHqpAhg3ij0Ob9oiNAAlN0DdVIviQWGBsXUtWgbBaBAY5hzSGzxjnvYnN5NIe0xSyuBQZZSJro3XjjjTJr1iy5\/vrrVeXlrVu3MojXFng+V0gESGAKqZb6CRUVA6MR8VpgtOtEx34EpVlXEdG8XUh+GJqNBhFngb\/ss7i86c95F5vLojmkC05xLTCfvuUW+exnP9to1ok6N7fffnujhQSIKLOQXDTBd4qAAAlMEbRAGawR8CMwcdKorScqwYNel1GWQbxhcJjWGBAYZH6lcfkF5Da72FwazSGTYBPXAtPjjjtUWrjZYdxs6Ik+WCAwvIhAGREggSmj1mossx+B0VYY20J2VYEvjzTqOFjpmIqkzSH9rC3eQziOXGk\/69IcMi0Z4lpgdAwMSAv0M23aNPWvbugJN1JWVrO01sxxiEAQAiQw3BtEoMQIZFnIzgUW1+aQUcXmXGTJ+h2b5pBpyxDXAoMsJJCWuXPnKlFAXNDMEzEwvIhA2REggSm7Bik\/ESggArbNIeMWmyvgUsXbHBLxOVldcS0wug6MDtoFgeFFBKqCAAlMVTTJdRCBAiLg1xwyabG5Ai5TZfogywf\/pt2OwEwX79Wrl8j27XK7iEwWkZ+IyGARGSMi6DeNEGqUMewQkc0i8qmvfY0uoiJuGMqUCgIkMKnAyEGIABEIQgBupY0bN6rbupgarBTNDsjNQmNpFcDTJE+njGu8YEHpc\/iwdTdqsxdSFuvlmESgmQiQwDQTfc5NBCqMgI6HQao1vkcRtdbWVvWvLqZWxeUnaUdgutSADYiLN4D5P\/TtK7PuuUc+HjFC\/t+NG2Xkrl2BFhjdjbqKOHNNRIAEhnuACBCBTBDQJeyR5aKzXYrcjiBtEGytMTYBzDqLaMKECQJC2Oudd2T\/rl0yaeNG6bdrl7x9uZXE57sW8LMuF9KKnj3lqquvTntZHI8IFAYBEpjCqIKCEIFqIYBDNyhotCjtCLJGPMgaE+QiQr0WuNb0ZWZ14f+A5zPPPCNj1q8PjYH5l5YWWdO\/v6iYGV5EoKIIkMBUVLFcFhEoOgJ1tMYgiwgE5ZprrlHq8aso7FezBVlEe\/fubQTk4r2FH34oP+xSsmmBeailRZ6+6irp27dv0bcA5SMCiRAggUkEH18mAkQgKQJVt8ZoawtcSuiyDqvIwYMHVWdobW3xuohsCs19\/PHH8klnp9x96VIjBub5lhb5VN++KhuKFxGoOgIkMFXXMNdHBEqAQDObQ2YBT5CLSBOWMWPGKOvLuXPn5E\/\/9E+7FZqDtQXF5uJcsOy0tLTEeYXPEoHSI0ACU3oVcgFEID8EdCuHAQMGqEmPHTsmCxculEOHDqUiRF7NIVMR1mcQm6aTIGtf\/vKXVQ8idIbWrRdMF1FW8nFcIlAlBEhgqqRNroUIZIiA7vR94MABWbx4seifz549K3PmzElt5qyaQ6YmoGcglywi7SLSrqM0m2BmtU6OSwSKhgAJTNE0QnmIQIkQQC+m0aNHp2qF0ctPqzlkFnCCtPTp00euvfZa9aULzenifHpOby0cv0aX5jrRyZvl\/rPQGMesIgIkMFXUKtdEBHJCIEsCgyW4NofMavm2LiJv52dvFpFXPr3OuLEvWa2T4xKBMiBAAlMGLVFGIlBABHQ8DNwfcClledk2h8xCBj8XkbcNAgiIaW3RLiLEtbDzcxZa4ZhEQIQEhruACBCB2Ajo+Be8mGYQb5Qgfs0ho95xva\/7ECV1EbnOz\/eIABEIR4AEhjuECBABXwRmzpwpiMno3bu3ur9+\/XpZuXJlI3g3b\/KihYSlY9WqVcrisWbNmtgpx2HqDnIRHT9+XMW54PIrNGe2S+B2IgJEIB8ESGDywZmzEIFKIJBV5pELOGkVwLNxEUE+HdeCVG+6iFw0xneIQLoIkMCkiydHIwKVRmDdunUycODAXN1GYYAmaUcA4oLquNpFdOHCBUEvInyZ1p5p06apBoq6txNiWtCcMs9sIQRLT5o0qSGXtoZVerNxcUQgAgESGG4RIkAErBDwFrHTL6GaLDpP79y502qcLB6ytcaU0UW0YMECmT59egNj\/HzvvfcqN9qGDRuygJNjEoFSIEACUwo1UUgiQASiEAiyxti6iPA+iJB2EcHCAktL0bKINJHctGmTikniRQTqigAJTF01z3UTgYoiYFpjnnnmGUGlYFw69bmILqI4qiCBiYMWn60yAiQwVdYu10YEaoaAN0MIy\/\/Nb36j3C9VySJCPIyui9NMt13NthaXW0AESGAKqBSKRASIQHwEzGJ3umT\/e++9p8gLfkaV2wkTJigXES6d+lw0F1HYyhH\/Mnv27EZKe3yU+AYRqA4CJDDV0SVXQgRqjwDcR96uzjo2BuDg+2ZkEaWhGJKXNFDkGFVCgASmStrkWogAEQhEAJYYWF\/SLn6XB+TMPMoDZc5RNgRIYMqmMcpLBIiAMwK26dbOE2Twoq6I\/Nxzz9Um60hbmy5evMh08Qz2VFWGJIGpiia5DiJABKwQ0OnSiInJsxidlXA+D3mL2OlHduzYkXkTTVeZk7yHas+PPfaYPPXUU3LdddfJlClT5P77708yJN+tKAIkMBVVLJdFBKqGACwR8+bNk9WrV7OAWwmUC33dd9998vjjjzeKHJrFEG0KIMISM2LEiEoStRKosPAiksAUXkUUkAgQAd2DafDgwXQplGA7aLcXXEBmlWa0osA1Z84cMb\/3WxIsTxMnTqS+S6DvZolIAtMs5DkvESAC1gjocvp4gRYYa9ia8qB2ee3fv1+GDx\/eIDDeAnymhebWW29V6eHemBfTncSaN01RZ6EnJYEptHooHBEgAjj4Fi1aJO3t7aonEAlMc\/eE162DnydPntxo8Pmd73xHfv3rX8vXv\/71bj2cvC7AIJeg1veTTz6pFgrd43sSmObqvYizk8AUUSuUiQgQgQYC+Ise1549e3KJgYkbp1E3VWl33tVXXy0HDx6Uz3zmM77dyb1NKL0xMcD50UcflV\/84hdXxDTpLCRgy87bddth9uslgbHHik8SASKQMwI49O655x555JFH5JZbbsmFwMSJ08gZjkJNFxXDkoTAFGqhFKawCJDAFFY1FIwI1AsBHfjZu3dvtXD85X3zzTcrywu6LueRhRQWp0EXxr\/txyQWGDOLLA991usTVL\/VksDUT+dcMREoBQKmK8crcFZuBds4jVIAmJGQUTEwelqvBcbrMvJLs85IZA5bUQRIYCqqWC6LCFQNgTz+Yo8Tp1E1fNNej5fAYHy659JGud7jkcDUW\/9cPREoDQIkMKVRlRLUj8AwQLpcOiy6tCQwRdcQ5SMCRCA3BOhCyg1qTkQEEiNAApMYQg5ABIhAVRBgnEZVNMl11AEBEpg6aJlrJAJEwBoBxmlYQ8UHiUBTESCBaSr8nJwIEIGiIcA4jaJphPIQAX8ESGC4M4gAESACRIAIEIHSIUACUzqVUWAiQASIABEgAkSABIZ7gAgQASJABIgAESgdAiQwpVMZBSYCzUMAjRVHjx5t1bwvrpTeMv5x34\/7POa77rrrrmgk6DdOlGxm80Hz\/XPnzsnSpUsbnZSjxvHOjbTu3\/\/+94XuxOztRh1XD3GeZ\/XeOGhV\/1kSmOrrmCskAqkhkCWBSU1Ii4HiEomo5\/2KtkEM4DVx4kRZtWqVFVEyRQ8a02J5uT0S7pu3DwAABNJJREFU1lE6KyGA6ZAhQ2TOnDlZTcFxS4IACUxJFEUxiUARECCB2aQaS3qvMLKBtOyBAwf6Wq3CdFoGAtMMMtEM0lSEzx5luBIBEhjuCiJABKwRiEtgcHiPHTtWjX\/x4sWGJUIfQgcPHpRbbrlF4Gp59tln5dvf\/rZs2rRJdu3aJcuWLZMBAwZ0k80cw9u9ev\/+\/Y2\/ynVF3d\/85jcyZcoUMTtce8fW73nHw8T6nqsFBmOY1X3hDsK6sEYQIW\/DSu1yuvXWW2X27NmNtevmlX6uqh07dsjixYsb8\/itWZMu3Ul62LBhV+gE\/+EdX4\/tt0H8MNH7A3r6\/Oc\/r14Dhk899VQ3fer1xH1ey9EM4mT9IeGDuSFAApMb1JyICJQfAVsCc+LECVmxYoVa8MKFC+XQoUPd3Cn6IP\/ggw8a98NIgj54z549q0iKPmj1Qajffeedd9R9TUZOnTrVbf4vfelLKiYFlx+R0MTCPMwxhyY95n1Tm2HWEnNd5jibN29WGB04cEAREFxmET3vmH69oEwXFd6fP3++BK05SCcakxkzZoj+fufOnQ1y9fLLLzfki1oz5Jk0aZKYpAoygdDoWCA8Y84Z53nIpXUzffr0bvFF5f90cQVxESCBiYsYnycCNUZAH1BBEGgLAoJj582bJ6tXr27EfmgSggN748aNikCYh2MYgfG6YcyDXstiHviYHwfnc88913D5hFlCoiwMWRAYPWYcguCV01yTJjBBa8Z9r070eFo37e3t3VxkYcTMj8ya5ARkw9S5JmnmmF7SFPW8JjB5NPas8ce8NEsngSmNqigoEWg+ArYWGLhA\/P5C1qZ\/7VIwLRpBBMbmUAQyZmyE32FtQ2C8LhaMk5UFBm4dkxB6M5aCyIPX7aTdalFrRuBrkNXCO6a507xy6XtBBMbMUosiJCAwcZ7XBCbKpdf8TwolyAMBEpg8UOYcRKAiCCQlMNpyYktgcIjfe++93bJ4\/A5FwOu1RnitDbaxKH5xL0ksMGbqL+Q0XVd6W\/i1L\/CSQDNGxxv3AktXEgKjxzatN1FblgQmCiHezxoBEpisEeb4RKBCCNgSGFsXUpgFJuxQtXEhxSEwftYOc\/4kBMaUNcpyYN7HtjEtJn7YmwQvisD43fcSqCB3lt8WbiaBoQupQr9UEiyFBCYBeHyVCNQNAVsCYxvEG0RgNGHQQblenG2CeMMIzO9+97tuAbReS49pEUniQvLWgfEL6DVJgzeOx1yD15WmCRawQZ2ZKAKj16wDob1Wq5tvvvmKmjVhKeB+ReW8+yMrF1IZUszr9ruhGeslgWkG6pyTCJQUAVsCo2MVwtKova4U83AfMWKEymbxu3TmkU0atRlE7P2rXcefHDt2TGUq3X\/\/\/d3mfP755wVuHBz4fi4vUzbXSrzeNZhp4mY8jnYZefGEywdWGpCgPXv2XBGk611zVBq1N0hbY4MsMu\/lF\/ibF4Hxs8CV9CNFsRMgQAKTADy+SgSIABGoMwLNqMfCQnZ13nHd104Cw71ABIgAESACTgg0g0w0gzQ5gcOXMkeABCZziDkBESACRKC6CLCZY3V1W\/SVkcAUXUOUjwgQASJABIgAEbgCARIYbgoiQASIABEgAkSgdAiQwJROZRSYCBABIkAEiAAR+P8BPnL1L9xY7kgAAAAASUVORK5CYII=","height":337,"width":560}}
%---
