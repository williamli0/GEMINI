%% GEMINI Parameter Optimizer
%% Housekeeping
close all
clearvars

%% Constants 
E = 100; % Beam energy (eV). Science: 100eV to 1keV (Like IMAP-Lo)
% Actually want 1.2 keV for buffer/expanded capabilities
V = linspace(100,10000); % 2*Deflector Plate Voltage
qArray = [-1 0 1 2]; % Array of possible particle charges
q = 1; % Particle charge (e)

% 12.6 degrees peak scattering angle
Lmax = 19/tand(12.6); % Maximum L so all vertical scattering lands on MCP
r = 20; % Radius of MCP (mm) FIXED
s = 67; % Distance from CF to rotational stage center (mm) IDEALLY FIXED
d = 14; % Distance between deflector plates (mm). Want to maximize d
% If using the AST250 Equations, change t = 4 and w = 1
t = 0.05; % Thickness of slit (mm)
w = 1; % Width of slit (mm)
l1 = 5.588; % Distance from slit plate to deflector plate (mm)
l3 = 5.588; % Distance from deflector plate to MCP front (mm)
% 5.588mm min safe distance for high voltage (with 10% safety margin)
l2 = Lmax-s-l1-l3; % Length of deflector plates (mm) based on Lmax,s,l1,l3
alpha = w/t; % Unitless. For ease of calculation later

% Arrays used for optimization of V and E

EArr = linspace(100, 2500, 25); % Array of beam energies used for GEMINI

% Creating 2D arrays of voltages at different energies and vice versa
[meshV, meshE] = meshgrid(V, EArr);

% Expanding constant r and d/2 values to properly sized arrays
meshr = ones(size(meshV))*r;
meshd2 = ones(size(meshV))*(d/2);


%% Graph of the constraining equations in terms of deflector plate voltage
% Each equation has a RHS and LHS. The intersection point is a Vmax or Vmin
% Each side of the equation is in length. Units: mm
figure(1)
y1 = yline(r,'--','Color',[0, 0.4470, 0.7410],'LineWidth',1.25);
hold on
%eq1 = alpha*(l1+l2+l3) + (V*q*l2/2/E/d)*(1+alpha^2)*(l2/2+l3); % AST250 Eqn
eq1 = (V*q*l2/2/E/d)*(l2/2+l3);
plot(V,eq1,'Color',[0, 0.4470, 0.7410],'LineWidth',2);
y2 = yline(d/2,'--','Color',[0.8500, 0.3250, 0.0980],'LineWidth',1.25);
%eq2 = alpha*(l1+l2)+(V*q*l2^2/4/E/d)*(1+alpha^2); % AST 250 Eqn
eq2 = (V*q*l2^2/4/E/d);
plot(V,eq2,'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
%eq3L = (V*abs(q)*l2/2/E/d)*(1+alpha^2)*(l2/2+l3); % AST 250 Eqn
%eq3R = 2*alpha*(l1+l2+l3); % AST 250 Eqn
eq3L = (V*abs(q)*l2/2/E/d)*(l2/2+l3);
eq3R = 8;
plot(V,eq3L,'--','Color',[0.9290, 0.6940, 0.1250],'LineWidth',1.25);
y3 = yline(eq3R,'Color',[0.9290, 0.6940, 0.1250],'LineWidth',2);
hold off
xlabel('V')
ylabel('Length scale (mm)')
ylim([0 30])
xlim([0 10000])
legend('Equation 1 LHS', 'Equation 1 RHS', 'Equation 2 LHS',['Equation 2 ' ...
    'RHS'],'Equation 3 LHS','Equation 3 RHS')
title('GEMINI Constraining Equations #1-#3')
subtitle('E = 100 eV. Requirement: LHS > RHS')
grid on
grid minor

% Solving constraining equations but in matrix form at different energies
%EQ1 = alpha*(l1+l2+l3) + (meshV.*q.*l2./2./meshE./d).*(1+alpha^2).*(l2/2+l3);
%EQ2 = alpha*(l1+l2)+(meshV.*q.*l2^2./4./meshE./d).*(1+alpha^2);
%EQ3L = (meshV.*abs(q).*l2./2./meshE./d).*(1+alpha^2).*(l2/2+l3);
%EQ3R = ones(size(meshV))*(2*alpha*(l1+l2+l3));
EQ1 = (meshV.*q.*l2./2./meshE./d).*(l2/2+l3);
EQ2 = (meshV.*q.*l2^2./4./meshE./d);
EQ3L = (meshV.*abs(q).*l2./2./meshE./d).*(l2/2+l3);
EQ3R = ones(size(meshV))*(8);

%% Mesh surface plots to graph Figure 1 but across a range of energies
%{
figure(4);
surf(meshV, meshE, EQ1);
hold on
surf(meshV, meshE, meshr,'FaceAlpha', 0.5);
xlabel('V');
ylabel('E');

figure(5);
surf(meshV, meshE, EQ2);
hold on
surf(meshV, meshE, meshd2,'FaceAlpha', 0.5);
hold off

figure(6);
surf(meshV, meshE, EQ3L);
hold on
surf(meshV, meshE, EQ3R,'FaceAlpha', 0.5);
hold off
%}

%% Creating plot of deflector plate voltage to use vs beam energy
% Requires 2D arrays because deflector plate voltage AND beam energy change
% This is NOT the final graph to be used for GEMINI

% Creating empty arrays (row vectors) to hold values of intersection points
% Each column is a different of beam energy
XEQ1X = zeros(1,size(meshV,1));
XEQ1Y = zeros(1,size(meshV,1));
XEQ2X = zeros(1,size(meshV,1));
XEQ2Y = zeros(1,size(meshV,1));
XEQ3X = zeros(1,size(meshV,1));
XEQ3Y = zeros(1,size(meshV,1));

% Iterates through energies for the three eqns to find the intersection pts
for i = 1:size(meshV,1)
    % If there is no intersection point, then Vmax >5000V
    % Manually set the 'intersection point' (Vmax) to 5000V
    if EQ1(i,end) < meshr(i,end)
        XEQ1X(i) = 10000;
        XEQ1Y(i) = meshr(i,end);
    else
        % polyxpoly calculates the intersection point of the RHS and LHS
        [XEQ1X(i), XEQ1Y(i)] = polyxpoly([meshV(i,1) meshV(i,end)], ...
            [EQ1(i,1) EQ1(i,end)],[meshV(i,1) meshV(i,end)],[meshr(i,1) ...
            meshr(i,end)]); % Vmax
        %XEQ1X(i) = XEQ1X(i)/2; % Dividing by 2 to go from V to V/2
        %XEQ1Y(i) = XEQ1Y(i)/2;
    end
    
    % Eqn 2 LHS and RHS intersection. This is a Vmax (see Figure 1)
    if EQ2(i,end) < meshd2(i,end)
        XEQ2X(i) = 10000;
        XEQ2Y(i) = meshd2(i,end);
    else
        [XEQ2X(i), XEQ2Y(i)] = polyxpoly([meshV(i,1) meshV(i,end)], ...
            [EQ2(i,1) EQ2(i,end)],[meshV(i,1) meshV(i,end)],[meshd2(i,1) ...
            meshd2(i,end)]);
        %XEQ2X(i) = XEQ2X(i)/2;
        %XEQ2Y(i) = XEQ2Y(i)/2;
    end

    % Eqn 3 LHS and RHS intersection. This is a Vmin (see Figure 1)
    if EQ3L(i,end) < EQ3R(i,end)
        XEQ3X(i) = 10000;
        XEQ3Y(i) = EQ3R(i,end);
    else
        [XEQ3X(i), XEQ3Y(i)] = polyxpoly([meshV(i,1) meshV(i,end)], ...
            [EQ3L(i,1) EQ3L(i,end)],[meshV(i,1) meshV(i,end)],[EQ3R(i,1) ...
            EQ3R(i,end)]);
        %XEQ3X(i) = XEQ3X(i)/2;
        %XEQ3Y(i) = XEQ3Y(i)/2;
    end    
end

% Each column is an energy, each row is an intersection point (Voltage)
% XEQ1X, XEQ2X are Vmax for eqns 1 and 2 (see Figure 1)
% Take the minimum to properly constrain
Vmax = min(XEQ1X, XEQ2X);
Vrange = [Vmax; XEQ3X];
Vmean = mean(Vrange, 1);

% Dividing by 2 because each deflector plate runs at V/2
Vmax = Vmax/2;
Vrange = Vrange/2;
Vmean = Vmean/2;

% errorbar plots errorbar length so subtract mean to/from Vmax or Vmin
errneg = Vmean-Vrange(2,:); % 2nd row of Vrange is Vmin
errpos = Vrange(1,:)-Vmean; % 1st row of Vrange is Vmax

figure(2);
errorbar(EArr, Vmean, errneg, errpos, '.', 'Color','b')
ylim([0, 5000]) % y limited to 5000V because of SHV-5 connectors
grid on
grid minor
title('Deflector Plate Required Voltage (V) vs. Beam Energy (eV)')
subtitle('Bars Show Max and Min Allowable Deflector Plate Voltage')
xlabel('Beam Energy (eV)')
ylabel('Deflector Plate Required Voltage (V)')


%% Recreating Figure 2 with tuning input from SIMION

% VmaxTuned values found manually in SIMION
VmaxTuned = [125, 250, 375, 500, 625, 750, 875, 1000, 1125, 1250, 1375, ...
    1525, 1650, 1775, 1900, 2025, 2150, 2275, 2375, 2500, 2625, 2750, ...
    2875, 3000, 3150];
% Divide XEQ3X by 2 once for V/2 (line 164), once for factor 1/2 correction
VrangeTuned = [VmaxTuned; XEQ3X./4];
VmeanTuned = mean(VrangeTuned, 1);

errnegTuned = VmeanTuned-VrangeTuned(2,:); % 2nd row of Vrange is Vmin
errposTuned = VrangeTuned(1,:)-VmeanTuned; % 1st row of Vrange is Vmax

figure(3);
errorbar(EArr, VmeanTuned, errnegTuned, errposTuned, '.', 'Color','b')
ylim([0, 5000]) % y limited to 5000V because of SHV-5 connectors
grid on
grid minor
title(['Final Tuned Deflector Plate Required Voltage (V) vs. ' ...
    'Beam Energy (eV)'])
subtitle('Bars Show Max and Min Allowable Deflector Plate Voltage')
xlabel('Beam Energy (eV)')
ylabel('Deflector Plate Required Voltage (V)')