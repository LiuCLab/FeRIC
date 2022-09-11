%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FeRIC Coil Magnetic and Electric Field Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Victor Han
% Last Modified: 9/11/22
% Based on: http://openems.de/index.php/Tutorial:_MRI_Loop_Coil
%
% This code simulates the magnetic and electric fields of a two-layer
% solenoid coil with a saline solution dish inside. All units are SI
% (Tesla, V/m, meter, etc).
%
% Please first install openEMS
% (https://openems.de/index.php/OpenEMS#Installation). Then, simply run
% this script. Close the coil geometry viewer when it appears to start the
% simulation.
%
% Inputs: None
%
% Outputs: 2D plots of xy and xz cross sections for the electric and
% magnetic fields of a 1-turn coil at 180 MHz. A temporary simulation
% folder is also created.
%
% Tested with
%  - openEMS v0.0.35
%  - Matlab R2019a
%  - Windows 10
%
% Note: Takes several hours to run on a laptop

close all
clear
clc

large_coil = false; % Set this to true if simulating a large coil

if large_coil
    coil_diameter = 0.05;
else
    coil_diameter = 0.028;
end

%% General Setup
physical_constants; % Sets some physical constants in SI units
unit = 1; % Sets length scale to meters
B_norm = 1.6e-6; % Sets magnetic field strength in center of simulation in Tesla

% Set the box size for saving and visualizing the fields
visualize_box.start = [-0.1 -0.1 -0.1];
visualize_box.stop  = [0.1 0.1 0.1];

% Set the box size for a high resolution mesh. The simulation has a
% variable spatial resolution, so this box is where we decide to use more
% computational power and a higher resolution for more accuracy and finer
% details.

mesh_box.start = [-0.050 -0.050 -0.010]; 
mesh_box.stop  = [+0.050 +0.050 +0.010];
mesh_box.resolution = 0.0004;

Air_Box = 0.250;      % Size of the surrounding air box

%% Setup Simulation Parameters
FDTD = InitFDTD( 'NrTS', 6000000, 'EndCriteria', 1e-8, 'CellConstantMaterial', 0);

% Define excitation signal pulse. Change these values to change simulated
% frequency. f0 is the frequency at which fields are calculated later, and
% fc determines the shape of the time-domain pulse. Note that because this
% is a time-domain simulation method, lower frequencies lead to longer
% time-domain pulses and thus often much longer simulation times for the
% same geometry.
f0 = 180e6; % Center frequency
fc = 200e6; % 20 dB corner frequency
FDTD = SetGaussExcite( FDTD, f0, fc );

% Setup boundary conditions
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'};
FDTD = SetBoundaryCond( FDTD, BC );

%% Setup CSXCAD geometry & mesh
CSX = InitCSX();

% We are modeling our coil as a helix.
Helix.radius = coil_diameter/2; % Change this number for the size of the coil
Helix.turns = 1; 
Helix.pitch = 0; 

% Excitation signal port
feed.height = -0.001; % This is the z-position where the coil starts
feed.R = 500;    % Port impedance. Arbitrary for our purposes because we normalize later
portlen = 0.010; % Length of wire from port to coil
start = [-0.001 Helix.radius+portlen feed.height+Helix.pitch*Helix.turns];
stop  = [0.001 Helix.radius+portlen feed.height];
[CSX, port] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [1 0 0], true);

CSX = AddMetal( CSX, 'helix' ); % Create helix as a perfect electric conductor (PEC)

% Generate the helix points
ang = linspace(pi/16,2*pi-pi/16,21) + pi/2;
coil_x = Helix.radius*cos(ang);
coil_y = Helix.radius*sin(ang);
coil_z = 0*ang + feed.height;

helix.x=[];
helix.y=[];
helix.z=[];
zpos = feed.height;

% Start with a point at the port
helix.y = [helix.y Helix.radius+portlen];
helix.x = [helix.x -0.001];
helix.z = [helix.z zpos];

% Add all of the helix points
for n=0:Helix.turns-1
    helix.x = [helix.x coil_x];
    helix.y = [helix.y coil_y];
    helix.z = [helix.z coil_z];
    zpos = zpos + Helix.pitch;
end

% Add point at the other side of the port
helix.y = [helix.y Helix.radius+portlen];
helix.x = [helix.x 0.001];
helix.z = [helix.z zpos];

clear p
p(1,:) = helix.x;
p(2,:) = helix.y;
p(3,:) = helix.z;
CSX = AddWire(CSX, 'helix', 0, p, 0.0005);


%% Add the saline solution
CSX = AddMaterial(CSX, 'saline_solution');
CSX = SetMaterialProperty( CSX,'saline_solution', 'Epsilon', 80, 'Kappa', 1.5, 'Density', 1000);

clear p
major_d = 0.02;
minor_d = 0.015;
p(1,1) = -major_d/2; p(2,1) = 0;
p(1,2) = -major_d/4; p(2,2) = minor_d/3.3;
p(1,3) = -major_d/6; p(2,3) = minor_d/2.5;
p(1,4) = -major_d/20; p(2,4) = minor_d/2;
p(1,5) = major_d/20; p(2,5) = minor_d/2;
p(1,6) = major_d/6; p(2,6) = minor_d/2.5;
p(1,7) = major_d/4; p(2,7) = minor_d/3.3;
p(1,8) = major_d/2; p(2,8) = 0;
p(1,9) = major_d/4; p(2,9) = -minor_d/3.3;
p(1,10) = major_d/6; p(2,10) = -minor_d/2.5;
p(1,11) = major_d/20; p(2,11) = -minor_d/2;
p(1,12) = -major_d/20; p(2,12) = -minor_d/2;
p(1,13) = -major_d/6; p(2,13) = -minor_d/2.5;
p(1,14) = -major_d/4; p(2,14) = -minor_d/3.3;
p(1,15) = -major_d/2; p(2,15) = 0;
CSX = AddLinPoly(CSX, 'saline_solution', 1, 2, 0, p, 0.003);

clear p
channel_l = 0.016;
channel_w = 0.002;
p(1,1) = major_d/2 - channel_w; p(2,1) = minor_d/3.3 * channel_w/(major_d/4);
p(1,2) = major_d/2; p(2,2) = 0;

angles = linspace(1.2*pi, 0.9*pi, 3);
for ii=1:2
    p(1,ii+2) = channel_l + 0.003*cos(angles(ii));
    p(2,ii+2) = 0.003 + 0.003*sin(angles(ii));
end

CSX = AddLinPoly(CSX, 'saline_solution', 2, 2, 0, p, 0.002);

CSX = AddCylinder(CSX, 'saline_solution', 1, [channel_l 0.003 0], [channel_l 0.003 0.003], 0.003);


%% Add headstage input resistance
port_R = 500e6;    % Port impedance. 
start = [0 -0.05 0.05];
stop  = [0 -0.05 0.049];
[CSX, port_receive] = AddLumpedPort(CSX, 5, 2, port_R, start, stop, [0 0 1], false);


%% Add electrode
CSX = AddMetal( CSX, 'electrode' ); % Create helix as a perfect electric conductor (PEC)

CSX = AddCylinder(CSX, 'electrode', 0, [0 -0.001 0.001], [0 -0.05 0.05], 0.0001);

clear p
p(1,1) = 0; p(2,1) = -0.05; p(3,1) = 0.049;
p(1,2) = 0; p(2,2) = -0.05; p(3,2) = 0.005;
p(1,3) = channel_l; p(2,3) = -0.05; p(3,3) = 0.005;
p(1,4) = channel_l; p(2,4) = 0.003; p(3,4) = 0.005;
p(1,5) = channel_l; p(2,5) = 0.003; p(3,5) = 0.001;

CSX = AddCylinder(CSX, 'electrode', 0, p(:,1).', p(:,2).', 0.0001);
CSX = AddCylinder(CSX, 'electrode', 0, p(:,2).', p(:,3).', 0.0001);
CSX = AddCylinder(CSX, 'electrode', 0, p(:,3).', p(:,4).', 0.0001);
CSX = AddCylinder(CSX, 'electrode', 0, p(:,4).', p(:,5).', 0.0001);

%% Finalize mesh
mesh = DetectEdges(CSX);

% Add a dense mesh defined earlier for the coil
mesh.x = [mesh.x mesh_box.start(1) mesh_box.stop(1)];
mesh.y = [mesh.y mesh_box.start(2) mesh_box.stop(2)];
mesh.z = [mesh.z mesh_box.start(3) mesh_box.stop(3)];

% Smooth the mesh
mesh = SmoothMesh(mesh, mesh_box.resolution);

% Add air spacer
mesh.x = [-Air_Box+mesh.x(1) mesh.x mesh.x(end)+Air_Box];
mesh.y = [-Air_Box+mesh.y(1) mesh.y mesh.y(end)+Air_Box];
mesh.z = [-Air_Box+mesh.z(1) mesh.z mesh.z(end)+Air_Box];

mesh = SmoothMesh(mesh, c0 / (f0+fc) / unit / 30, 1.5, 'algorithm', 1);

%% Add Dump boxes (2D boxes) for H and E on xy- and xz-plane
offset = [0 0 0]; % Offset to the visualization plane if the saline solution is off-center
CSX = AddDump(CSX,'Hf_xy_0','DumpType',11,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Hf_xy_0',0, visualize_box.start.*[1 1 0] + offset, visualize_box.stop.*[1 1 0] + offset);
CSX = AddDump(CSX,'Hf_xz_0','DumpType',11,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Hf_xz_0',0, visualize_box.start.*[1 0 1], visualize_box.stop.*[1 0 1]);

CSX = AddDump(CSX,'Ef_xy_0','DumpType',10,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Ef_xy_0',0, visualize_box.start.*[1 1 0] + offset, visualize_box.stop.*[1 1 0] + offset);
CSX = AddDump(CSX,'Ef_xz_0','DumpType',10,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Ef_xz_0',0, visualize_box.start.*[1 0 1], visualize_box.stop.*[1 0 1]);


%% Add 10 lines in all direction to make space for PML or MUR absorbing
mesh = AddPML(mesh, 10);

%% Finaly define the FDTD mesh grid
disp(['number of cells: ' num2str(1e-6*numel(mesh.x)*numel(mesh.y)*numel(mesh.z)) ' Mcells'])
CSX = DefineRectGrid( CSX, unit, mesh );

%% Prepare simulation folder
Sim_Path = ['tmp_' mfilename];
Sim_CSX = [mfilename '.xml'];

[~, ~, ~] = rmdir( Sim_Path, 's' ); % Clear previous directory
[status, message, messageid] = mkdir( Sim_Path ); % Create empty simulation folder

%% Write openEMS compatible xml-file
WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

%% Show the structure and export as vtk data automatically
CSXGeomPlot( [Sim_Path '/' Sim_CSX] , ['--export-polydata-vtk=' Sim_Path ' --RenderDiscMaterial -v']);

%% Run openEMS
RunOpenEMS( Sim_Path, Sim_CSX);

%% Calculate scaling factor
% Now the simulation is done and it is time to postprocess.
% Read the H field and calculate scaling factor for simulation to make the
% center B field equal to the value set earlier
physical_constants; % Sets some physical constants in SI units
unit = 1; % Sets length scale to meters
[H_field, H_mesh] = ReadHDF5Dump([Sim_Path '/Hf_xz_0.h5']);

Bx = MUE0*H_field.FD.values{1}(:,:,:,1);
By = MUE0*H_field.FD.values{1}(:,:,:,2);
Bz = MUE0*H_field.FD.values{1}(:,:,:,3);
Btot = sqrt(abs(Bx).^2 + abs(By).^2 + abs(Bz).^2); % B field magnitude
scale = Btot(find(H_mesh.lines{1}==0,1), find(H_mesh.lines{3}==0,1)) / B_norm; % Scaling factor for E and B fields.

%% Plot E field magnitude
% First plot in the xz plane
[E_field, E_mesh] = ReadHDF5Dump([Sim_Path '/Ef_xz_0.h5']);
Ex = E_field.FD.values{1}(:,:,:,1)/scale;
Ey = E_field.FD.values{1}(:,:,:,2)/scale;
Ez = E_field.FD.values{1}(:,:,:,3)/scale;

% Create a 2D grid to plot on
[X, Z] = ndgrid(E_mesh.lines{1},E_mesh.lines{3});
% Get E field magnitude
E = sqrt(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2);
% Plot
figure()
subplot(1,2,1);
h = pcolor(X,Z,squeeze(abs(E)));
set(h,'EdgeColor','none');
xlabel('x (m)');
ylabel('z (m)');
title('E field in xz (V/m), y = 0');
axis equal tight
colorbar

% Now plot in the xy plane
[E_field, E_mesh] = ReadHDF5Dump([Sim_Path '/Ef_xy_0.h5']);
Ex = E_field.FD.values{1}(:,:,:,1)/scale;
Ey = E_field.FD.values{1}(:,:,:,2)/scale;
Ez = E_field.FD.values{1}(:,:,:,3)/scale;

% create a 2D grid to plot on
[X, Y] = ndgrid(E_mesh.lines{1},E_mesh.lines{2});
% Get E field magnitude
E = sqrt(abs(Ex).^2 + abs(Ey).^2 + abs(Ez).^2);
% Plot
subplot(1,2,2);
h = pcolor(X,Y,abs(E));
set(h,'EdgeColor','none');
xlabel('x (m)');
ylabel('y (m)');
title('E field in xy (V/m), z = 0 mm');
axis equal tight
colorbar

%% Plot B field magnitude
% First plot in the xz plane
[B_field, B_mesh] = ReadHDF5Dump([Sim_Path '/Hf_xz_0.h5']);
Bx = MUE0*B_field.FD.values{1}(:,:,:,1)/scale;
By = MUE0*B_field.FD.values{1}(:,:,:,2)/scale;
Bz = MUE0*B_field.FD.values{1}(:,:,:,3)/scale;

% Create a 2D grid to plot on
[X, Z] = ndgrid(B_mesh.lines{1},B_mesh.lines{3});
% Get B field magnitude
B = sqrt(abs(Bx).^2 + abs(By).^2 + abs(Bz).^2);
% Plot
figure()
subplot(1,2,1);
h = pcolor(X,Z,squeeze(abs(B)));
set(h,'EdgeColor','none');
xlabel('x (m)');
ylabel('z (m)');
title('B field in xz (T), y = 0');
axis equal tight
colorbar

% Now plot in the xy plane
[B_field, B_mesh] = ReadHDF5Dump([Sim_Path '/Hf_xy_0.h5']);
Bx = MUE0*B_field.FD.values{1}(:,:,:,1)/scale;
By = MUE0*B_field.FD.values{1}(:,:,:,2)/scale;
Bz = MUE0*B_field.FD.values{1}(:,:,:,3)/scale;

% Create a 2D grid to plot on
[X, Y] = ndgrid(B_mesh.lines{1},B_mesh.lines{2});
% Get B field magnitude
B = sqrt(abs(Bx).^2 + abs(By).^2 + abs(Bz).^2);
% Plot
subplot(1,2,2);
h = pcolor(X,Y,abs(B));
set(h,'EdgeColor','none');
xlabel('x (m)');
ylabel('y (m)');
title('B field in xy (T), z = 0 mm');
axis equal tight
colorbar

%% Print out some values of interest
disp('E field at center:');
disp(E(find(E_mesh.lines{1}==0,1),find(E_mesh.lines{2}==0,1)));

freq = 180e6;
port_receive = calcPort(port_receive, Sim_Path, freq);
disp('Voltage across headstage impedance:');
disp(abs(port_receive.uf.tot/scale));