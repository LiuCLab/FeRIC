%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FeRIC Coil Magnetic and Electric Field Simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: Victor Han
% Last Modified: 6/25/20
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
% magnetic fields of a 2-turn coil at 180 MHz. A temporary simulation
% folder is also created.
%
% Tested with
%  - openEMS v0.0.35
%  - Matlab R2019a

close all
clear
clc

%% General Setup
physical_constants; % Sets some physical constants in SI units
unit = 1; % Sets length scale to meters
B_norm = 12e-6; % Sets magnetic field strength in center of simulation in Tesla

% Set the box size for saving and visualizing the fields
visualize_box.start = [-0.1 -0.1 -0.1];
visualize_box.stop  = [0.1 0.1 0.1];

% Set the box size for a high resolution mesh. The simulation has a
% variable spatial resolution, so this box is where we decide to use more
% computational power and a higher resolution for more accuracy and finer
% details.

mesh_box.start = [-0.030 -0.030 -0.010];
mesh_box.stop  = [+0.030 +0.030 +0.010];
mesh_box.resolution = 0.0005;

Air_Box = 0.250;      % Size of the surrounding air box (150mm)

%% Setup Simulation Parameters
FDTD = InitFDTD( 'EndCriteria', 5e-7, 'CellConstantMaterial', 0);

% Define excitation signal pulse
f0 = 180e6; % Center frequency
fc = 170e6; % 20 dB corner frequency
FDTD = SetGaussExcite( FDTD, f0, fc );

% Setup boundary conditions
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'MUR'};
FDTD = SetBoundaryCond( FDTD, BC );

%% Setup CSXCAD geometry & mesh
CSX = InitCSX();

% We are modeling our coil as a helix.
Helix.radius = 0.025;
Helix.turns = 2; 
Helix.pitch = 0.002; 

% Excitation signal port
feed.height = -Helix.pitch*Helix.turns/2; % This is the z-position where the coil starts
feed.R = 500;    % Port impedance. Arbitrary for our purposes because we normalize later
portlen = 0.010; % Length of wire from port to coil
start = [Helix.radius+portlen 0 feed.height+Helix.pitch*Helix.turns];
stop  = [Helix.radius+portlen 0 feed.height];
[CSX, port] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 0 1], true);

CSX = AddMetal( CSX, 'helix' ); % Create helix as a perfect electric conductor (PEC)

% Generate the helix points
ang = linspace(0,2*pi,21);
coil_x = Helix.radius*cos(ang);
coil_y = Helix.radius*sin(ang);
coil_z = ang/2/pi*Helix.pitch;

helix.x=[];
helix.y=[];
helix.z=[];
zpos = feed.height;

% Start with a point at the port
helix.x = [helix.x Helix.radius+portlen];
helix.y = [helix.y 0];
helix.z = [helix.z zpos];

% Add all of the helix points
for n=0:Helix.turns-1
    helix.x = [helix.x coil_x];
    helix.y = [helix.y coil_y];
    helix.z = [helix.z coil_z+zpos];
    zpos = zpos + Helix.pitch;
end

% Add point at the other side of the port
helix.x = [helix.x Helix.radius+portlen];
helix.y = [helix.y 0];
helix.z = [helix.z zpos];

clear p
p(1,:) = helix.x;
p(2,:) = helix.y;
p(3,:) = helix.z;
CSX = AddCurve(CSX, 'helix', 0, p);

%% Add the saline solution
CSX = AddMaterial(CSX, 'saline_solution');
CSX = SetMaterialProperty( CSX,'saline_solution', 'Epsilon', 80, 'Kappa', 1.5, 'Density', 1000);
CSX = AddCylinder(CSX, 'saline_solution', 0, [0 0 -0.005], [0 0 0], 0.0175);

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
offset = [0 0 -0.0025]; % Offset to the visualization plane because the saline solution is off-center
CSX = AddDump(CSX,'Hf_xy','DumpType',11,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Hf_xy',0, visualize_box.start.*[1 1 0] + offset, visualize_box.stop.*[1 1 0] + offset);
CSX = AddDump(CSX,'Hf_xz','DumpType',11,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Hf_xz',0, visualize_box.start.*[1 0 1], visualize_box.stop.*[1 0 1]);

CSX = AddDump(CSX,'Ef_xy','DumpType',10,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Ef_xy',0, visualize_box.start.*[1 1 0] + offset, visualize_box.stop.*[1 1 0] + offset);
CSX = AddDump(CSX,'Ef_xz','DumpType',10,'FileType',1,'Frequency',f0);
CSX = AddBox(CSX,'Ef_xz',0, visualize_box.start.*[1 0 1], visualize_box.stop.*[1 0 1]);

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
[H_field, H_mesh] = ReadHDF5Dump([Sim_Path '/Hf_xy.h5']);

Bx = MUE0*H_field.FD.values{1}(:,:,:,1);
By = MUE0*H_field.FD.values{1}(:,:,:,2);
Bz = MUE0*H_field.FD.values{1}(:,:,:,3);
Btot = sqrt(abs(Bx).^2 + abs(By).^2 + abs(Bz).^2); % B field magnitude
ind = ceil(size(Btot)/2); % Index for center of simulation
scale = Btot(ind(1), ind(2)) / B_norm; % Scaling factor for E and B fields.

%% Plot E field magnitude
% First plot in the xz plane
[E_field, E_mesh] = ReadHDF5Dump([Sim_Path '/Ef_xz.h5']);
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
title('E field in xz (V/m), y = 0 m');
axis equal tight
colorbar

% Now plot in the xy plane
[E_field, E_mesh] = ReadHDF5Dump([Sim_Path '/Ef_xy.h5']);
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
title('E field in xy (V/m), z = -0.0025 m');
axis equal tight
colorbar

%% Plot B field magnitude
% First plot in the xz plane
[B_field, B_mesh] = ReadHDF5Dump([Sim_Path '/Hf_xz.h5']);
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
title('B field in xz (T), y = 0 m');
axis equal tight
colorbar

% Now plot in the xy plane
[B_field, B_mesh] = ReadHDF5Dump([Sim_Path '/Hf_xy.h5']);
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
title('B field in xy (T), z = -0.0025 m');
axis equal tight
colorbar
