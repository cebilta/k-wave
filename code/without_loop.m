% Simulating B-mode Ultrasound Images Example
%
% This example illustrates how k-Wave can be used for the simulation of
% B-mode ultrasound images using a phased-array or sector transducer. It
% builds on the Simulating B-mode Ultrasound Images Example.
%
% To allow the simulated scan line data to be processed multiple times with
% different settings, the simulated RF data is saved to disk. This can be
% reloaded by setting RUN_SIMULATION = false within the example m-file. The
% data can also be downloaded from
% http://www.k-wave.org/datasets/example_us_phased_array_scan_lines.mat 
%
% author: Bradley Treeby
% date: 7th September 2012
% last update:  26th August 2014
%  
% This function is part of the k-Wave Toolbox (http://www.k-wave.org)
% Copyright (C) 2009-2014 Bradley Treeby and Ben Cox

% This file is part of k-Wave. k-Wave is free software: you can
% redistribute it and/or modify it under the terms of the GNU Lesser
% General Public License as published by the Free Software Foundation,
% either version 3 of the License, or (at your option) any later version.
% 
% k-Wave is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for
% more details. 
% 
% You should have received a copy of the GNU Lesser General Public License
% along with k-Wave. If not, see <http://www.gnu.org/licenses/>. 

%#ok<*UNRCH>

clear all;
% addpath('D:\Hareesh\YSC\K-Wave-toolbox-1.1\k-Wave')
% simulation settings
DATA_CAST       = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = true;  % set to false to reload previous results instead of running simulation

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
pml_x_size = 15;            % [grid points]
pml_y_size = 10;            % [grid points]
pml_z_size = 10;            % [grid points]

% set total number of grid points not including the PML
sc = 1;
Nx = 256/sc - 2*pml_x_size;   % [grid points]
Ny = 256/sc - 2*pml_y_size;   % [grid points]
Nz = 128/sc - 2*pml_z_size;   % [grid points]

% set desired grid size in the x-direction not including the PML
x = 50e-3;                  % [m]

% calculate the spacing between the grid points
dx = x/Nx;                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
c0 = 1540;                  % [m/s]
rho0 = 1000;                % [kg/m^3]
medium.alpha_coeff = 0.75; 	% [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% create the time array
t_end = (Nx*dx)*2.2/c0;     % [s]
kgrid.t_array = makeTime(kgrid, c0, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;    	% [Pa]
tone_burst_freq = 1e6/sc;  	% [Hz]
tone_burst_cycles = 4;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength./(c0*rho0)).*input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% define the physical properties of the phased array transducer
transducer.number_elements = 64/sc;         % total number of transducer elements
transducer.element_width = 1;               % width of each element [grid points]
transducer.element_length = 40/sc;          % length of each element [grid points]
transducer.element_spacing = 0;             % spacing (kerf  width) between the elements [grid points]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements*transducer.element_width ...
    + (transducer.number_elements - 1)*transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = c0;                % sound speed [m/s]
transducer.focus_distance = 30e-3;          % focus distance [m]
transducer.elevation_focus_distance = 30e-3;% focus distance in the elevation plane [m]
transducer.steering_angle = 0;              % steering angle [degrees]
transducer.steering_angle_max = 32;         % maximum steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Hanning';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
transducer.active_elements = ones(transducer.number_elements, 1);

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = makeTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% =========================================================================
% DEFINE THE MEDIUM PROPERTIES
% =========================================================================

% define a random distribution of scatterers for the medium
background_map_mean = 1;
background_map_std = 0.008;
background_map = background_map_mean + background_map_std*randn([Nx, Ny, Nz]);

% define a random distribution of scatterers for the highly scattering
% region
scattering_map = randn([Nx, Ny, Nz]);
scattering_c0 = c0 + 150 + 75*scattering_map;
scattering_c0(scattering_c0 > 1800) = 1800;
scattering_c0(scattering_c0 < 1500) = 1500;
scattering_rho0 = scattering_c0/1.5;

% define properties
sound_speed_map = c0*ones(Nx, Ny, Nz).*background_map;
density_map = rho0*ones(Nx, Ny, Nz).*background_map;
% scattering_c0 = 1600* ones(Nx, Ny, Nz);
% scattering_rho0 = scattering_c0 / 1.5;
% 
% % define properties
% sound_speed_map = c0 * ones(Nx, Ny, Nz);
% density_map = rho0 * ones(Nx, Ny, Nz);
% define a sphere for a highly scattering region
% radius = 8e-3;
% x_pos = 32e-3;
% y_pos = dy*Ny/2;
% scattering_region1 = makeBall(Nx, Ny, Nz, round(x_pos/dx), round(y_pos/dx), Nz/2, round(radius/dx));
% 
% % assign region
% sound_speed_map(scattering_region1 == 1) = scattering_c0(scattering_region1 == 1);
% density_map(scattering_region1 == 1) = scattering_rho0(scattering_region1 == 1);
%assign region 0
cylinder = zeros(226,236);
p = zeros(226,236,108);
for r =45:50
    cylinder_1 = repmat(makeCircle(226,236,100,150,r),[1,1,108]);
    cylinder =cylinder+cylinder_1;
end
cylinder(:,:,1:50) = p(:,:,1:50);
cylinder(:,:,58:108) = p(:,:,58:108);
voxelPlot(cylinder);
% assign region cylinder
sound_speed_map(cylinder == 1) = scattering_c0(cylinder == 1);
density_map(cylinder == 1) = scattering_rho0(cylinder == 1);
% % assign region 1
% sound_speed_map(80,140:180,50:60) = scattering_c0(80,140:180,50:60);
% density_map(80,140:180,50:60) = scattering_rho0(80,140:180,50:60);
% % assign region 2
% sound_speed_map(120,140:180,50:60) = scattering_c0(120,140:180,50:60);
% density_map(120,140:180,50:60) = scattering_rho0(120,140:180,50:60);
% % assign region 3
% sound_speed_map(90,140:180,50:60) = scattering_c0(90,140:180,50:60);
% density_map(90,140:180,50:60) = scattering_rho0(90,140:180,50:60);
% % assign region 4
% sound_speed_map(100,140:180,50:60) = scattering_c0(100,140:180,50:60);
% density_map(100,140:180,50:60) = scattering_rho0(100,140:180,50:60);
% % assign region 5
% sound_speed_map(110,140:180,50:60) = scattering_c0(110,140:180,50:60);
% density_map(110,140:180,50:60) = scattering_rho0(110,140:180,50:60);
% % assign region 6
% sound_speed_map(80:120,138,50:60) = scattering_c0(80:120,138,50:60);
% density_map(80:120,138,50:60) = scattering_rho0(80:120,138,50:60);
% % assign region 7
% sound_speed_map(80:120,148,50:60) = scattering_c0(80:120,148,50:60);
% density_map(80:120,148,50:60) = scattering_rho0(80:120,148,50:60);
% % assign region 8
% sound_speed_map(80:120,158,50:60) = scattering_c0(80:120,158,50:60);
% density_map(80:120,158,50:60) = scattering_rho0(80:120,158,50:60);
% % assign region 9
% sound_speed_map(80:120,168,50:60) = scattering_c0(80:120,168,50:60);
% density_map(80:120,168,50:60) = scattering_rho0(80:120,168,50:60);
% % assign region 9
% sound_speed_map(80:120,178,50:60) = scattering_c0(80:120,178,50:60);
% density_map(80:120,178,50:60) = scattering_rho0(80:120,178,50:60);
% voxelPlot(sound_speed_map);
% p = ones(216,300,108);
% line_1(175,40:50,48:52) = p(75,40:50,48:52);
% line = zeros(216,300,108);
% line(:,180,:) = makeLine(216, 108, [50,100], 0.6, 150);
% line(:,179,:) = line(:,180,:);
% line(:,181,:) = line(:,180,:);
% line(:,200,:) = makeLine(216, 108, [50,100], 0.6, 150);
% line(:,220,:) = makeLine(216, 108, [50,100], 0.6, 150);
% line(:,240,:) = makeLine(216, 108, [50,100], 0.6, 150);
%  line_1 = zeros(226,236,108);
% line_2 = repmat(makeLine(226, 236, [50,100], 0.6, 15),[1,1,108]);
% line_3 = repmat(makeLine(226, 236, [50,130], 0.6, 15),[1,1,108]);
% line_4 = repmat(makeLine(226, 236, [50,160], 0.6, 15),[1,1,108]);
% line_5 = repmat(makeLine(226, 236, [50,190], 0.6, 15),[1,1,108]);
% line_1 = line_2 + line_3 + line_4 + line_5;
% voxelPlot(line_1);
% p = ones(216,300,108);
% line_1(175,40:50,48:52) = p(75,40:50,48:52);
% line = zeros(216,300,108);
% line(:,180,:) = makeLine(216, 108, [50,100], 0.6, 150);
% line(:,179,:) = line(:,180,:);
% line(:,181,:) = line(:,180,:);
% line(:,200,:) = makeLine(216, 108, [50,100], 0.6, 150);
% line(:,220,:) = makeLine(216, 108, [50,100], 0.6, 150);
% line(:,240,:) = makeLine(216, 108, [50,100], 0.6, 150);
% voxelPlot(line_1);
% % assign region
% sound_speed_map(line == 1) = scattering_c0(line == 1);
% density_map(line == 1) = scattering_rho0(line == 1);
% % assign region
% sound_speed_map(line_1 == 1) = scattering_c0(line_1 == 1);
% density_map(line_1 == 1) = scattering_rho0(line_1 == 1);
% 
% 
% 
% assign to the medium inputs
medium.sound_speed = sound_speed_map;
medium.density = density_map;


% =========================================================================
% RUN THE SIMULATION
% =========================================================================

% range of steering angles to test
steering_angles = -32:0.5:32;

% preallocate the storage
number_scan_lines = length(steering_angles);
scan_lines = zeros(number_scan_lines, kgrid.Nt);

% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

% run the simulation if set to true, otherwise, load previous results
if RUN_SIMULATION
     
    % loop through the range of angles to test
    for angle_index = 1:number_scan_lines
        
        % update the command line status
        disp('');
        disp(['Computing scan line ' num2str(angle_index) ' of ' num2str(number_scan_lines)]);

        % update the current steering angle
        transducer.steering_angle = steering_angles(angle_index);
        
       % run the simulation
        sensor_data = kspaceFirstOrder3DC(kgrid, medium, transducer, transducer, input_args{:});
    
         %sensor_data=randn(64,1658);
        % extract the scan line from the sensor data
        scan_lines(angle_index, :) = transducer.scan_line(sensor_data);
        
    end

    % save the scan lines to disk
    save example_us_phased_array_scan_lines scan_lines;
    
else
    
    % load the scan lines from disk
    %load example_us_phased_array_scan_lines
    load scan_lines_64.mat
end

% trim the delay offset from the scan line data
t0_offset = round(length(input_signal)/2) + (transducer.appended_zeros - transducer.beamforming_delays_offset);
scan_lines = scan_lines(:, t0_offset:end);

% get the new length of the scan lines
Nt = length(scan_lines(1, :));

% =========================================================================
% PROCESS THE RESULTS
% =========================================================================

% -----------------------------
% Remove Input Signal
% -----------------------------

% create a window to set the first part of each scan line to zero to remove
% interference from the input signal
scan_line_win = getWin(Nt*2, 'Tukey', 'Param', 0.05).';
scan_line_win = [zeros(1, t0_offset*2), scan_line_win(1:end/2 - t0_offset*2)];

% apply the window to each of the scan lines
scan_lines = bsxfun(@times, scan_line_win, scan_lines);

% -----------------------------
% Time Gain Compensation
% -----------------------------

% create radius variable
r = c0*(1:Nt)*kgrid.dt/2;    % [m]

% create time gain compensation function based on attenuation value,
% transmit frequency, and round trip distance
tgc_alpha = 0.4;       % [dB/(MHz cm)]
tgc = exp(2*tgc_alpha*tone_burst_freq/1e6*r*100);

% apply the time gain compensation to each of the scan lines
scan_lines = bsxfun(@times, tgc, scan_lines);

% -----------------------------
% Frequency Filtering
% -----------------------------

% filter the scan lines using both the transmit frequency and the second
% harmonic
scan_lines_fund = gaussianFilter(scan_lines, 1/kgrid.dt, tone_burst_freq, 100, true);
scan_lines_harm = gaussianFilter(scan_lines, 1/kgrid.dt, 2*tone_burst_freq, 30, true);

% -----------------------------
% Envelope Detection
% -----------------------------

% envelope detection
scan_lines_fund = envelopeDetection(scan_lines_fund);
scan_lines_harm = envelopeDetection(scan_lines_harm);

% -----------------------------
% Log Compression
% -----------------------------

% normalised log compression
compression_ratio = 3;
scan_lines_fund = logCompression(scan_lines_fund, compression_ratio, true);
scan_lines_harm = logCompression(scan_lines_harm, compression_ratio, true);

% -----------------------------
% Scan Conversion
% -----------------------------

% set the desired size of the image
image_size = [Nx*dx, Ny*dy];

% convert the data from polar coordinates to cartesian coordinates for
% display
b_mode_fund = scanConversion(scan_lines_fund, steering_angles, image_size, c0, kgrid.dt);
b_mode_harm = scanConversion(scan_lines_harm, steering_angles, image_size, c0, kgrid.dt);

% =========================================================================
% VISUALISATION
% =========================================================================

% create the axis variables
x_axis = [0, Nx*dx*1e3];    % [mm]
y_axis = [0, Ny*dy*1e3];    % [mm]

% plot the raw scan lines
figure;
imagesc(steering_angles, x_axis, scan_lines.');
colormap(gray);
xlabel('Steering angle [deg]');
ylabel('Depth [mm]');
title('Raw Scan-Line Data');

% plot the medium and the B-modes images
figure;

imagesc(y_axis, x_axis, medium.sound_speed(:, :, end/2));
axis image;
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');
title('Scattering Phantom');

figure;
imagesc(y_axis, x_axis, b_mode_fund);
axis image;
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');
title('B-Mode Image');

figure;
imagesc(y_axis, x_axis, b_mode_harm);
colormap(gray);
axis image;
xlabel('Horizontal Position [mm]');
ylabel('Depth [mm]');
title('Harmonic Image');

scaleFig(2, 1);