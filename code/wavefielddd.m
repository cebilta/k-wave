clear all;

% simulation settings
DATA_CAST = 'gpuArray-single';

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 20;            % [grid points]
PML_Y_SIZE = 10;            % [grid points]
PML_Z_SIZE = 10;            % [grid points]

% set total number of grid points not including the PML
Nx = 1200 - 2*PML_X_SIZE;    % [grid points]
Ny = 256 - 2*PML_Y_SIZE;    % [grid points]
Nz = 128 - 2*PML_Z_SIZE;     % [grid points]

% set desired grid size in the x-direction not including the PML
dx = 1e-3;                  % [m]
x = Nx*dx;
% calculate the spacing between the grid points                  % [m]
dy = dx;                    % [m]
dz = dx;                    % [m]

% create the k-space grid
kgrid = makeGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
medium.sound_speed = 1500;      % [m/s]
medium.density = 1000;          % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
% medium.BonA = 6;

% create the time array
t_end = x*(2.2/1500);                  % [s]
kgrid.t_array = makeTime(kgrid, medium.sound_speed, [], t_end);

% =========================================================================
% DEFINE THE Sensor
% =========================================================================

% create sensor mask
circle = makeCircle(Nx, Ny, 1, Ny/2, round(100e-2/dx), (360/180)*pi, true);
figure
imagesc(circle);
p1 = zeros(Nx,Ny,Nz);
p1(:,:,Nz/2) = circle;
p2(:,:,Nz/2) = circle;
sensor.mask = zeros(Nx,Ny,Nz);
sensor.mask(:,:,Nz/2) =circle;




% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 1000;        % [Hz]
tone_burst_cycles = 10;

% create the input signal using toneBurst 
source.p = source_strength.*toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 96;    % total number of transducer elements
transducer.element_width = 1;       % width of each element [grid points/voxels]
transducer.element_length = 10;     % length of each element [grid points/voxels]
transducer.element_spacing = 1;     % spacing (kerf  width) between the elements [grid points/voxels]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements*transducer.element_width ...
    + (transducer.number_elements - 1)*transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);

% properties used to derive the beamforming delays
transducer.sound_speed = 1500;              % sound speed [m/s]
transducer.focus_distance = 100;          % focus distance [m]
transducer.elevation_focus_distance = 10;% focus distance in the elevation plane [m]
transducer.steering_angle = 0;              % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% % define the transducer elements that are currently active
 transducer.active_elements = zeros(transducer.number_elements, 1);
 transducer.active_elements(:,:) = 1;
% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 0.7e6;        % [Hz]
tone_burst_cycles = 10;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength./(medium.sound_speed*medium.density)).*input_signal;

% append input signal used to drive the transducer
transducer.input_signal = input_signal;


% create the transducer using the defined settings
transducer = makeTransducer(kgrid, transducer);
%% transducer and the sensor
trans_index = transducer.all_elements_mask;
trans_index = uint8(p1) + trans_index;
voxelPlot(double(trans_index));

% print out transducer properties
transducer.properties;

% =========================================================================
% RUN THE SIMULATION
% =========================================================================
sensor.record = {'p_rms', 'p_max','p'};
% set the input settings
input_args = {'DisplayMask', transducer.active_elements_mask, ...
    'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    'DataCast', DATA_CAST, 'PlotScale', [-source_strength/4, source_strength/4]};

% run the simulation
sensor_data = kspaceFirstOrder3DC(kgrid, medium,transducer,sensor, input_args{:});
  
% % extract a single scan line from the sensor data using the current
% % beamforming settings
% scan_line = transducer.scan_line(sensor_data);

% % =========================================================================
% % VISUALISATION
% % =========================================================================
% 
% % plot the recorded time series
% figure;
% stackedPlot(kgrid.t_array*1e6, sensor_data);
% xlabel('Time [\mus]');
% ylabel('sensor Element');
% title('Recorded Pressure');
%% p_max plot
 figure;   
 plot(sensor_data.p_rms/1e6);
 xlabel('sensor Element(1:100)');
 ylabel('p_rms');
 title('Recorded Pressure');
% % reshape the returned rms and max fields to their original position
%     sensor_data.p_rms = reshape(sensor_data.p_rms, [Nx, Nj]);
%     sensor_data.p_max = reshape(sensor_data.p_max, [Nx, Nj]);
%     
%     % plot the beam pattern using the pressure maximum
%     figure;
%     imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, sensor_data.p_max/1e6);
%     xlabel([j_label '-position [mm]']);
%     ylabel('x-position [mm]');
%     title('Total Beam Pattern Using Maximum Of Recorded Pressure');
%     colormap(jet(256));
%     c = colorbar;
%     ylabel(c, 'Pressure [MPa]');
%     axis image;
%     
%     % plot the beam pattern using the pressure rms
%     figure;
%     imagesc(j_vec*1e3, (kgrid.x_vec - min(kgrid.x_vec(:)))*1e3, sensor_data.p_rms/1e6);
%     xlabel([j_label '-position [mm]']);
%     ylabel('x-position [mm]');
%     title('Total Beam Pattern Using RMS Of Recorded Pressure');
%     colormap(jet(256));
%     c = colorbar;
%     ylabel(c, 'Pressure [MPa]');
%     axis image;
