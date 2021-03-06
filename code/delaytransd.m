clear all;

% simulation settings
DATA_CAST = 'gpuArray-single';

% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
PML_X_SIZE = 30;            % [grid points]
PML_Y_SIZE = 15;            % [grid points]
PML_Z_SIZE = 10;            % [grid points]

% set total number of grid points not including the PML
Nx = 512 - 2*PML_X_SIZE;    % [grid points]
Ny = 256 - 2*PML_Y_SIZE;    % [grid points]
Nz = 128 - 2*PML_Z_SIZE;     % [grid points]

% set desired grid size in the x-direction not including the PML
                 % [m]
dx = 0.5e-3;                  % [m]
x = Nx*dx;
% calculate the spacing between the grid points
                % [m]
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
medium.BonA = 6;

% create the time array
t_end = x*(2.2/1500);                  % [s]
kgrid.t_array = makeTime(kgrid, medium.sound_speed, [], t_end);

% =========================================================================
% DEFINE THE SOURCE
% =========================================================================

% create source mask
source.p_mask = makeBall(Nx, Ny, Nz,round(2e-1/dx), Ny/2, Nz/2, 3);%...
%     + makeBall(Nx, Ny, Nz, round(8e-3/dx), Ny/4, Nz/2, 3);
p1 = makeBall(Nx, Ny, Nz,round(2e-1/dx), Ny/2, Nz/2, 3);
% =========================================================================
% DEFINE THE SOURCE
% =========================================================================

% % create source mask
% circle = makeCircle(Nx, Ny, 1, Ny/2, 300, (360/180)*pi, true);
% % figure
% % imagesc(circle);
% source.p_mask = zeros(Nx,Ny,Nz);
% source.p_mask(:,:,Nz/2) =circle;
% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 0.7e6;        % [Hz]
tone_burst_cycles = 5;

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
transducer.sound_speed = 1540;              % sound speed [m/s]
transducer.focus_distance = 2.5e-1;          % focus distance [m]
transducer.elevation_focus_distance = 19e-3;% focus distance in the elevation plane [m]
transducer.steering_angle = 0;              % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Rectangular';    
transducer.receive_apodization = 'Rectangular';

% % define the transducer elements that are currently active
% transducer.active_elements = zeros(transducer.number_elements, 1);
% transducer.active_elements(21:52) = 1;

% create the transducer using the defined settings
transducer = makeTransducer(kgrid, transducer);
%% transducer and the sensor
trans_index = transducer.all_elements_mask;
trans_index = uint8(p1) + trans_index;
voxelPlot(double(trans_index));


% print out transducer properties
transducer.properties;
%plot the transducer
transducer.plot;

% =========================================================================
% RUN THE SIMULATION
% =========================================================================
sensor.record = {'p_rms', 'p_max','p'};
% set the input settings
input_args = {'DisplayMask', transducer.active_elements_mask, ...
    'PMLInside', false, 'PlotPML', false, 'PMLSize', [PML_X_SIZE, PML_Y_SIZE, PML_Z_SIZE], ...
    'DataCast', DATA_CAST, 'PlotScale', [-source_strength/4, source_strength/4]};

% run the simulation
sensor_data = kspaceFirstOrder3DC(kgrid, medium, source, transducer, input_args{:});

 k = sensor_data.p_rms;
 save('transballsource0.2mpressure.mat','k');
% extract a single scan line from the sensor data using the current
% beamforming settings
scan_line = transducer.scan_line(sensor_data);

% =========================================================================
% VISUALISATION
% =========================================================================

% plot the recorded time series
figure;
stackedPlot(kgrid.t_array*1e6, sensor_data);
xlabel('Time [\mus]');
ylabel('Transducer Element');
title('Recorded Pressure');

% % plot the scan line
% figure;
% plot(kgrid.t_array*1e6, scan_line, 'k-');
% xlabel('Time [\mus]');
% ylabel('Pressure [au]');
% title('Scan Line After Beamforming');
% %% plot fft of scanline
% figure
% plot(abs(fft(scan_line)));
% figure
% [pxx,f] = pwelch(sensor_data(500,:),500,300,500,1e7);
% plot(f,10*log10(pxx));
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');


