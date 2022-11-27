function visualize_velocity_ratio(truck_data_filename,PIV_data_filename,output_folder)
    %
    % compute_velocity_ratio
    %
    % created Hibiya Haraki 2022
    % All risks of running this script is always with you.
    %

    % Input check
    if (~exist(truck_data_filename,'file'))
        error("Error: Cannot find %s",truck_data_filename);
        return;
    end

    if (~exist(PIV_data_filename,'file'))
        error("Error: Cannot find %s",PIV_data_filename);
        return;
    end

    if (~exist(output_folder,'dir'))
        error("Error: Cannot find %s",output_folder);
        return;
    end

    % Load data
    truck_data = load(truck_data_filename, ...
        'sfreq', ...
        'truck_acceleration_mean_vx', ...
        'determined_pixel_size','convert_point');
    PIV_data = load(PIV_data_filename,'meanValue_u_filtered');

    % Define time step
    time_step = min([length(truck_data.truck_acceleration_mean_vx),length(PIV_data.meanValue_u_filtered)]);

    % Compute velocity ratio
    reference_velocity_ratio = PIV_data.meanValue_u_filtered(1:time_step) ./ truck_data.truck_acceleration_mean_vx(1:time_step);

    % Visualize both velocity
    figure
    hold on
    plot((0:time_step-1)./truck_data.sfreq,truck_data.truck_acceleration_mean_vx(1:time_step));
    plot((0:time_step-1)./truck_data.sfreq,PIV_data.meanValue_u_filtered(1:time_step));
    hold off
    grid on
    legend("Truck x velocity","PIV flow velocity u");
    xlabel('t [s]');
    ylabel('Velocity ratio');
    title('Velocity (truck v_x & PIV u)');
    saveas(gcf,output_folder+"/"+"reference_velocity_all.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"reference_velocity_all.png"));

    % Visualize all
    figure
    plot((0:time_step-1)./truck_data.sfreq,reference_velocity_ratio);
    hold on
    for ii = 1:length(truck_data.convert_point)
        xline(truck_data.convert_point(ii),'k-');
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('Velocity ratio');
    xlim([0 (time_step-1)./truck_data.sfreq]);
    title('Velocity ratio (PIV u / truck v_x)');
    saveas(gcf,output_folder+"/"+"reference_velocity_ratio_all.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"reference_velocity_ratio_all.png"));

    % Visualize part
    figure
    plot((0:time_step-1)./truck_data.sfreq,reference_velocity_ratio);
    hold on
    for ii = 1:length(truck_data.convert_point)
        xline(truck_data.convert_point(ii)/truck_data.sfreq,'k-');
    end
    hold off
    grid on
    xlabel('t [s]');
    ylabel('Velocity ratio');
    xlim([truck_data.convert_point(2)./truck_data.sfreq (time_step-1)./truck_data.sfreq]);
    ylim([min(reference_velocity_ratio(truck_data.convert_point(2):end)) max(reference_velocity_ratio(truck_data.convert_point(2):end))]);
    title('Velocity ratio (PIV u / truck v_x)');
    saveas(gcf,output_folder+"/"+"reference_velocity_ratio_part.png");
    logging_func(sprintf("Save %s",output_folder+"/"+"reference_velocity_ratio_part.png"));

end