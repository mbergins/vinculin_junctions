function correlate_golgi_movement(exp_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.parse(exp_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cell_angle = csvread(fullfile(exp_dir,'properties','cell_angle.csv'));
cell_mag = csvread(fullfile(exp_dir,'properties','cell_mag.csv'));

centroid_x = csvread(fullfile(exp_dir,'properties','centroid_x.csv'));
centroid_y = csvread(fullfile(exp_dir,'properties','centroid_y.csv'));

predicted_angle = [];
actual_angle = [];
this_magnitude = [];

for cell_num = 1:size(cell_angle,1)
    if (sum(not(isnan(cell_angle(cell_num,:)))) > 1)
        this_cell_angle = cell_angle(cell_num,:);
        this_cent_x = centroid_x(cell_num,:);
        this_cent_y = centroid_y(cell_num,:);
        this_cell_mag = cell_mag(cell_num,:);
        
        for i_num = 1:length(this_cell_angle)
            if (isnan(this_cell_angle(i_num)))
                continue;
            end
            
            if (i_num == length(this_cell_angle))
                continue;
            end
            
            if (isnan(this_cent_x(i_num + 1)))
                continue;
            end
            
            x_movement = this_cent_x(i_num + 1) - this_cent_x(i_num);
            y_movement = this_cent_y(i_num) - this_cent_y(i_num + 1);
            
            actual_angle = [actual_angle, atan2(y_movement,x_movement)*(180/pi)];
            predicted_angle = [predicted_angle,this_cell_angle(i_num)];
            this_magnitude = [this_magnitude, this_cell_mag(i_num)];
        end
    end
end

plot(actual_angle,predicted_angle,'o');
xlabel('Actual Movement Angle')
ylabel('Golgi Angle');
saveas(gcf,fullfile(exp_dir,'golgi_vs_actual_angle.png'))
close gcf;

min_mag_vals = 0:15;
correlation = [];

for min_mag = min_mag_vals
    correlation = [correlation, ...
        corr(actual_angle(this_magnitude >= min_mag)',predicted_angle(this_magnitude >= min_mag)')];
end

plot(min_mag_vals,correlation,'o');
xlabel('Min Magnitude of Golgi Direction (pixels)')
ylabel('Correlation');
ylim([0,1]);
saveas(gcf,fullfile(exp_dir,'min_mag_vs_correlation.png'))
close gcf;

hist(this_magnitude);
xlabel('Magnitude of Golgi Direction (pixels)')
saveas(gcf,fullfile(exp_dir,'golgi_direction_magnitude.png'))
close gcf;