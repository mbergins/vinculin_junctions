function analyze_cell_polarity_results(exp_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.parse(exp_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_angle_to_center = csvread(fullfile(exp_dir,'properties','cell_angle_to_center.csv'));
cell_mag = csvread(fullfile(exp_dir,'properties','cell_mag.csv'));
cell_close_to_wound = csvread(fullfile(exp_dir,'properties','cell_close_to_wound.csv'));

cell_angle_high_mag = cell_angle_to_center;
cell_angle_high_mag(cell_mag < 6) = NaN;

high_mag = struct('mean',zeros(size(cell_angle_to_center,2),1),...
    'upper',zeros(size(cell_angle_to_center,2),1),...
    'lower', zeros(size(cell_angle_to_center,2),1),...
    'measurements',0);
all_mag = high_mag;
close_to_wound = high_mag;

for i_num = 1:size(cell_angle_to_center,2)
    high_mag.mean(i_num) = nanmean(cell_angle_high_mag(:,i_num));
    [~,~,this_ci] = ttest(cell_angle_high_mag(:,i_num));
    high_mag.upper(i_num) = this_ci(2);
    high_mag.lower(i_num) = this_ci(1);
    high_mag.measurements = high_mag.measurements + ...
        sum(not(isnan(cell_angle_high_mag(:,i_num))));
        
    all_mag.mean(i_num) = nanmean(cell_angle_to_center(:,i_num));
    [~,~,this_ci] = ttest(cell_angle_to_center(:,i_num));
    all_mag.upper(i_num) = this_ci(2);
    all_mag.lower(i_num) = this_ci(1);
    all_mag.measurements = all_mag.measurements + ...
        sum(not(isnan(cell_angle_to_center(:,i_num))));
    
    close_angles = cell_angle_to_center(cell_close_to_wound(:,i_num) == 1,i_num);
    
    close_to_wound.mean(i_num) = nanmean(close_angles);
%     [~,~,this_ci] = ttest(cell_angle_to_center(:,i_num));
%     close_to_wound.upper(i_num) = this_ci(2);
%     close_to_wound.lower(i_num) = this_ci(1);
%     close_to_wound.measurements = close_to_wound.measurements + ...
%         sum(not(isnan(cell_angle_to_center(:,i_num))));
end

plot(all_mag.mean,'LineWidth',2); hold on;
plot(all_mag.upper,'r','LineWidth',2); plot(all_mag.lower,'r','LineWidth',2); 
ylim([0,180]); 
xlabel('Image Number','FontSize',16); 
ylabel('Average Cell Direction Deviance From Center','FontSize',16);
plot([1 97],[90 90],'g');
plot(close_to_wound.mean,'k','LineWidth',2);
hold off;
saveas(gcf,fullfile(exp_dir,'all_angles_from center.png'))
close gcf;

% plot(high_mag.mean); hold on; plot(high_mag.upper,'r'); plot(high_mag.lower,'r'); 
% ylim([0,180]); xlabel('Image Number'); ylabel('Average Cell Direction Deviance');
% plot([1 97],[90 90],'g');
% hold off;