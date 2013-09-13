function process_golgi_signal_nuclei(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.parse(exp_dir,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


start_direction = tic;
tracking_mat = csvread(fullfile(exp_dir,'tracking_mat.csv'));

tdTom_files = dir(fullfile(exp_dir,'tdTom'));
tdTom_files = tdTom_files(3:end);

props = struct('centroid_x', NaN*ones(size(tracking_mat)), ...
    'centroid_y', NaN*ones(size(tracking_mat)),...
    'cell_angle', NaN*ones(size(tracking_mat)),...
    'cell_angle_to_center', NaN*ones(size(tracking_mat)),...
    'cell_mag', NaN*ones(size(tracking_mat)),...
    'golgi_signal', NaN*ones(size(tracking_mat)),...
    'cell_close_to_wound', zeros(size(tracking_mat)));

golgi_files = dir(fullfile(exp_dir,'golgi'));
golgi_files = golgi_files(3:end);

labeled_nuclei_files = dir(fullfile(exp_dir,'labeled_nuclei'));
labeled_nuclei_files = labeled_nuclei_files(3:end);

cell_region_binary_files = dir(fullfile(exp_dir,'cell_region_binary'));
cell_region_binary_files = cell_region_binary_files(3:end);

direction_cmap = cool(180);

center_line_file = fullfile(exp_dir,'center_line.txt');
if (exist(center_line_file,'file') == 2)
   center_line = csvread(center_line_file);
else
    error('Need a defined center line file, ''center_line.txt''');
end

mkdir(fullfile(exp_dir,'golgi_direction'));
for i_num = 1:length(tdTom_files)
    tdTom = double(imread(fullfile(exp_dir,'tdTom',tdTom_files(i_num).name)));
    tdTom_norm = (tdTom - min(tdTom(:)))/range(tdTom(:));

    golgi = double(imread(fullfile(exp_dir,'golgi',golgi_files(i_num).name)));
    golgi_norm = (golgi - min(golgi(:)))/range(golgi(:));

    composite_norm = cat(3,tdTom_norm,golgi_norm,zeros(size(golgi)));

    %Center Line - Highlighting
    center_highlight = zeros(size(tdTom));
    center_highlight(center_line,:) = 1;
    composite_norm = create_highlighted_image(composite_norm,center_highlight,'color_map',[0,0,1]);
    
    nuclei_labeled = imread(fullfile(exp_dir,'labeled_nuclei',labeled_nuclei_files(i_num).name));
    
    %Adjacent to center region finding
    cell_region = imread(fullfile(exp_dir,'cell_region_binary',cell_region_binary_files(i_num).name));
    no_cell_region_labeled = bwlabel(not(cell_region));
    touching_center_nums = unique(no_cell_region_labeled(center_line,:));
    touching_center_nums = touching_center_nums(touching_center_nums ~= 0);
    center_region = ismember(no_cell_region_labeled,touching_center_nums);
    adjacent_to_center = imdilate(bwperim(center_region),strel('disk',40)) & not(center_region);
    
    close_nuclei_nums = sort(unique(nuclei_labeled(adjacent_to_center)));
    
    highlight_images = struct('direction',zeros(size(tdTom)),'nucleus',zeros(size(tdTom)),...
        'golgi_cent',zeros(size(tdTom)));
    
    for obj_num = 1:max(nuclei_labeled(:))
        tracking_col = tracking_mat(:,i_num);
        nucleus_tracking_row = find(tracking_col == obj_num);
        
        if (any(close_nuclei_nums == obj_num))
            props.cell_close_to_wound(nucleus_tracking_row,i_num) = 1;
        end
        
        this_nucleus = nuclei_labeled == obj_num;
        surrounding = imdilate(this_nucleus,strel('disk',10)) & not(this_nucleus);
        props.golgi_signal(nucleus_tracking_row,i_num) = mean(golgi(surrounding));

%         highlights = zeros(size(this_nucleus));
%         highlights(bwperim(this_nucleus)) = 1;
%         highlights(bwperim(surrounding)) = 2;
%         imshow(create_highlighted_image(golgi_norm,highlights))

        golgi_prop = regionprops(double(surrounding),golgi,'WeightedCentroid');
        nucleus_prop = regionprops(this_nucleus,'Centroid');

        direction = [nucleus_prop.Centroid(2) - golgi_prop.WeightedCentroid(2) ...
            golgi_prop.WeightedCentroid(1) - nucleus_prop.Centroid(1)];

        props.centroid_x(nucleus_tracking_row,i_num) = nucleus_prop.Centroid(1);
        props.centroid_y(nucleus_tracking_row,i_num) = nucleus_prop.Centroid(2);

        angle = atan2(direction(1),direction(2))*(180/pi);
        props.cell_angle(nucleus_tracking_row,i_num) = angle;
        props.cell_mag(nucleus_tracking_row,i_num) = sqrt(direction(1)^2 + direction(2)^2);

        if (nucleus_prop.Centroid(2) < center_line)
            direction_to_center = angle + 90;
        else
            direction_to_center = angle - 90;
        end
        if (direction_to_center > 180), direction_to_center = direction_to_center - 360; end
        if (direction_to_center < -180), direction_to_center = direction_to_center + 360; end
        direction_to_center = abs(direction_to_center);
        props.cell_angle_to_center(nucleus_tracking_row,i_num) = direction_to_center;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Image Highlights
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Nucleus - Direction Color Map
        round_direction = round(direction_to_center);
        if(round_direction == 0), round_direction = 1; end
        if(round_direction > 180), round_direction = 180; end
        highlight_images.direction(this_nucleus) = round_direction;

        %Nucleus - Border + Centroid Mark
        highlight_images.nucleus(bwperim(this_nucleus)) = 1;

        nucleus_mark = zeros(size(this_nucleus));
        nucleus_mark(round(nucleus_prop.Centroid(2)),round(nucleus_prop.Centroid(1))) = 1;
        nucleus_mark = imdilate(nucleus_mark,strel('disk',2));
        highlight_images.nucleus(logical(nucleus_mark)) = 1;

        %Golgi - Centroid Mark
        golgi_mark = zeros(size(this_nucleus));
        golgi_mark(round(golgi_prop.WeightedCentroid(2)),round(golgi_prop.WeightedCentroid(1))) = 1;
        golgi_mark = imdilate(golgi_mark,strel('disk',2));
        highlight_images.golgi_cent(logical(golgi_mark)) = 1;
    end
    
    composite_norm = create_highlighted_image(composite_norm,highlight_images.direction, ...
        'color_map',direction_cmap);
    composite_norm = create_highlighted_image(composite_norm,highlight_images.nucleus,'color_map',[1,1,1]);
    composite_norm = create_highlighted_image(composite_norm,highlight_images.golgi_cent,'color_map',[0,0,1]);
    composite_norm = create_highlighted_image(composite_norm,bwperim(adjacent_to_center), ...
        'color_map',[246,127,37]/255);
    imwrite(composite_norm,fullfile(exp_dir,'golgi_direction',sprintf('%03d.png',i_num)));
end

prop_dir = fullfile(exp_dir,'properties');
if (exist(prop_dir,'dir') == 0), mkdir(prop_dir); end

for prop_name = fieldnames(props)'
    csvwrite(fullfile(prop_dir,[prop_name{1},'.csv']),props.(prop_name{1}));  
end

toc(start_direction);