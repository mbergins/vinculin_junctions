function find_nuclei(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('terraced_thresh_debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

terraced_thresh_debug = i_p.Results.terraced_thresh_debug;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_nuc = tic;

tdTom_files = dir(fullfile(exp_dir,'tdTom'));
tdTom_files = tdTom_files(3:end);

golgi_files = dir(fullfile(exp_dir,'golgi'));
golgi_files = golgi_files(3:end);

cell_region_binary_files = dir(fullfile(exp_dir,'cell_region_binary'));
cell_region_binary_files = cell_region_binary_files(3:end);

mkdir(fullfile(exp_dir,'labeled_nuclei'));
mkdir(fullfile(exp_dir,'terraced_thresh_nuclei'));

if (terraced_thresh_debug)
    mkdir(fullfile(exp_dir,'terraced_thresh_stages'));
end

for i = 1:length(tdTom_files)
    tdTom = double(imread(fullfile(exp_dir,'tdTom',tdTom_files(i).name)));
    tdTom_norm = (tdTom - min(tdTom(:)))/range(tdTom(:));
    
    tdTom = medfilt2(tdTom,[5,5],'symmetric');
    
    golgi = double(imread(fullfile(exp_dir,'golgi',golgi_files(i).name)));
    golgi_labeled = bwlabel(imfill(golgi > mean(golgi(:)) + 2*std(golgi(:)),'holes'));
    golgi_props = regionprops(golgi_labeled,'Area');
    golgi_binary = ismember(golgi_labeled,find([golgi_props.Area] > 20 & [golgi_props.Area] < 200));
    golgi_perim = bwperim(golgi_binary);
    
    cell_region = imread(fullfile(exp_dir,'cell_region_binary',cell_region_binary_files(i).name));
    cell_region_perim = bwperim(cell_region);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Terraced Segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    net_seg = zeros(size(cell_region));
    net_vote = zeros(size(cell_region));
    upper_search_range = ceil(quantile(tdTom(:),0.99)/100)*100;
    for upper_thresh = 250:5:upper_search_range
        temp_orig_binary = tdTom < upper_thresh & cell_region;
        temp_orig_binary = remove_edge_objects(temp_orig_binary);
        
        temp_label = bwlabel(temp_orig_binary,4);
        
        props = regionprops(temp_label,'FilledArea','MajorAxisLength','MinorAxisLength');
        
        eclipse_area_ratio = [props.FilledArea]./(pi .* [props.MajorAxisLength]/2 .* [props.MinorAxisLength]/2);
        ratio = [props.MajorAxisLength]./[props.MinorAxisLength];
        temp_binary = ismember(temp_label, ...
            find(eclipse_area_ratio >= 0.95 & ratio <= 3 & [props.FilledArea] > 200 & [props.FilledArea] < 2500));
        temp_binary = imfill(temp_binary,'holes');
        
        net_seg = net_seg | temp_binary;
        net_vote = net_vote + temp_binary;
        
        if (terraced_thresh_debug)
            highlights = double(temp_orig_binary);
            highlights(net_seg > 0) = 2;
            highlights(cell_region_perim) = 3;
            highlights(golgi_perim) = 4;
            
            imwrite(create_highlighted_image(tdTom_norm,highlights,'mix_percent',0.25),...
                fullfile(exp_dir,'terraced_thresh_stages',sprintf('%05d.png',upper_thresh)));
        end
    end
    
    nuclei_binary = imfill(net_seg,'holes');
    nuclei_labeled = bwlabel(net_seg,4);

    nuclei_filt_highlight = create_highlighted_image(tdTom_norm,nuclei_binary,'mix_percent',0.25);
    output_file = fullfile(exp_dir,'terraced_thresh_nuclei',sprintf('%03d.png',i));
    
    imwrite(nuclei_filt_highlight,output_file);
    
    imwrite(uint16(nuclei_labeled),fullfile(exp_dir,'labeled_nuclei',sprintf('%03d.png',i)),'BitDepth',16);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % High Pass Segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     I_filt = fspecial('disk',21);
%     blurred_image = imfilter(tdTom,I_filt,'same','replicate');
% 
%     high_passed_tdTom = tdTom - blurred_image;
% 
%     nuclei = high_passed_tdTom < -.25*std(high_passed_tdTom(cell_region)) & cell_region;
% 
%     nuclei = imfill(nuclei,'holes');
% 
%     nuclei_labeled = bwlabel(nuclei,4);
%     props = regionprops(nuclei_labeled,'Area'); %#ok<MRPBW>
% 
%     nuclei = ismember(nuclei_labeled,find([props.Area] > 200 & [props.Area] < 5000));
%     nuclei_area_filt = nuclei;
% 
%     nuclei_highlight = create_highlighted_image(tdTom_norm,nuclei_area_filt,'mix_percent',0.25);
%     output_file = fullfile(exp_dir,'no_filt_nuclei',sprintf('%03d.png',i));
%     imwrite(nuclei_highlight,output_file);
% 
%     nuclei_labeled = bwlabel(nuclei,4);
%     props = regionprops(nuclei_labeled,'Eccentricity','Area','Solidity','BoundingBox'); %#ok<MRPBW>
% 
%     nuclei_nums = find([props.Solidity] > 0.9234 & [props.Area] > 200 & [props.Area] < 5000);
%     nuclei_filt = ismember(nuclei_labeled,nuclei_nums);
%     nuclei_labeled = bwlabel(nuclei_filt,4);
% 
%     nuclei_filt_highlight = create_highlighted_image(tdTom_norm,nuclei_filt,'mix_percent',0.25);
%     output_file = fullfile(exp_dir,'filt_nuclei',sprintf('%03d.png',i));
% 
%     imwrite(nuclei_filt_highlight,output_file);
%     imwrite(uint16(nuclei_labeled),fullfile(exp_dir,'labeled_nuclei',sprintf('%03d.png',i)),'BitDepth',16);

%     differential_highlight = zeros(size(nuclei));
%     differential_highlight(nuclei_area_filt) = 1;
%     differential_highlight(nuclei_labeled > 0) = 2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Manual Classification Code
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     props_mat = zeros(max(nuclei_labeled(:)),length(fieldnames(props)) + 1 - 1);
%     props_mat(:,2) = [props.Eccentricity]';
%     props_mat(:,3) = [props.Area]';
%     props_mat(:,4) = [props.Solidity]';
%     csvwrite('obj_props.csv',props_mat);
%
%     for obj_num = 1:max(nuclei_labeled(:))
%         this_obj = nuclei_labeled == obj_num;
%         this_perim = bwperim(this_obj);
%
%         bbox = round(props(obj_num).BoundingBox);
%         rows = [bbox(2)-20000,bbox(2)+bbox(4)+20000];
%         if (rows(1) < 1), rows(1) = 1; end
%         if (rows(2) > size(this_obj,1)), rows(2) = size(this_obj,1); end
%         cols = [bbox(1)-200000,bbox(1)+bbox(3)+200000];
%         if (cols(1) < 1), cols(1) = 1; end
%         if (cols(2) > size(this_obj,2)), cols(2) = size(this_obj,2); end
%
%         this_region = tdTom_norm(rows(1):rows(2),cols(1):cols(2));
%         this_region_perim = this_perim(rows(1):rows(2),cols(1):cols(2));
%         this_highlight = create_highlighted_image(this_region,this_region_perim,'mix_percent',0.5);
%         imwrite(this_highlight,fullfile('single_objs',sprintf('%04d.png',obj_num)));
%     end
end
toc(start_nuc);