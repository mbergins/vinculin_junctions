function find_positive_nuclei_terrace(exp_dir,varargin)

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

mkdir(fullfile(exp_dir,'labeled_nuclei'));
mkdir(fullfile(exp_dir,'nuclei_highlight'));

if (terraced_thresh_debug)
    mkdir(fullfile(exp_dir,'terraced_thresh_stages'));
end

pix_vals = [];
for i = 1:length(tdTom_files)
    pix_vals = [pix_vals,double(imread(fullfile(exp_dir,'tdTom',tdTom_files(i).name)))]; %#ok<AGROW>
end
search_vals = linspace(quantile(pix_vals(:),0.7),quantile(pix_vals(:),0.99),100);
clear pix_vals;

for i = 1:length(tdTom_files)
    tdTom = double(imread(fullfile(exp_dir,'tdTom',tdTom_files(i).name)));
    tdTom_norm = (tdTom - min(tdTom(:)))/range(tdTom(:));
    
    tdTom = medfilt2(tdTom,[5,5],'symmetric');
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Terraced Segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    net_seg = zeros(size(tdTom));
    net_vote = zeros(size(tdTom));
    
    for thresh_index = 1:length(search_vals)
        temp_orig_binary = tdTom > search_vals(thresh_index);
        temp_orig_binary = remove_edge_objects(temp_orig_binary);

        temp_label = bwlabel(temp_orig_binary,4);
        
        props = regionprops(temp_label,'FilledArea','MajorAxisLength','MinorAxisLength');
        
        eclipse_area_ratio = [props.FilledArea]./(pi .* [props.MajorAxisLength]/2 .* [props.MinorAxisLength]/2);
        ratio = [props.MajorAxisLength]./[props.MinorAxisLength];
        temp_binary = ismember(temp_label, ...
            find(eclipse_area_ratio >= 0.95 & ratio <= 3 & [props.FilledArea] > 500 & [props.FilledArea] < 2500));
        temp_binary = imfill(temp_binary,'holes');
        
        net_seg = net_seg | temp_binary;
        net_vote = net_vote + temp_binary;
        
        if (terraced_thresh_debug)
            highlights = double(temp_orig_binary);
            highlights(temp_orig_binary) = 1;
            highlights(net_seg > 0) = 2;
            
            cmap = [1,0,0;0,1,0];
            
            imwrite(create_highlighted_image(tdTom_norm,highlights,'mix_percent',0.25,'color_map',cmap),...
                fullfile(exp_dir,'terraced_thresh_stages',sprintf('%05d.png',round(search_vals(thresh_index)))));
        end
    end

    nuclei_binary = imfill(net_seg,'holes');
    nuclei_labeled = bwlabel(net_seg,4);

    nuclei_filt_highlight = create_highlighted_image(tdTom_norm,nuclei_binary,'mix_percent',0.25);
    output_file = fullfile(exp_dir,'nuclei_highlight',sprintf('%03d.png',i));
    
    imwrite(nuclei_filt_highlight,output_file);
    
    imwrite(uint16(nuclei_labeled),fullfile(exp_dir,'labeled_nuclei',sprintf('%03d.png',i)),'BitDepth',16);
end
toc(start_nuc);