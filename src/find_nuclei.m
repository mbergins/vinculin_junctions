function find_nuclei(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('terraced_thresh_debug',0,@(x)x == 1 || x == 0);

i_p.addParamValue('min_nuclei_size',500,@(x)x == 1 || x == 0);
i_p.addParamValue('max_nuclei_size',20000,@(x)x == 1 || x == 0);
i_p.addParamValue('min_nuclei_intensity',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

terraced_thresh_debug = i_p.Results.terraced_thresh_debug;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_nuc = tic;

vinc_files = dir(fullfile(exp_dir,'Vinculin'));
vinc_files = vinc_files(3:end);

mkdir_no_err(fullfile(exp_dir,'labeled_nuclei'));
mkdir_no_err(fullfile(exp_dir,'terraced_thresh_nuclei'));

if (terraced_thresh_debug)
    mkdir_no_err(fullfile(exp_dir,'terraced_thresh_stages'));
end

for i = 1:length(vinc_files)
    vinc = double(imread(fullfile(exp_dir,'Vinculin',vinc_files(i).name)));
    vinc_orig = vinc;
    vinc_norm = (vinc - min(vinc(:)))/range(vinc(:));
    
    vinc = medfilt2(vinc,[5,5],'symmetric');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Terraced Segmentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    net_seg = zeros(size(vinc));
    net_vote = zeros(size(vinc));
    search_range = linspace(quantile(vinc(:),0.01),quantile(vinc(:),0.99));
    for upper_thresh = search_range
        temp_orig_binary = vinc < upper_thresh;
        temp_orig_binary = remove_edge_objects(temp_orig_binary);
        
        %Quick filter to toss out small/large objects
        temp_label = bwlabel(temp_orig_binary,4);
        props = regionprops(temp_label,'Area'); %#ok<MRPBW>
        temp_orig_binary = ismember(temp_label,find([props.Area] > 200 & [props.Area] < 20000));
        temp_label = bwlabel(temp_orig_binary,4);
        
        %Better filters to pick out appropriate size/shaped objects
        props = regionprops(temp_label,'FilledArea','MajorAxisLength','MinorAxisLength'); %#ok<MRPBW>
        eclipse_area_ratio = [props.FilledArea]./(pi .* [props.MajorAxisLength]/2 .* [props.MinorAxisLength]/2);
        ratio = [props.MajorAxisLength]./[props.MinorAxisLength];
        temp_binary = ismember(temp_label, ...
            find(eclipse_area_ratio >= 0.95 & ratio <= 3 & ...
            [props.FilledArea] > i_p.Results.min_nuclei_size & ...
            [props.FilledArea] < i_p.Results.max_nuclei_size));
        
        temp_binary = imfill(temp_binary,'holes');
        
        net_seg = net_seg | temp_binary;
        net_vote = net_vote + temp_binary;
        
        if (terraced_thresh_debug)
            highlights = double(temp_orig_binary);
            highlights(net_seg > 0) = 2;
            
            target_folder = fullfile(exp_dir,'terraced_thresh_stages',sprintf('%02d',i));
            if (not(exist(target_folder,'dir')))
                mkdir(target_folder);
            end
            highlight_image = create_highlighted_image(vinc_norm,highlights,'mix_percent',0.25);
            highlight_image = imresize(highlight_image,[512,NaN]);
            imwrite(highlight_image,...
                fullfile(target_folder,sprintf('%010d.png',round(upper_thresh))));
        end
    end
    
    nuclei_binary = imfill(net_seg,'holes');
    nuclei_labeled = bwlabel(nuclei_binary,4);
    
    props = regionprops(nuclei_labeled,vinc_orig,'MeanIntensity'); %#ok<MRPBW>
    nuclei_binary = ismember(nuclei_labeled, ...
        find([props.MeanIntensity] >= i_p.Results.min_nuclei_intensity));
    nuclei_labeled = bwlabel(nuclei_binary,4);
    
    nuclei_filt_highlight = create_highlighted_image(vinc_norm,nuclei_binary,'mix_percent',0.25);
    output_file = fullfile(exp_dir,'terraced_thresh_nuclei',sprintf('%03d.png',i));
    
    imwrite(nuclei_filt_highlight,output_file);
    
    imwrite(uint16(nuclei_labeled),...
        fullfile(exp_dir,'labeled_nuclei',sprintf('%03d.png',i)),'BitDepth',16);
end
toc(start_nuc);