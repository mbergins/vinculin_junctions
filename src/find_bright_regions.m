function find_bright_regions(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.parse(exp_dir,varargin{:});

min_size = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_nuc = tic;

vinc_files = dir(fullfile(exp_dir,'Vinculin'));
vinc_files = vinc_files(3:end);

mkdir_no_err(fullfile(exp_dir,'bright_regions_highlight'));
mkdir_no_err(fullfile(exp_dir,'bright_regions'));

for i = 1:length(vinc_files)
    vinc = double(imread(fullfile(exp_dir,'Vinculin',vinc_files(i).name)));
    vinc_ranges = quantile(vinc(:),[0.01,0.99]);
    vinc_norm = (vinc - vinc_ranges(1))/vinc_ranges(2);
    
    I_filt = fspecial('disk',11);
    vinc_blur = imfilter(vinc,I_filt,'replicate');
    high_passed_image = vinc - vinc_blur;
    
    bright_regions = high_passed_image >= mean(high_passed_image(:)) + ...
        2*std(high_passed_image(:));
    bright_regions = imfill(bright_regions,'holes');
    bright_regions = bwareaopen(bright_regions,min_size);
    
    bright_highlight = create_highlighted_image(vinc_norm,bright_regions, ...
        'mix_percent',0.5,'color_map',[0,1,0]);
    
    imwrite(bright_regions,fullfile(exp_dir,'bright_regions',sprintf('%03d.png',i)));
    imwrite(bright_highlight,fullfile(exp_dir,'bright_regions_highlight',sprintf('%03d.png',i)));
end
toc(start_nuc);