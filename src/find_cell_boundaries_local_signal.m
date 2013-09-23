function find_cell_boundaries_local_signal(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('no_cell_threshold',2000,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_nuc = tic;

vinc_files = dir(fullfile(exp_dir,'Vinculin'));
vinc_files = vinc_files(3:end);

nuclei_label_files = dir(fullfile(exp_dir,'labeled_nuclei'));
nuclei_label_files = nuclei_label_files(3:end);

for i = 1:length(vinc_files)
    vinc = double(imread(fullfile(exp_dir,'Vinculin',vinc_files(i).name)));
    vinc_norm = (vinc - min(vinc(:)))/range(vinc(:));
    
    nuclei_label = imread(fullfile(exp_dir,'labeled_nuclei',nuclei_label_files(i).name));
    vinc_med = medfilt2(vinc,[7,7],'symmetric');
    
    for nuc_num = 1:max(nuclei_label(:))
        diag_highlights = zeros(size(vinc));
        this_nuc = nuclei_label == nuc_num;
        not_this_nuc = nuclei_label ~= nuc_num;
        not_this_nuc(nuclei_label == 0) = 0;
        
        diag_highlights(this_nuc) = 1;
        diag_highlights(not_this_nuc) = 2;
        this_nuc_border = imdilate(this_nuc,strel('disk',3)) & not(nuclei_label);
        diag_highlights(this_nuc_border) = 3;
        
        mean_val_in_border = mean(vinc_med(this_nuc_border));
        
        within_thresh = vinc_med < 1.1*mean_val_in_border & vinc_med > 0.9*mean_val_in_border;
        diag_highlights(within_thresh & not(this_nuc) & not(this_nuc_border)) = 4;
        
        region_high = create_highlighted_image(vinc_norm,diag_highlights,'mix_percent',0.25);
        
        1;
    end
    
end
toc(start_nuc);