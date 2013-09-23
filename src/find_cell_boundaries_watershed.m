function find_cell_boundaries_watershed(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('no_cell_threshold',0,@(x)x == 1 || x == 0);

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
    
    %Marker based watershed
    vinc_watershed = imimposemin(vinc_med,nuclei_label);
    
    cell_bodies = watershed(vinc_watershed);    
    cell_edges = cell_bodies == 0;
    cell_bodies_label = bwlabel(cell_bodies);
    
    thick_cell_edges = imdilate(cell_edges,strel('disk',3));
    
    %Block out Regions not inside cells
    no_cells = vinc < i_p.Results.no_cell_threshold;
    no_cells = imfill(no_cells,'holes');
    
    results_highlight = zeros(size(vinc));
    results_highlight(no_cells) = 3;
    results_highlight(thick_cell_edges) = 1;
    results_highlight(nuclei_label > 0) = 2;
    
    
    cell_edge_highlight = create_highlighted_image(vinc_norm,results_highlight,'mix_percent',0.25);
    mkdir_no_err(fullfile(exp_dir,'cell_edge_highlight'));
    imwrite(cell_edge_highlight,...
        fullfile(exp_dir,'cell_edge_highlight',sprintf('%03d.png',i)));
    
    cell_body_highlight = create_highlighted_image(vinc_norm,cell_bodies_label,'mix_percent',0.25);
    
    
    
    
    1;
end
toc(start_nuc);