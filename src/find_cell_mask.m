function find_cell_mask(exp_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.parse(exp_dir);

start_boundary = tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vinc_files = dir(fullfile(exp_dir,'Vinculin'));
vinc_files = vinc_files(3:end);

mkdir(fullfile(exp_dir,'cell_region_highlight'));
mkdir(fullfile(exp_dir,'cell_region_binary'));

vinc_norm = cell(length(vinc_files),1);
vinc_med_filt = cell(length(vinc_files),1);
threshold_set = zeros(length(vinc_files),1);
for i = 1:length(vinc_files)
    vinc = double(imread(fullfile(exp_dir,'Vinculin',vinc_files(i).name)));
    vinc_norm{i} = (vinc - min(vinc(:)))/range(vinc(:));

    vinc_med_filt{i} = medfilt2(vinc,[11,11],'symmetric');
    
    threshold_set(i) = determine_threshold(vinc_med_filt{i}(:));
    
    if (mod(i,10) == 0)
        disp(['Done finding threshold in image number: ',num2str(i)]);
    end
end

for i = 1:length(vinc_files)
    background_region = vinc_med_filt{i} < threshold_set(i);
    background_region = imfill(background_region,'holes');

    cell_region = not(background_region);
    output_file = fullfile(exp_dir,'cell_region_binary',sprintf('%03d.png',i));
    imwrite(cell_region,output_file);

    cell_highlight = create_highlighted_image(vinc_norm{i},cell_region,'mix_percent',0.25);
    output_file = fullfile(exp_dir,'cell_region_highlight',sprintf('%03d.png',i));
    imwrite(cell_highlight,output_file);
end

fprintf('Done finding the cell regions, total time: %d, each image took: %d\n', ...
    toc(start_boundary),toc(start_boundary)/length(vinc_files))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function threshold = determine_threshold(pixels)

[heights, intensity] = hist(pixels,1000);

smoothed_heights = smooth(heights,0.05,'loess');
[~,imax,~,imin]= extrema(smoothed_heights);

%keep in mind that the zmax is sorted by value, so the highest peak is
%first and the corresponding index is also first in imax, the same pattern
%hold for zmin and imin

sorted_max_indexes = sort(imax);
first_max_index = find(sorted_max_indexes == imax(1));

%locate the index between the first two maximums
min_index = find(imin > sorted_max_indexes(first_max_index) & imin < sorted_max_indexes(first_max_index + 1));
assert(length(min_index) == 1, 'Error: expected to only find one minimum index between the first two max indexes, instead found %d', length(min_index));
threshold = intensity(imin(min_index));