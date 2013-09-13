function track_nuclei(exp_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.parse(exp_dir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_track = tic;
labeled_nuclei_files = dir(fullfile(exp_dir,'labeled_nuclei'));
labeled_nuclei_files = labeled_nuclei_files(3:end);
all_tracking_props = cell(length(labeled_nuclei_files),1);

for i_num = 1:length(labeled_nuclei_files)
    nuclei_labeled = imread(fullfile(exp_dir,'labeled_nuclei',labeled_nuclei_files(i_num).name));
    current_props = regionprops(nuclei_labeled,'Area');

    %if we are looking at the last image, there isn't a next image to
    %compare to, so we save the areas, so the tracking algorithm knows to
    %put some placeholders in, then cut out of the loop
    if (i_num == length(labeled_nuclei_files))
        all_tracking_props{i_num} = current_props;
        continue;
    end

    %we know there is another image ahead of this one, calculate pixel
    %similarity for the next image and this one
    nuclei_labeled_next = imread(fullfile(exp_dir,'labeled_nuclei',labeled_nuclei_files(i_num+1).name));

    for obj_num = 1:max(nuclei_labeled(:))
        this_nucleus = nuclei_labeled == obj_num;

        overlap_pix = nuclei_labeled_next(this_nucleus);
        overlap_pix = overlap_pix(overlap_pix ~= 0);

        current_props(obj_num).Pix_sim = zeros(1,max(nuclei_labeled_next(:)));

        for overlap_obj = unique(overlap_pix)'
            current_props(obj_num).Pix_sim(overlap_obj) = sum(overlap_pix == overlap_obj)/current_props(obj_num).Area;
        end
    end
    all_tracking_props{i_num} = current_props;
end
save(fullfile(exp_dir,'all_tracking_props.mat'),'all_tracking_props');

track_objects(exp_dir);

toc(start_track);