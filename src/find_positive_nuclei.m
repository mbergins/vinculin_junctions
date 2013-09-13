function find_positive_nuclei(exp_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.parse(exp_dir,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_nuc = tic;

tdTom_files = dir(fullfile(exp_dir,'tdTom'));
tdTom_files = tdTom_files(3:end);

mkdir(fullfile(exp_dir,'labeled_nuclei'));
mkdir(fullfile(exp_dir,'nuclei_highlight'));

for i = 1:length(tdTom_files)
    tdTom = double(imread(fullfile(exp_dir,'tdTom',tdTom_files(i).name)));
    tdTom_norm = (tdTom - min(tdTom(:)))/range(tdTom(:));
    
    tdTom = medfilt2(tdTom,[5,5],'symmetric');
    
    nuclei = tdTom > 300;
    nuclei_label = bwlabel(nuclei,4);
    
    props = regionprops(nuclei_label,'FilledArea','MajorAxisLength','MinorAxisLength');

    eclipse_area_ratio = [props.FilledArea]./(pi .* [props.MajorAxisLength]/2 .* [props.MinorAxisLength]/2);
    ratio = [props.MajorAxisLength]./[props.MinorAxisLength];
    nuclei = ismember(nuclei_label, ...
        find(eclipse_area_ratio >= 0.95 & ratio <= 3 & [props.FilledArea] > 200));
%     nuclei = ismember(nuclei_label, ...
%         find([props.FilledArea] > 200));
    nuclei = imfill(nuclei,'holes');
    nuclei_labeled = bwlabel(nuclei,4);
    
    highlighted_image = create_highlighted_image(tdTom_norm, bwperim(nuclei));
    imwrite(highlighted_image,fullfile(exp_dir,'nuclei_highlight',sprintf('%03d.png',i)));
    
    imwrite(uint16(nuclei_labeled),fullfile(exp_dir,'labeled_nuclei',sprintf('%03d.png',i)),'BitDepth',16);
end
toc(start_nuc);