function track_objects(exp_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(fullfile(exp_dir,'all_tracking_props.mat'),'all_tracking_props');

%fill in place holder values
for i_num=1:length(all_tracking_props) %#ok<NODEF>
    area = [all_tracking_props{i_num}.Area];
    
    %this variable will be used when filling out the tracking matrix
    for i=1:size(area,2)
        all_tracking_props{i_num}(i).assigned = 0; %#ok<AGROW>
        all_tracking_props{i_num}(i).next_obj = []; %#ok<AGROW>
    end
end

%make tracking decisions based on pixel similarity
for i_num=1:(length(all_tracking_props)-1)
    %if Pix_sim isn't present there either aren't any cells in the current
    %image or there aren't any in the next image, either way, we don't need
    %to do any tracking, jump to next image
    if (not(any(strcmp(fieldnames(all_tracking_props{i_num}),'Pix_sim'))))
        continue;
    end
    
    %this reshape gives the pixel similarity with all the starting cells in
    %the rows and ending cells in the columns
    pix_sim = reshape([all_tracking_props{i_num}.Pix_sim],[],length(all_tracking_props{i_num}))';
    
    %Start by searching for reciprical high pixel similarity matches,
    %defined as those cells that overlap a single cell in the next frame by
    %50% or more
    [start_cell_hits, end_cell_hits] = find(pix_sim > 0.5);    
    for i = 1:length(start_cell_hits)
        start_cell = start_cell_hits(i);
        end_cell = end_cell_hits(i);
        
        %check to make sure these two cells are unique in their lists, if
        %so, make the connection and block the row (start cell) and column
        %(end cell)
        if (sum(start_cell == start_cell_hits) == 1 && ...
            sum(end_cell == end_cell_hits) == 1)
            
            all_tracking_props{i_num}(start_cell).next_obj = end_cell; %#ok<AGROW>
            
            pix_sim(start_cell,:) = NaN;
            pix_sim(:,end_cell) = NaN;
        end
    end
    
    %This loop finds any remaining pixel similarity measures above 20% and
    %matches them, one-by-one to their corresponding end cell. This code
    %mostly helps clear out the cells that overlap two objects with fairly
    %high similarity, but not one single object completely.
    while (any(any(pix_sim > 0.2)))
        [start_cell,end_cell] = find(pix_sim == max(pix_sim(:)),1,'first');
        
        all_tracking_props{i_num}(start_cell).next_obj = end_cell; %#ok<AGROW>
        
        pix_sim(start_cell,:) = NaN;
        pix_sim(:,end_cell) = NaN;
    end
    
%     disp(['Done tracking in image number: ',i_num,'/',(length(all_tracking_props)-1)]);
end

tracking_mat = convert_tracking_props_to_matrix(all_tracking_props);

output_file = fullfile(exp_dir,'tracking_mat.csv');

%If the tracking matrix is empty, don't output anything and don't make a
%folder for the empty matrix
if (not(any(size(tracking_mat) == 0)))
    %check each column of the tracking matrix to make sure each cell number
    %occurs once and only once per column
    for col_num = 1:size(tracking_mat,2)
        all_tracking_nums = sort(tracking_mat(:,col_num));
        for cell_num = 1:max(all_tracking_nums)
            assert(sum(all_tracking_nums == cell_num) == 1);
        end
    end
    
    if (not(exist(fileparts(output_file),'dir'))), mkdir(fileparts(output_file)); end
    
    csvwrite(output_file,tracking_mat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tracking_matrix = convert_tracking_props_to_matrix(all_tracking_props)

%this structure will contain all the individual sequences of tracking
%numbers for each cell
cells = struct('start',{},'sequence',{});

%each tracked cell will be stored in a single row in the tracking matrix,
%this value will hold the current tracking number for this cell
tracking_num = 1;

%we start by initializing the image number and the object number to the
%first object in the first frame
[i_num,obj_num] = find_unassigned_obj(all_tracking_props);

%i_num of 0 and obj_num of 0 will indicate that there aren't any more
%objects to follow
while (i_num ~= 0 && obj_num ~= 0)
    %if tracking_num has been incremented, we are looking at a new cell, so
    %we initialize the start value
    if (length(cells) < tracking_num)
        cells(tracking_num).start = i_num;
    end
    
    %add the current object number to the sequence for the tracking num
    cells(tracking_num).sequence = [cells(tracking_num).sequence, obj_num];
    %tell future interations that this object is assigned
    all_tracking_props{i_num}(obj_num).assigned = 1;
    
    %pick out the next object to follow
    [i_num, obj_num] = follow_to_next_obj(all_tracking_props,i_num,obj_num);
    
    %an i_num of 0 indicates that no next object was found, see if there
    %are any unassigned objects left
    if (i_num == 0)
        [i_num, obj_num] = find_unassigned_obj(all_tracking_props);
        tracking_num = tracking_num + 1;
    end
end

tracking_matrix = zeros(length(cells),length(all_tracking_props));

for cell_num = 1:length(cells)
    col_range = cells(cell_num).start:(cells(cell_num).start + length(cells(cell_num).sequence) - 1);
    tracking_matrix(cell_num,col_range) = cells(cell_num).sequence;
end

for col_num = 1:size(tracking_matrix,2)
    this_col = tracking_matrix(:,col_num);
    this_col = this_col(this_col ~= 0);
    
    this_col = sort(unique(this_col))';
    
    %empty columns mean there weren't any adhesions in that time step, so
    %check for that in the following assert first, then if there were
    %cells, make sure all were accounted for
    assert(isempty(this_col) || all(this_col == 1:max(this_col)))
    
    assert((isempty(all_tracking_props{col_num}) && isempty(this_col)) || ...
        (length(all_tracking_props{col_num}) == max(this_col)));
end

function [i_num,obj_num] = find_unassigned_obj(all_tracking_props)

%this is a bit nasty, but basically is searches through the assigned
%variable in each image number from 1-last image, seeing if there are any
%unassigned objects, if we hit an unassigned object, we immediately return
for i_num=1:length(all_tracking_props)
    for obj_num=1:size(all_tracking_props{i_num},1)
        try
            %this combinationof image number of object number is already
            %assigned, continue searching
            if (all_tracking_props{i_num}(obj_num).assigned)
                continue;
            else
                return;
            end
        catch
            continue;
        end
    end
end

i_num = 0;
obj_num = 0;

return;

function [i_num,obj_num] = follow_to_next_obj(all_tracking_props,i_num,obj_num)

%if next object is empty, set i_num and obj_num to zero to indicate that we
%need to search for a new object follow, otherwise, increment i_num and
%update obj_num to the next value
try %#ok<*TRYNC>
    if (size(all_tracking_props{i_num}(obj_num).next_obj,2) == 0)
        i_num = 0;
        obj_num = 0;
        return;
    else
        obj_num = all_tracking_props{i_num}(obj_num).next_obj;
        i_num = i_num + 1;
        return;
    end
end

i_num = 0;
obj_num = 0;

return;