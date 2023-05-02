function [xc_1,xc_2] = FindTwoAggregates(xc,NC,internal_stresses,indices_of_internal_stresses,internal_faces_and_cubes_index_array,internal_faces_and_cubes_index_array_no_double_counting)

% for now we assume that the aggregate will break only in two parts
% based on which face is subject to the max internal stress
max_stress = 0;
for nn = 1:length(internal_stresses)
    max_found = max(norm(internal_stresses(:,nn)));
    if max_found>max_stress
        max_stress = max_found;
        index = nn;
        breaking_face = indices_of_internal_stresses(index);
    end
end

index_1 = find(internal_faces_and_cubes_index_array(:,1)==breaking_face);

cubes_to_break = internal_faces_and_cubes_index_array(index_1,2:3);

index_2 = find(internal_faces_and_cubes_index_array(:,2)==cubes_to_break(2))

% get index of which cubes break in the unrepeated array
for ii=1:length(internal_faces_and_cubes_index_array_no_double_counting)
    cube_bool = (internal_faces_and_cubes_index_array_no_double_counting(ii,2:3)==cubes_to_break(1:2));
    if sum(cube_bool)==2
        where = ii
        break;
    end
end

cube_copy = zeros(length(internal_faces_and_cubes_index_array_no_double_counting),1);
cube_copy(where) = -1; % this is where they break, so it should be discounted from the two aggregates built

% use second cube (this is the greater number by default because of sorting of second column)
first_instance = find(internal_faces_and_cubes_index_array_no_double_counting(:,2)==cubes_to_break(2));

% need special case if we have that the second cube is the only one that breaks from original aggregate
% this was an edge case found for NC=20 and Seed=30
if isempty(first_instance)==1
    first_cube_instead = find(internal_faces_and_cubes_index_array_no_double_counting(:,2)==cubes_to_break(1));
    cube_copy(first_cube_instead(:)) = NC+1;
    
    cubes_to_break = internal_faces_and_cubes_index_array(index_1,3);

    cubes_all = [internal_faces_and_cubes_index_array_no_double_counting cube_copy];

    % okay this seems to work now
    % so what is left to do is to store the indices of the NC+1 in one array
    % and the 0 in another one
    % those will make the two aggregates
    indices_NC = find(cubes_all(:,4)==(NC+1));
    cubes_all_xc_1 = cubes_all(indices_NC,3);
    cubes_all_xc_1 = cubes_all_xc_1(:);
    cubes_all_xc_1 = sort(cubes_all_xc_1);

    agg_idx_list = [];
    agg_idx_list(1) = cubes_all_xc_1(1);
    for kk=2:length(cubes_all_xc_1)
        if cubes_all_xc_1(kk)~=cubes_all_xc_1(kk-1)
            agg_idx_list = [agg_idx_list;cubes_all_xc_1(kk)];
        end
    end
    agg_idx_list;
    agg_idx_list = sort(agg_idx_list);
    agg_idx_list;

    xc;
    xc_1 = xc(agg_idx_list,:);
    xc_2 = xc;
    xc_2(agg_idx_list,:) = [];
end

% here we are away from the edge case, so the second cube has more than one
% cube connected to it
if isempty(first_instance)==0
    cube_copy(first_instance(:)) = NC+1;

    first_instance_old = first_instance;

    % the idea is to have two nested y loops to go up toward beginning of array, and down 1
    % until the end of the array to check all possible links
    % it is similar to a two-pointer method
    
    % jj will go down the array, and so it will increase by 1 each time
    jj = first_instance(end)+1;
    counter = length(internal_faces_and_cubes_index_array_no_double_counting);
    while jj<=length(internal_faces_and_cubes_index_array_no_double_counting)
        % ii will go up the array and so it will decrease by 1 each time
        % note that after each jj iteration, ii needs to be reset
        ii = first_instance_old(end);
        while ii>1
            % edge case
            if internal_faces_and_cubes_index_array_no_double_counting(ii-1,3)==internal_faces_and_cubes_index_array_no_double_counting(ii,3) && cube_copy(ii-1)==(NC+1)
               cube_copy(ii-1) = NC+1;
               disp("Case 1")
            end
            if internal_faces_and_cubes_index_array_no_double_counting(jj,2)==internal_faces_and_cubes_index_array_no_double_counting(jj-1,2) && cube_copy(jj-1)==(NC+1)  && cube_copy(jj)==0
               cube_copy(jj) = NC+1;
               disp("Case 2")
            end
            if internal_faces_and_cubes_index_array_no_double_counting(jj,2)==internal_faces_and_cubes_index_array_no_double_counting(ii,3) && cube_copy(ii)==(NC+1)
               cube_copy(jj) = NC+1;
               disp("Case 3")
            end
            first_instance(end) = first_instance(end)+1;
           
           % update starting index of ii index, as long as we are not at
           % as long as it is not the last element of the array
           if first_instance(end)<length(internal_faces_and_cubes_index_array_no_double_counting)
              first_instance_old(end) = first_instance(end);
           end
           
           ii = ii - 1;
        end
        % possible edge case (unneccesary probably now that things seem to work)
        if jj==length(internal_faces_and_cubes_index_array_no_double_counting)
            while counter<2
                if internal_faces_and_cubes_index_array_no_double_counting(jj,2)==internal_faces_and_cubes_index_array_no_double_counting(counter-1,3)
                    cube_copy(jj) = (NC+1);
                    disp("Case 4")
                end
                counter = counter - 1;
            end
        end
        jj = jj+1;
    end
    cubes_to_break = internal_faces_and_cubes_index_array(index_1,2:3);

    cubes_all = [internal_faces_and_cubes_index_array_no_double_counting cube_copy];

    % okay this seems to work now
    % so what is left to do is to store the indices of the NC+1 in one array
    % and the 0 in another one
    % those will make the two aggregates
    indices_NC = find(cubes_all(:,4)==(NC+1));
    cubes_all_xc_1 = cubes_all(indices_NC,2:3);
    cubes_all_xc_1 = cubes_all_xc_1(:);
    cubes_all_xc_1 = sort(cubes_all_xc_1);

    agg_idx_list = [];
    agg_idx_list(1) = cubes_all_xc_1(1);
    for kk=2:length(cubes_all_xc_1)
        if cubes_all_xc_1(kk)~=cubes_all_xc_1(kk-1)
            agg_idx_list = [agg_idx_list;cubes_all_xc_1(kk)];
        end
    end
    agg_idx_list;
    agg_idx_list = sort(agg_idx_list);
    agg_idx_list; 
    
    xc;
    xc_1 = xc(agg_idx_list,:);
    xc_2 = xc;
    xc_2(agg_idx_list,:) = [];
end

end
