function pos_dumb = build_dumb(NC_dumb,res)

% if NC_dumb~=4
%     msg = 'Input must be 4';
%     error(msg)
% end

pos_dumb = zeros(NC_dumb,3);

% bell is 4 cubes aligned with y-axis: go from y=-2 to y=4
count = 1;
for ii=-2:2:4
    pos_dumb(count,:) = [0,ii,0];
    count = count + 1;
end


% I want each cube face to become 4 cube faces

if res==2
    % increase resolution
%     count_top = 1;
%     for top=-2:2:4
%         pos_dumb(count,:) = pos_dumb(count_top,:) + [0,0,2];
%         count_top = count_top + 1;
%         count = count + 1;
%     end
    count = 1;
    for ii=-8:2:6
        pos_dumb(count,:) = [-2,ii,-2];
        count = count + 1;
    end
    count_bottom = 1;
    for bottom=-8:2:6
        pos_dumb(count,:) = pos_dumb(count_bottom,:) + [0,0,-2];
        count_bottom = count_bottom + 1;
        count = count + 1;
    end

    count_left = 1;
    for left=-8:2:6
        pos_dumb(count,:) = pos_dumb(count_left,:) + [-2,0,0];
        count_left = count_left + 1;
        count = count + 1;
    end

%     count_right = 1;
%     for right=-2:2:4
%         pos_dumb(count,:) = pos_dumb(count_right,:) + [2,0,0];
%         count_right = count_right + 1;
%         count = count + 1;
%     end

%     count_left_up = 1;
%     for left_up=-2:2:4
%         pos_dumb(count,:) = pos_dumb(count_left_up,:) + [-2,0,2];
%         count_left_up = count_left_up + 1;
%         count = count + 1;
%     end
    
%     count_right_up = 1;
%     for left=-2:2:4
%         pos_dumb(count,:) = pos_dumb(count_right_up,:) + [2,0,2];
%         count_right_up = count_right_up + 1;
%         count = count + 1;
%     end

    count_left_down = 1;
    for left_up=-8:2:6
        pos_dumb(count,:) = pos_dumb(count_left_down,:) + [-2,0,-2];
        count_left_down = count_left_down + 1;
        count = count + 1;
    end
% 
%     count_right_down = 1;
%     for left_up=-2:2:4
%         pos_dumb(count,:) = pos_dumb(count_right_down,:) + [2,0,-2];
%         count_right_down = count_right_down + 1;
%         count = count + 1;
%     end

end



% these are just lines to be copied and pasted into the command window to
% verify that the dumbell is built correctly
% [finalposint, finalndir, finalori,Nf] = build_faces(pos_dumb, 4);
% plot_faces(finalposint, finalndir, finalori,Nf,1,'c');
