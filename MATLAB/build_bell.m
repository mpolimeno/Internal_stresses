function [pos_bell_l,pos_bell_r] = build_bell(NC_bell,res)

% since the firt cube is at the origin, then I can just build one bell next
% to it and then copy each cube to the side of the dumb by shifting them by
% +10,+14,+18 in the y direction

% if NC_bell~=54
%     msg = 'Input must be 54';
%     error(msg)
% end

% build left bell
NC_bell_l = NC_bell/2;
%pos_bell_l = zeros(NC_bell_l+3,3);

% 3 cubes connecting dumb to bell
% for ii = 1:3
%     jj = ii*2;
%     pos_bell_l(ii,:) = [0,-jj,0];
% end

pos_bell_l = zeros(NC_bell_l,3);
% left-most x coordinate is y=-2, so for the left bell we will start at
% y=-4, and go from x=-2 to x=2
if res==1
    count = 1;
    for jj=-4:-2:-8
        for ii=-2:2:2
            pos_bell_l(count,:) = [ii,jj,0];
            count = count+1;
        end
    end
    % count=10; fill the next 9 elements in positive z direction
    count_zp = count;
    for kk = 1:count-1
        pos_bell_l(count_zp,:) = pos_bell_l(kk,:) + [0,0,2];
        count_zp = count_zp + 1;
    end

    %count_zp=19; fill the next 9 elements in negative z direction
    count_zm = count_zp;
    for mm = 1:count-1
        pos_bell_l(count_zm,:) = pos_bell_l(mm,:) + [0,0,-2];
        count_zm = count_zm + 1;
    end

end

if res==2
    
    count = 1;
    for jj=-10:-2:-14
        for ii=-2:2:2
            pos_bell_l(count,:) = [ii,jj,0];
            count = count+1;
        end
    end
    % count=10; fill the next 9 elements in positive z direction
    count_zp = count;
    for kk = 1:count-1
        pos_bell_l(count_zp,:) = pos_bell_l(kk,:) + [0,0,2];
        count_zp = count_zp + 1;
    end

    %count_zp=19; fill the next 9 elements in negative z direction
    count_zm = count_zp;
    for mm = 1:count-1
        pos_bell_l(count_zm,:) = pos_bell_l(mm,:) + [0,0,-2];
        count_zm = count_zm + 1;
    end
    
    count_xp = count_zm;
    for nn = 1:count-1
        pos_bell_l(count_xp,:) = pos_bell_l(nn,:) + [-6,0,0];
        count_xp = count_xp + 1;
    end

    count_xp_zp = count_xp;
    pos_bell_l(count_xp_zp:count_xp_zp+8,:) = pos_bell_l(count_xp_zp-9:count_xp_zp-1,:) + [0,0,2];

    count_xp_zm = count_xp_zp+9;
    pos_bell_l(count_xp_zm:count_xp_zm+8,:) = pos_bell_l(count_xp_zm-9:count_xp_zm-1,:) + [0,0,-4];
    
    % now I can copy this whole thing and go 6 units in -y
%     count_rest = length(pos_bell_l)/2;
%     pos_bell_l(count_rest+1:end,:) = pos_bell_l(1:count_rest,:) + [0,-6,0];
    
    % now I can copy this whole thing and go 6 units in -y
    count_rest_half = length(pos_bell_l)/4;
    pos_bell_l(count_rest_half+1:length(pos_bell_l)/2,:) = pos_bell_l(1:count_rest_half,:) + [0,-6,0];
    
    %copy in -z down 6
    count_rest = length(pos_bell_l)/2;
    pos_bell_l(count_rest+1:end,:) = pos_bell_l(1:count_rest,:) + [0,0,-6];
end


% build top 9 cubes by shifting pos_bell_l by +2 in the z-direction
% for mm=1:9
%     pos_bell_l(mm+ii,:) = pos_bell_l(mm,:) + [0,0,2];
% end
% 
% % build bottom 9 cubes by shifting pos_bell_l by -2 in the z-direction
% pos_bell_l(mm+ii+1:end,:) = pos_bell_l(ii+1:mm+ii,:);
% pos_bell_l(mm+ii+1:end,3) = pos_bell_l(mm+ii+1:end,3)*(-1);
% 
% % build right bell
pos_bell_r = pos_bell_l;
% 
% for nn=1:length(pos_bell_r)
%     if pos_bell_l(nn,2)==-2
%         pos_bell_r(nn,2) = pos_bell_l(nn,2) + 10;
%     elseif pos_bell_l(nn,2)==-4
%         pos_bell_r(nn,2) = pos_bell_l(nn,2) + 14;
%     else
%         pos_bell_r(nn,2) = pos_bell_l(nn,2) + 18;
%     end
% end


% now mirror the left bell to build the right bell;
% note that the rightmost y-coordinate is y = 4

if res==1
    for nn=1:length(pos_bell_r)
        if pos_bell_l(nn,2)==-4
            pos_bell_r(nn,2) = pos_bell_l(nn,2) + 10;
        elseif pos_bell_l(nn,2)==-6
            pos_bell_r(nn,2) = pos_bell_l(nn,2) + 14;
        elseif pos_bell_l(nn,2)==-8
            pos_bell_r(nn,2) = pos_bell_l(nn,2) + 18;
        end
    end
else
    for nn=1:length(pos_bell_r)
        if pos_bell_l(nn,2)==-10
            pos_bell_r(nn,2) = pos_bell_l(nn,2) + 18;
        elseif pos_bell_l(nn,2)==-12
            pos_bell_r(nn,2) = pos_bell_l(nn,2) + 22;
        elseif pos_bell_l(nn,2)==-14
            pos_bell_r(nn,2) = pos_bell_l(nn,2) + 26;
        elseif pos_bell_l(nn,2)==-16
            pos_bell_r(nn,2) = pos_bell_l(nn,2) + 30;
        elseif pos_bell_l(nn,2)==-18
            pos_bell_r(nn,2) = pos_bell_l(nn,2) + 34;
        elseif pos_bell_l(nn,2)==-20
            pos_bell_r(nn,2) = pos_bell_l(nn,2) + 38;
        end
    end
end
    

% these are just lines to be copied and pasted into the command window to
% verify that the dumbell is built correctly
% [finalposint, finalndir, finalori,Nf] = build_faces(pos_bell_l, 36);
% plot_faces(finalposint, finalndir, finalori,Nf,1,'c');
end