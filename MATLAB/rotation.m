function Q = rotation(n_dir,n_ori)

    i_hat = [1,0,0];
    j_hat = [0,1,0];
    k_hat = [0,0,1];
    
    if n_dir==1 && n_ori>0
        v3 = i_hat; % normal
        
        v1 = j_hat;
        v2 = k_hat;
    end
    if n_dir==1 && n_ori<0
        v3 = -i_hat; % normal
        
        v1 = -j_hat;
        v2 = k_hat;
    end
    if n_dir==2 && n_ori>0
        v3 = j_hat; % normal
        
        v1 = -i_hat;
        v2 = k_hat;
    end
    if n_dir==2 && n_ori<0
        v3 = -j_hat; % normal
        
        v1 = i_hat;
        v2 = k_hat;
    end
    if n_dir==3 && n_ori>0
        v3 = k_hat; % normal
        
        v1 = i_hat;
        v2 = j_hat;
    end
    if n_dir==3 && n_ori<0
        v3 = -k_hat; % normal
        
        v1 = -i_hat;
        v2 = j_hat;
    end
    
    Q = [v1;v2;v3];
end
