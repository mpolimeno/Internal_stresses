function DL = DoubleLayer_DL(ai,bj,myposint,myndir,myori,pos0,M,flow)

    % k_hat
    k_hat = [0;0;1];

%     I_list = [];
    DL = [0;0;0];
    
    xc = myposint; % center of each external face
    n_dir = myndir; % 1,2,3 for i,j,k
    n_ori = myori; % 1 or -1
    x0 = pos0;

    % get Q for ff face
    Q = rotation(n_dir,n_ori);
    Qt = Q';

    n_hat = Qt*k_hat; % normal vector for the mapped system
    %integrand
    Intg_ij = [0,0,0];
    for ii=1:length(ai)
        for jj=1:length(bj)

            % build gamma
            gammaij = [ai(ii,jj);bj(ii,jj);0];

            % shift to origin
            shift = xc;

            % build Tensor T
            Tl = (Qt*gammaij+shift')-x0';
            Ti = Tl';
            Ts = (shift'-x0');

            % background flow:
            if flow==1
                Ubg = ([0;0;1]);
            elseif flow==2
                Ubg = (cross([0;0;1]',(Qt*gammaij+shift')))';
            elseif flow==3
                Ubg = ([0;0;1]' + cross([0;0;1]',(Qt*gammaij+shift')))';
            elseif flow==4
                x_ext = Qt*gammaij+shift';
                Ubg = M*x_ext;
            else
                x_pos = (Qt*gammaij+shift');
                if x_pos(:,2)>=0
                    Ubg = (cross([0;0;1]',(Qt*gammaij+shift')))';
                else
                    Ubg = (cross([0;0;-1]',(Qt*gammaij+shift')))';
                end
            end
            % for extensional flow
            %M = [0,0,0;0,1,0;0,0,-1];
            %x_ext = Qt*gammaij+shift';
            %Ubg = M*x_ext;
            
            % numerator
            num = (Ubg'*Tl)*Ti*(Ts'*n_hat);

            % denominator
            den = norm((Qt*gammaij+shift')-x0')^5;

            % evaluate integrand
            if norm((Qt*gammaij+shift')-x0')<1e-3
                disp("BOOM") % we blow up
                Integrand = [0,0,0]; % at the singularity
            else
                Integrand = num/den; % elsewhere
            end

            % Riemann sum
            Intg_ij = Intg_ij + Integrand;
        end
    end
    Integral = Intg_ij;
%     I_list = [I_list Integral']; % to list faces
    %DL = DL + I_list(:,ff);
    DL = DL + Integral';
end