 function Q = get_quan_metric_v2(obj, k, dku, dkv, Uk, Ek, delta) %nodegerate form but with delta
            %Uk: eigenstate at k Norbitals*Norbitals
            %Ek: eigenvalue at k Norbitals
            %dku,dkv: min vector in K-mesh
            u0=Uk;
            Ek=ones(1,size(Ek,1));
            dE = (Ek-Ek').^2;
            %the left states
            % dE = (E0(1) - Ek').^2;
            tem=find(dE<delta);
            dE(tem)=dE(tem)+delta;
            Uk_re = Uk ./ dE;
            Hku = (obj.get_hk(k + dku) - obj.get_hk(k)) ./ norm(dku);
            Hkv = (obj.get_hk(k + dkv) - obj.get_hk(k)) ./ norm(dkv);
            Q = (u0' * Hku * Uk_re) * (Uk' * Hkv * u0);
            
 end