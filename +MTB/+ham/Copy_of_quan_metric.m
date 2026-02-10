 function Q = quan_metric(obj, k, dku, dkv, Uk, Ek, band_index, delta) %nodegerate form but with delta
            %Uk: eigenstate at k
            %Ek: eigenvalue at k
            %dku,dkv: min vector in K-mesh
            u0 = Uk(:, band_index); 
            E0 = Ek(band_index);
            %the left states
            Uk(:, band_index) = []; Ek(band_index) = [];
            dE = (E0(1) - Ek').^2;
            tem=find(dE<delta);
            dE(tem)=dE(tem)+delta;
            Uk_re = Uk ./ dE;
            Hku = (obj.get_hk(k + dku) - obj.get_hk(k)) ./ norm(dku);
            Hkv = (obj.get_hk(k + dkv) - obj.get_hk(k)) ./ norm(dkv);
            Q = (u0' * Hku * Uk_re) * (Uk' * Hkv * u0);
 end
