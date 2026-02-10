function Q = get_slab_quan_metric_v2(obj,del_index, k, dku, dkv, nslab, Uk, Ek, band_index, delta) %nodegerate form but with delta
            %Uk: eigenstate at k
            %Ek: eigenvalue at k
            %dku,dkv: min vector in K-mesh
            u0 = Uk(:, band_index); E0 = Ek(band_index);
            %the left states
            Uk(:, band_index) = []; Ek(band_index) = [];
            dE = (E0(1) - Ek').^2;
            tem=find(dE<delta);
            dE(tem)=dE(tem)+delta;
            Uk_re = Uk ./ dE;
            nbands=size(obj.ham,1);
            nrpts=size(obj.ham,3);
            Hku = (MTB.ham.get_slab_hk_v2(obj.ham,obj.hopr2,del_index,nslab,nbands,nrpts,k+dku,obj.a2) - MTB.ham.get_slab_hk_v2(obj.ham,obj.hopr2,del_index,nslab,nbands,nrpts,k,obj.a2)) ./ norm(dku);
            Hkv = (MTB.ham.get_slab_hk_v2(obj.ham,obj.hopr2,del_index,nslab,nbands,nrpts,k+dkv,obj.a2) - MTB.ham.get_slab_hk_v2(obj.ham,obj.hopr2,del_index,nslab,nbands,nrpts,k,obj.a2)) ./ norm(dkv);
            Q = (u0' * Hku * Uk_re) * (Uk' * Hkv * u0);
 end
