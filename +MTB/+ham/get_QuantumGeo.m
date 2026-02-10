function [Geo,Kx,Ky,berrycur,traceg,T,Sigma,K]=QuantumGeo(obj,Kx,Ky,Enk,Unk,bandindex,delta)
     %{2023-11-11. cal the quantum geometry tensor, whose real part give quantum metric
     %tensor and imag part give berry curvature
     %Kx,Ky: kmesh; obj: the continuum model objetc; 
     %Enk: knum*knum*Nbands
     %Unk: (eigenvector) Norbits*Nband*knum*knum
     %bandindex: eg: 1:3; delta: small term to avoid the gapless points
    % the output
    % Geo: the initial quantum geo tensor: knum*knum*nband*3, 3 for {xy,xx,yy}
    % berrycur: the intergral of berrycurvature on each small loop, so in
    % fact it should not be called berrycurvature here
    % traceg: the trace of real part of quantum geo tensor: knum*knum*nband
    % T,Sigma,K:voliation of trace condition, flucutation of
    % berrycurvature,quantum weigth of each band: 1*nband
%}       
            dimH=size(Unk,1);            
            knumx=size(Unk,3);
            knumy=size(Unk,4);
            nband=numel(bandindex);                                               
            dkx=[Kx(1,2)-Kx(1,1),0];dky=[0,Ky(2,1)-Ky(1,1)];                               
            dS=norm(dkx)*norm(dky);
            qm=@quan_metric;            
            dimG=[knumx,knumy,nband];
            itotal=prod(dimG);
            geoxx=zeros(1,itotal);
            geoyy=zeros(1,itotal);
            geoxy=zeros(1,itotal);
            parameter=cell(1,itotal);
            for l=1:itotal
                [i,j,band]=ind2sub(dimG,l);
                unk=Unk(:,:,i,j);
                enk=reshape(Enk(i,j,:),[dimH,1]);
                bandhere=bandindex(band);
                  k=[Kx(i,j),Ky(i,j)]; 
                parameter{l}={k,bandhere,unk,enk};
            end            
            parfor i=1:itotal
                k=parameter{i}{1};
                bandhere=parameter{i}{2};
                unk=parameter{i}{3};
                enk=parameter{i}{4};
                 geoxx(i)=qm(obj,k,dkx,dkx,unk,enk,bandhere,delta);
                 geoyy(i)=qm(obj,k,dky,dky,unk,enk,bandhere,delta);
                 geoxy(i)=qm(obj,k,dkx,dky,unk,enk,bandhere,delta);
            end
            Geo_xx=reshape(geoxx,dimG);
            Geo_yy=reshape(geoyy,dimG);
            Geo_xy=reshape(geoxy,dimG);
            Omega_k=zeros(knumx,knumy,nband,3);            
            Omega_k(:,:,:,1)=Geo_xy;
            Omega_k(:,:,:,2)=Geo_xx;
            Omega_k(:,:,:,3)=Geo_yy;                            
            Geo=Omega_k.*dS;
            berrycur=zeros(knumx,knumy,nband);
            traceg=zeros(knumx,knumy,nband);
            T=zeros(1,nband);
            Sigma=zeros(1,nband);
            K=zeros(1,nband);
            for i=1:nband
                 berrycur(:,:,i)=-2.*imag(Geo(:,:,i,1));
                 traceg(:,:,i)=real(Geo(:,:,i,2)+Geo(:,:,i,3));
                 T(i)=sum(traceg(:,:,i)-abs(berrycur(:,:,i)),'all');
                 Ci=round(sum(berrycur(:,:,i),'all')./2./pi);
                 Sigma(i)=sqrt(sum((berrycur(:,:,i)./2./pi.*knumx*knumy-Ci).^2,'all')./knumx./knumy);
                 K(i)=sum(abs(traceg),'all')./2./pi;
            end


                

end

 function Q = quan_metric(obj, k, dku, dkv, Uk, Ek, band_index, delta) %nodegerate form but with delta
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
            % Hku = (obj.geneHK(k + dku) - obj.geneHK(k)) ./ norm(dku);
            % Hkv = (obj.geneHK(k + dkv) - obj.geneHK(k)) ./ norm(dkv);
            Hku = (obj.get_hk([k + dku,0]) - obj.get_hk([k,0])) ./ norm(dku);
            Hkv = (obj.get_hk([k + dkv,0]) - obj.get_hk([k,0])) ./ norm(dkv);
            % Hku = (obj.HK(k + dku) - obj.HK(k)) ./ norm(dku);
            % Hkv = (obj.HK(k + dkv) - obj.HK(k)) ./ norm(dkv);
            Q = (u0' * Hku * Uk_re) * (Uk' * Hkv * u0);
 end
