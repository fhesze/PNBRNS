function spec = getChemotaxis(FlowData, GeometryData, TransportCoeffs, spec)

n           = GeometryData.GeometryCoeffs.xIncr;
m           = GeometryData.GeometryCoeffs.yIncr;
pore_len    = GeometryData.GeometryCoeffs.LengthOfPore;
node_no     = GeometryData.NodeData.NumberOfNodes;
nodes       = GeometryData.NodeData.NodesOfPores;
pore_pos    = GeometryData.PoreData.PorePosition;

% flux        = FlowData.value;
radi        = FlowData.radii;
mean_radi   = mean(radi);

pores_no    = n/2*(m-1) + (n/2-1)*fix(m/2);                                % number of pores
pore_CS     = pi.*radi.^2;                                                 % pore cross section
% pore_vol    = pi.*radi.^2.*pore_len;                                       % volume of pores
vol         = pore_CS.*pore_len;                                           % volume of pores in network

dt          = TransportCoeffs.dt;

%%-------------------------------------------------------------------------

conc1       = spec.chem;
bac1        = spec.bac;

%%-- bacterial and chemotaxis parameters ----------------------------------

i_case      = 6;
chi_sub     = [3e-7 6e-7 9e-7 3e-6 6e-6 9e-6];
chi_bac     = mean_radi^2*pi*pore_len.*[6e-13 6e-13 6e-13 6e-13 6e-13 6e-13];    % volumetric chi
file_order  = {'higherXs1' 'higherXs2' 'higherXs3' 'higherXs4' 'higherXs5' 'higherXs6'};

bac1_dens   = bac1./vol;
max_bac_dens=10000/(mean_radi^2*pi*pore_len);                              % means 10000 cells can be in a pore with avg radius = 1.2434e8 cells/ml
diffb1_coff = 5.19e-10;
chib1c1     = chi_sub(i_case);
chib1b1     = chi_bac(i_case);

fact1       = diffb1_coff*dt/pore_len^2;
fact2       = chib1c1*dt/pore_len;
fact3       = chib1b1*dt/pore_len;

% bac_stay    = zeros(pores_no,1);
% bac_movin   = zeros(pores_no,1);
bac_movin_sum=zeros(pores_no,1);
prob        = zeros(1, pores_no);

    for pore_i=1:pores_no
    
        node_fard=pore_pos(pore_i,1);
        node_zoj=pore_pos(pore_i,2);
    
        pores_fard_tot=nodes(node_fard, [1,2,3]);
        pores_fard_tot(pores_fard_tot==0)=[];
    
        pores_fard=pores_fard_tot;
        pores_fard(pores_fard==pore_i)=[];
    
    
        pores_zoj_tot=nodes(node_zoj, [1,2,3]);
        pores_zoj_tot(pores_zoj_tot==0)=[];
    
        pores_zoj=pores_zoj_tot;
        pores_zoj(pores_zoj==pore_i)=[];
       
        %finding the values for attractant gradients, deltak and delta^2k
        if nodes(node_zoj,7)~=0
            %BC conditions here
           if nodes(node_zoj,7)==1
               %concBC=conc1_bc1;%BC input
               concBC=conc1(pore_i);
               bac1BC=bac1_dens(pore_i);
           elseif nodes(node_zoj,7)==2
               concBC=conc1(pore_i);
               bac1BC=bac1_dens(pore_i);
           end
            drive1=(concBC-sum(conc1(pores_fard_tot).*pore_CS(pores_fard_tot))/sum(pore_CS(pores_fard_tot)))...
               /pore_len;
          driveb1=(bac1BC-sum(bac1_dens(pores_fard_tot).*pore_CS(pores_fard_tot))/sum(pore_CS(pores_fard_tot)))...
               /pore_len;
%             if isempty(pores_zoj)
%                 concBC2=concBC;
%             else
%                concBC2=(conc1(pores_zoj)+concBC)/2;%only when there are 3plix connections
%            end
           drive2=(concBC-conc1(pore_i)...
               -(conc1(pore_i)-sum(conc1(pores_fard).*pore_CS(pores_fard))/sum(pore_CS(pores_fard))))...
               /pore_len^2;
            driveb2=(bac1BC-bac1_dens(pore_i)...
                -(bac1_dens(pore_i)-sum(bac1_dens(pores_fard).*pore_CS(pores_fard))/sum(pore_CS(pores_fard))))...
                /pore_len^2;
        elseif nodes(node_fard,7)~=0
           %BC conditions here
           if nodes(node_fard,7)==1
               %concBC=conc1_bc1;%BC input
               concBC=conc1(pore_i);
               bac1BC=bac1_dens(pore_i);
           elseif nodes(node_fard,7)==2
               concBC=conc1(pore_i);
               bac1BC=bac1_dens(pore_i);
           end
           drive1=(sum(conc1(pores_zoj_tot).*pore_CS(pores_zoj_tot))/sum(pore_CS(pores_zoj_tot))-concBC)...
               /pore_len;
            driveb1=(sum(bac1_dens(pores_zoj_tot).*pore_CS(pores_zoj_tot))/sum(pore_CS(pores_zoj_tot))-bac1BC)...
               /pore_len;
%            if isempty(pores_fard)
%                 concBC2=concBC;
%            else
%                concBC2=(conc1(pores_fard)+concBC)/2;%only when there are 3plix connections
%            end
            drive2=((sum(conc1(pores_zoj).*pore_CS(pores_zoj))/sum(pore_CS(pores_zoj))-conc1(pore_i))...
                -(conc1(pore_i)-concBC))...
                /pore_len^2;
            driveb2=((sum(bac1_dens(pores_zoj).*pore_CS(pores_zoj))/sum(pore_CS(pores_zoj))-bac1_dens(pore_i))...
                -(bac1_dens(pore_i)-bac1BC))...
                /pore_len^2;
        else
          drive1=(sum(conc1(pores_zoj_tot).*pore_CS(pores_zoj_tot))/sum(pore_CS(pores_zoj_tot))...
               -sum(conc1(pores_fard_tot).*pore_CS(pores_fard_tot))/sum(pore_CS(pores_fard_tot)))...
            /pore_len;
            driveb1=(sum(bac1_dens(pores_zoj_tot).*pore_CS(pores_zoj_tot))/sum(pore_CS(pores_zoj_tot))...
            -sum(bac1_dens(pores_fard_tot).*pore_CS(pores_fard_tot))/sum(pore_CS(pores_fard_tot)))...
             /pore_len;
         drive2=((sum(conc1(pores_zoj).*pore_CS(pores_zoj))/sum(pore_CS(pores_zoj))-conc1(pore_i))...
             -(conc1(pore_i)-sum(conc1(pores_fard).*pore_CS(pores_fard))/sum(pore_CS(pores_fard))))...
             /pore_len^2;
         driveb2=((sum(bac1_dens(pores_zoj).*pore_CS(pores_zoj))/sum(pore_CS(pores_zoj))-bac1_dens(pore_i))...
             -(bac1_dens(pore_i)-sum(bac1_dens(pores_fard).*pore_CS(pores_fard))/sum(pore_CS(pores_fard))))...
             /pore_len^2;
        end
   
        if isempty(pores_zoj) || isempty(pores_fard)
          prob(pore_i)=1-fact1+fact2*(drive2*pore_len+drive1*(pore_CS(pore_i)/sum(pore_CS(pores_zoj_tot))...
          -pore_CS(pore_i)/sum(pore_CS(pores_fard_tot))))...
            +fact3*(driveb2*pore_len+driveb1*(pore_CS(pore_i)/sum(pore_CS(pores_zoj_tot))...
         -pore_CS(pore_i)/sum(pore_CS(pores_fard_tot))));
        else
            prob(pore_i)=1-2*fact1+fact2*(drive2*pore_len+drive1*(pore_CS(pore_i)/sum(pore_CS(pores_zoj_tot))...
         -pore_CS(pore_i)/sum(pore_CS(pores_fard_tot))))...
         +fact3*(driveb2*pore_len+driveb1*(pore_CS(pore_i)/sum(pore_CS(pores_zoj_tot))...
         -pore_CS(pore_i)/sum(pore_CS(pores_fard_tot))));
        end

        if ~isempty(pores_zoj)
            prob(pores_zoj)=fact1.*pore_CS(pores_zoj)./sum(pore_CS(pores_zoj))...
                +fact2*drive1.*pore_CS(pores_zoj)./sum(pore_CS(pores_zoj_tot))...
                +fact3*driveb1.*pore_CS(pores_zoj)./sum(pore_CS(pores_zoj_tot));
        end
        if ~isempty(pores_fard)
            prob(pores_fard)=fact1.*pore_CS(pores_fard)./sum(pore_CS(pores_fard))...
                -fact2*drive1.*pore_CS(pores_fard)./sum(pore_CS(pores_fard_tot))...
                -fact3*driveb1.*pore_CS(pores_fard)./sum(pore_CS(pores_fard_tot));
        end
    
       bacno=bac1(pore_i);
        if isempty(pores_zoj)
            bac_mov_no = bac_prob2no( bacno,[prob(pore_i) prob(pores_fard)] );
            bac_stay(pore_i)= bac_mov_no(1);
             bac_movin(pores_fard)= bac_mov_no(2:2+length(pores_fard)-1);
           bac_movin_sum(pores_fard)=bac_movin_sum(pores_fard) + bac_movin(pores_fard)';
        elseif isempty(pores_fard)
            bac_mov_no = bac_prob2no( bacno,[prob(pore_i) prob(pores_zoj)] );
            bac_stay(pore_i)= bac_mov_no(1);
            bac_movin(pores_zoj)= bac_mov_no(2:2+length(pores_zoj)-1);
%             size(bac_movin_sum(pores_zoj))
%             size(bac_movin(pores_zoj))
            bac_movin_sum(pores_zoj)=bac_movin_sum(pores_zoj) + bac_movin(pores_zoj)';
        else
            bac_mov_no = bac_prob2no( bacno,[prob(pore_i) prob(pores_zoj) prob(pores_fard)] );
            bac_stay(pore_i)= bac_mov_no(1);
            bac_movin(pores_zoj)= bac_mov_no(2:2+length(pores_zoj)-1);
            bac_movin(pores_fard)= bac_mov_no(2+length(pores_zoj):2+length(pores_zoj)+length(pores_fard)-1);
            bac_movin_sum([pores_fard pores_zoj])=bac_movin_sum([pores_fard pores_zoj]) + bac_movin([pores_fard pores_zoj])';
        end

    
    end

bac1        = bac_stay' + bac_movin_sum;
bac1_dens   = bac1./vol;

%------------ Over flow of bactria in a pore ------------------------------
    
    while bac1_dens(bac1_dens>max_bac_dens)
        
        bac_movin_sum=zeros(pores_no,1);
        ind_bac_over=find(bac1_dens>max_bac_dens);
        
        for pore_i= ind_bac_over'
            node_fard=pore_pos(pore_i,1);
            node_zoj=pore_pos(pore_i,2);
            
            pores_fard_tot=nodes(node_fard, [1,2,3]);
            pores_fard_tot(pores_fard_tot==0)=[];
            
            pores_fard=pores_fard_tot;
            pores_fard(pores_fard==pore_i)=[];
            
            
            pores_zoj_tot=nodes(node_zoj, [1,2,3]);
            pores_zoj_tot(pores_zoj_tot==0)=[];
            
            pores_zoj=pores_zoj_tot;
            pores_zoj(pores_zoj==pore_i)=[];
            
            if ~isempty(pores_zoj)
                if bac1(pores_zoj)
                    prob(pores_zoj)=1./bac1(pores_zoj);
                else
                    prob(pores_zoj)=.5;
                end
            end
            if ~isempty(pores_fard)
                if bac1(pores_fard)
                    prob(pores_fard)=1./bac1(pores_fard);
                else
                    prob(pores_fard)=.5;
                end
            end
            
            
            %prob(pore_i)=1/bac1(pore_i);
            bacno=bac1(pore_i)-fix(max_bac_dens*vol(pore_i));
            bac1(pore_i)=bac1(pore_i)-bacno;
            if isempty(pores_zoj)
                bac_mov_no = bac_prob2no( bacno, prob(pores_fard) );
                %bac_stay(pore_i)= bac_mov_no(1);
                bac_movin(pores_fard)= bac_mov_no(1:length(pores_fard));
                bac_movin_sum(pores_fard)=bac_movin_sum(pores_fard)+bac_movin(pores_fard);
            elseif isempty(pores_fard)
                bac_mov_no = bac_prob2no( bacno, prob(pores_zoj) );
                %bac_stay(pore_i)= bac_mov_no(1);
                bac_movin(pores_zoj)= bac_mov_no(1:length(pores_zoj));
                bac_movin_sum(pores_zoj)=bac_movin_sum(pores_zoj)+bac_movin(pores_zoj);
            else
                bac_mov_no = bac_prob2no( bacno,[prob(pores_zoj) prob(pores_fard)] );
                %bac_stay(pore_i)= bac_mov_no(1);
                bac_movin(pores_zoj)= bac_mov_no(1:length(pores_zoj));
                bac_movin(pores_fard)= bac_mov_no(1+length(pores_zoj):length(pores_zoj)+length(pores_fard));
                bac_movin_sum([pores_fard pores_zoj])=bac_movin_sum([pores_fard pores_zoj])+bac_movin([pores_fard pores_zoj]);
            end
        end
        bac1=bac1+bac_movin_sum;
        bac1_dens=bac1./vol;
        %sum(bac1)
    end
    
%-------------- /end Over flow of bactria in a pore -----------------------

spec = struct('chem', conc1, 'bac', bac1);

end
