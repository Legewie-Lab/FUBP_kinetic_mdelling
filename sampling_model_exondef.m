% clear all;


initial_parameters=ones(1,15);

initial_parameters(4)=0.01;
initial_parameters(11:12)=0.01;

chan=7;
perturbation=[1.7,3.5];
p = 10000;   % Number of points
N = length(initial_parameters);   % Number of dimensions
lb = initial_parameters./chan;
ub = chan.*initial_parameters;
X = lhsdesign(p,N,'criterion','correlation');%latin hypercube sampling

for pert = perturbation

    random_param = bsxfun(@plus,lb,bsxfun(@times,X,(ub-lb)));%random parameters
     
     
    
    %solve and store the starting values      
    val=[];
    
    for j =1:size(random_param,1);
        val=[val;solv_ode_exondef(random_param(j,:))];
    end
        
       
    
    %calculation of PSI parfor the starting parameters
    skip = val(:,2);
    incl=val(:,1);
    PSI=incl./(skip+incl);
    
    
    
    
    %%% Model 1
    
    %affecting first long intron k2
    dperbk2=random_param;
    
    dperbk2(:,2)= random_param(:,2)./pert;%k2
    
    %affecting second long intron k3
    dperbk3=random_param;
    dperbk3(:,3)=random_param(:,3)./pert;%k3
    
    %affecting both introns k2&k3
    
    dperbk2k3=random_param;
    
    dperbk2k3(:,2:3)= random_param(:,2:3)./pert;%k2&k3
    
    %%Model2
    %affecting first long intron i1&i2
    
    dperbi1i12=random_param;
    
    dperbi1i12(:,9:10)= random_param(:,9:10)./pert;
    
    %affecting second long intron i2&i2

    dperbi2i12=random_param;
    
    
    dperbi2i12(:,8:9)= random_param(:,8:9)./pert;
    
    
    %affecting first and second long intron i1&i2&i2

    dperbi1i2i12=random_param;
    
    dperbi1i2i12(:,8:10)= random_param(:,8:10)./pert;
    
    
    
    
    %calculate the isoparforms parfor the different conditions
    valperbk2=[];
    valperbk3=[];
    valperbk2k3=[];
    
    valperbi1i12=[];
    valperbi2i12=[];
    valperbi1i2i12=[];
    
    for j =1:size(random_param,1);
    
        
        valperbk2=[valperbk2;solv_ode_exondef(dperbk2(j,:))];
        valperbk3=[valperbk3;solv_ode_exondef(dperbk3(j,:))];
        valperbk2k3=[valperbk2k3;solv_ode_exondef(dperbk2k3(j,:))];
    
        valperbi1i12=[valperbi1i12;solv_ode_exondef(dperbi1i12(j,:))];
        valperbi2i12=[valperbi2i12;solv_ode_exondef(dperbi2i12(j,:))];
        valperbi1i2i12=[valperbi1i2i12;solv_ode_exondef(dperbi1i2i12(j,:))];
    
    end
    
    
    
    
    %calculate PSI parfor the different conditions
    PSIperbk2=valperbk2(:,1)./(valperbk2(:,2)+valperbk2(:,1));
    PSIperbk3=valperbk3(:,1)./(valperbk3(:,2)+valperbk3(:,1));
    PSIperbk2k3=valperbk2k3(:,1)./(valperbk2k3(:,2)+valperbk2k3(:,1));
    
    PSIperbi1i12=valperbi1i12(:,1)./(valperbi1i12(:,2)+valperbi1i12(:,1));
    PSIperbi2i12=valperbi2i12(:,1)./(valperbi2i12(:,2)+valperbi2i12(:,1));
    PSIperbi1i2i12=valperbi1i2i12(:,1)./(valperbi1i2i12(:,2)+valperbi1i2i12(:,1));
    
    
    %calculate delta PSI
    dPSIk2=PSIperbk2-PSI;
    dPSIk3=PSIperbk3-PSI;
    dPSIk2k3=PSIperbk2k3-PSI;
    dPSIi1i12=PSIperbi1i12-PSI;
    dPSIi2i12=PSIperbi2i12-PSI;
    dPSIi1i2i12=PSIperbi1i2i12-PSI;
    
    
    
    
    
    %Plot the data parfor delta PSI
    figure()
    boxplot([dPSIk2,dPSIk3,dPSIk2k3],'Labels',["k2","k3","k2&k3"]);
    ylabel('\DeltaPSI')
    title(sprintf('FUBP perturbation binding of %G', pert));
    
    %plot the data parfor intron retention
    figure()
    boxplot([valperbk2(:,3:5)-val(:,3:5),valperbk3(:,3:5)-val(:,3:5),valperbk2k3(:,3:5)-val(:,3:5)],'Labels',["fullIrik2","firstIrk2","secondIrk2","fullIrik3","firstIrk3","secondIrk3","fullIrik2k3","firstIrk2k3","secondIrk2k3"]);
    ylabel('\intronretention')
    title(sprintf('FUBP perturbation binding of %G', pert));
    
    
    
    
    figure()
    boxplot([dPSIi1i12,dPSIi2i12,dPSIi1i2i12],'Labels',["i1&i12","i2&12","i1&i2&i12"]);
    ylabel('\DeltaPSI')
    title(sprintf('FUBP perturbation splicing of %G', pert));
    
    figure()
    boxplot([valperbi1i12(:,3:5)-val(:,3:5),valperbi2i12(:,3:5)-val(:,3:5),valperbi1i2i12(:,3:5)-val(:,3:5)],'Labels',["fullIri1&12","firstIri1&i12k3","secondIri1&i12","fullIri2&12","firstIri2&i12k3","secondIri2&i12","fullIri&i21&12","firstIri1&i2&i12k3","secondIri1&i2&i12"]);
    ylabel('\intronretention')
    title(sprintf('FUBP perturbation splicing of %G', pert));





% saving the data

    table_title=[{'k1','k2','k3','kret','k4','k5','k6','i1','i12','i2','kincl','kskip','kdr1','kdr2','s','inclusion','skipping','Full Intron retention','First intron retetion','Second Intron retention','PSI','Starting PSI','\DeltaPSI','initial inclusion','initial skipping',' initial Full Intron retention','initial First intron retetion','initial Second Intron retention'}];
    samplingk2=array2table([dperbk2,valperbk2,PSIperbk2,PSI,dPSIk2,val],'VariableNames',table_title);
    samplingk3=array2table([dperbk3,valperbk3,PSIperbk3,PSI,dPSIk3,val],'VariableNames',table_title);
    samplingk2k3=array2table([dperbk2k3,valperbk2k3,PSIperbk2k3,PSI,dPSIk2k3,val],'VariableNames',table_title);
    
    samplingi1i12=array2table([dperbi1i12,valperbi1i12,PSIperbi1i12,PSI,dPSIi1i12,val],'VariableNames',table_title);
    samplingi2i12=array2table([dperbi2i12,valperbi2i12,PSIperbi2i12,PSI,dPSIi2i12,val],'VariableNames',table_title);
    samplingi1i2i12=array2table([dperbi1i2i12,valperbi1i2i12,PSIperbi1i2i12,PSI,dPSIi1i2i12,val],'VariableNames',table_title);
    

    names = {'samplingk2', 'samplingk3', 'samplingk2&k3','samplingi1&i12','samplingi2&i12','samplingi1&i2&i12'};%name the sheets in the excel file
    
    % Create a new Excel file
    filename = ['exon_definition_data' num2str(pert) '.xlsx'];
    if exist(filename, 'file')
        delete(filename); % Delete the file if it already exists
    end
    
    % Write each table to a separate sheet in the Excel file

    writetable(samplingk2, filename, 'Sheet', names{1});
    writetable(samplingk3, filename, 'Sheet', names{2});
    writetable(samplingk2k3, filename, 'Sheet', names{3});
    writetable(samplingi1i12, filename, 'Sheet', names{4});
    writetable(samplingi2i12, filename, 'Sheet', names{5});
    writetable(samplingi1i2i12, filename, 'Sheet', names{6});

end





