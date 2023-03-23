function Y = solv_ode_exondef(parameter)


    ve=16;%number of parameters
    x0= zeros(1,ve);%initial vector
    
    t=linspace(0,100000,100);%time course. very high to be in equilibrium
    [T,X]=ode15s(@(t,x) exon_definition_ode(t,x,parameter),t,x0);%solve the ode
    
    Y1=X(end,:);%get the solution for the equilibrium
    
    %calculate the different isoform ratios

    incl=Y1(:,9);
    skip=Y1(:,10);
    fullIR=[Y1(:,11)+sum(Y1(:,1:8),2)];
    fIR=Y1(:,12)+Y1(:,13);
    seIR=Y1(:,14)+Y1(:,15);
    su=[incl+skip+fullIR+fIR+seIR];
    Y=[incl./su skip./su fullIR./su fIR./su seIR./su];% results with the 5 isoforms


end






