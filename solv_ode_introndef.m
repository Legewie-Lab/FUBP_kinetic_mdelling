function Y1 = solv_ode_introndef(parameter)


    ve=28;%number of parameters
    x0= zeros(1,ve);%initial vector

    t=linspace(0,100000,100);%time course. very high to be in equilibrium
    [T,X]=ode15s(@(t,x) intron_definition_ode(t,x,parameter),t,x0);%solve the ode

    Y1=X(end,:);%get the solution for the equilibrium

    %calculate the different isoform ratios
    incl=Y1(:,17);
    skip=Y1(:,18);
    fullIR=[Y1(:,19)+sum(Y1(:,1:16),2)];
    fIR=Y1(:,20)+Y1(:,21)+Y1(:,22)+Y1(:,23);
    seIR=Y1(:,24)+Y1(:,25)+Y1(:,26)+Y1(:,27);
    su=[incl+skip+fullIR+fIR+seIR];
    Y1=[incl./su skip./su fullIR./su fIR./su seIR./su];% results with the 5 isoforms


end






