function dxdt = exon_definition_ode(t,x,param)


    


    %different parameters
    k1=param(1);
    k2=param(2);
    k3=param(3);
    kret=param(4);
    k4=param(5);
    k5=param(6);
    k6=param(7);
    i1=param(8);
    i12=param(9);
    i2=param(10);
    kincl=param(11);
    kskip=param(12);
    kdr1=param(13);
    kdr2=param(14);
    s=param(15);
    
    
    %different variables
    P0_0_0 = x(1);
    P1_0_0= x(2);
    P0_1_0=x(3);
    P0_0_1=x(4);
    P1_1_0=x(5);
    P1_0_1=x(6);
    P0_1_1=x(7);
    P1_1_1=x(8);
    incl=x(9);
    skip=x(10);
    ret=x(11);
    P0_11=x(12);
    P1_11=x(13);
    P11_0=x(14);
    P11_1=x(15);
    deg=x(16);


    %Ordinary differential equations
    dP0_0_0dt = (k4 * P1_0_0) + (k5 * P0_1_0) + (k6 * P0_0_1) +s - (k1 + k2 + k3 + kret) *P0_0_0;
    dP1_0_0dt = k1 * P0_0_0 + k5 * P1_1_0 + k6 * P1_0_1 - (k2 + k3 + k4 + kret) * P1_0_0;
    dP0_1_0dt = k2 * P0_0_0 + k4 * P1_1_0 + k6 * P0_1_1 - (k1 + k3 + k5 + kret) * P0_1_0;
    dP0_0_1dt = k3 * P0_0_0 + k4 * P1_0_1 + k5 * P0_1_1 - (k1 + k2 + k6 + kret) * P0_0_1;
    dP1_1_0dt = k1 * P0_1_0 + k2 * P1_0_0 + k6 * P1_1_1 - (k3 + k4 + k5 + kret + i1) * P1_1_0;
    dP1_0_1dt = k1 * P0_0_1 + k3 * P1_0_0 + k5 * P1_1_1 - (k2 + k4 + k6 + kret + i12) * P1_0_1;
    dP0_1_1dt = k2 * P0_0_1 + k3 * P0_1_0 + k4 * P1_1_1 - (k1 + k6 + k5 + kret + i2) * P0_1_1;
    dP1_1_1dt = k1 * P0_1_1 + k2 * P1_0_1 + k3 * P1_1_0 - (k4 + k5 + k6 + kret + i1 + i2) * P1_1_1;
    dincldt = (P1_11*i1)+(P11_1*i2)-(incl*kincl);
    dskipdt =i12*P1_0_1-(skip*kskip);
    dretdt =kret*(P0_0_0+P1_0_0+P0_1_0+P0_0_1+P1_1_0+P1_0_1+ P0_1_1 +P1_1_1)-(ret*(kdr1+kdr2));
    dP0_11dt =i2*P0_1_1+k6*P1_11-(k3+kdr1)*P0_11;
    dP1_11dt =k3*P0_11+i2*P1_1_1-(k6+i1+kdr1)*P1_11;
    dP11_0dt =i1*P1_1_0+k4*P11_1-(k1+kdr2)*P11_0;
    dP11_1dt =k1*P11_0+i1*P1_1_1-(k4+i2+kdr1)*P11_1;
    deg = kincl*incl+kskip*skip+kdr1*(P0_11+P1_11)+kdr2*(P11_0+P11_1)+((kdr1+kdr2)*(ret));

    dxdt = [dP0_0_0dt;dP1_0_0dt;dP0_1_0dt;dP0_0_1dt;dP1_1_0dt;dP1_0_1dt;dP0_1_1dt;dP1_1_1dt;dincldt;dskipdt;dretdt;dP0_11dt;dP1_11dt;dP11_0dt;dP11_1dt;deg];

end


