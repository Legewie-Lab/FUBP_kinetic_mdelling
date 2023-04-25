function dxdt = intron_definition_ode(t,x,param)


    


    %different parameters
    k1=param(1);
    k2a=param(2);
    k2b=param(3);
    k3=param(4);
    kret=param(5);
    k4=param(6);
    k5a=param(7);
    k5b=param(8);
    k6=param(9);
    i1=param(10);
    i12=param(11);
    i2=param(12);
    kincl=param(13);
    kskip=param(14);
    kdr1=param(15);
    kdr2=param(16);
    s=param(17);
    
    

    %different variables
    P0_0_0=x(1);
    P1_0_0=x(2);
    P0_a_0=x(3);
    P0_b_0=x(4);
    P0_1_0=x(5);
    P0_0_1=x(6);
    P1_a_0=x(7);
    P1_b_0=x(8);
    P1_1_0=x(9);
    P1_0_1=x(10);
    P0_a_1=x(11);
    P0_b_1=x(12);
    P0_1_1=x(13);
    P1_a_1=x(14);
    P1_b_1=x(15);
    P1_1_1=x(16);
    incl=x(17);
    skip=x(18);
    ret=x(19);
    P0_b1=x(20);
    P1_b1=x(21);
    P0_11=x(22);
    P1_11=x(23);
    P1a_0=x(24);
    P1a_1=x(25);
    P11_0=x(26);
    P11_1=x(27);
    deg=x(28);
    
    
    %Ordinary differential equations
    dP0_0_0dt = k4*P1_0_0 +k5a*P0_a_0 + k5b*P0_b_0 +k6* P0_0_1 + s -(k1+k2a+k2b+k3+kret) * P0_0_0;
    dP1_0_0dt = k1*P0_0_0 + k5a*P1_a_0 + k5b*P1_b_0 + k6*P1_0_1 - (k2a + k2b + k3 +k4+kret)*P1_0_0 ;
    dP0_a_0dt = k2a*P0_0_0 + k4*P1_a_0 + k5b*P0_1_0 + k6*P0_a_1 - (k1 + k2b + k3 + k5a + kret)*P0_a_0;
    dP0_b_0dt = k2b*P0_0_0 + k4*P1_b_0 + k5a*P0_1_0 + k6*P0_b_1 - (k1+k2a+k3+k5b+kret)*P0_b_0;    
    dP0_1_0dt = k2b*P0_a_0 + k2a*P0_b_0 + k4*P1_1_0 + k6*P0_1_1 - (k1+k3+k5a+k5b+kret)*P0_1_0;
    dP0_0_1dt = k3*P0_0_0 + k4*P1_0_1  + k5a*P0_a_1 + k5b*P0_b_1- (k1+k2a+k2b+k6+kret)*P0_0_1;
    dP1_a_0dt = k1*P0_a_0 + k2a*P1_0_0 + k5b*P1_1_0 + k6*P1_a_1 - (k2b+k3+k4+k5a+kret+i1)*P1_a_0;
    dP1_b_0dt = k1*P0_b_0 + k2b*P1_0_0 + k5a*P1_1_0 + k6*P1_b_1 - (k2a+k3+k4+k5b+kret)*P1_b_0;
    dP1_1_0dt = k1*P0_1_0 + k2a*P1_b_0 + k2b*P1_a_0 + k6*P1_1_1 - (k3+k4+k5a+k5b+kret+i1)*P1_1_0;
    dP1_0_1dt = k1*P0_0_1 + k3*P1_0_0 + k5a*P1_a_1 + k5b*P1_b_1 - (k2a + k2b + k4 + k6 + kret + i12)*P1_0_1;
    dP0_a_1dt = k2a*P0_0_1 + k3*P0_a_0 + k4*P1_a_1 + k5b*P0_1_1 - (k1 + k2b + k5a +k6 + kret)*P0_a_1;
    dP0_b_1dt = k2b*P0_0_1 + k3*P0_b_0 + k4*P1_b_1 + k5a*P0_1_1 - (k1 + k2a + k5b +k6 + kret+i2)*P0_b_1;
    dP0_1_1dt = k2a*P0_b_1 + k2b*P0_a_1 + k3*P0_1_0 + k4*P1_1_1 - (k1 +k5a + k5b + k6 + kret + i2 )*P0_1_1;
    dP1_a_1dt = k1*P0_a_1 + k2a*P1_0_1 + k3*P1_a_0 + k5b*P1_1_1 - (k2b+k4+k5a+k6+kret+i1)*P1_a_1;
    dP1_b_1dt = k1*P0_b_1 + k2b*P1_0_1 + k3*P1_b_0 + k5a*P1_1_1 - (k2a+k4+k5b+k6+kret+i2)*P1_b_1;
    dP1_1_1dt = k1*P0_1_1 + k2a*P1_b_1  + k2b*P1_a_1+ k3*P1_1_0 - (k4+k5a+k5b+k6+kret+i1+i2)*P1_1_1;
    dincldt = i1*P1_11 + i2*P11_1 - kincl*incl;
    dskipdt = i12 * P1_0_1  - kskip * skip;
    dretdt = kret*(P0_0_0 + P1_0_0 + P0_a_0 + P0_b_0 + P0_1_0 + P0_0_1 + P1_1_0 + P1_a_0 + P1_b_0 + P1_0_1 +  P0_a_1 + P0_b_1 + P0_1_1 + P1_a_1 + P1_b_1 + P1_1_1) - (kdr1 + kdr2)*ret;
    dP0_b1dt = i2 * P0_b_1 + k4 * P1_b1 + k5a*P0_11 - (k1 + k2a + kdr1)*P0_b1;
    dP1_b1dt = i2 * P1_b_1 + k1 * P0_b1 + k5a * P1_11 - (k2a + k4 + kdr1 ) * P1_b1;
    dP0_11dt = i2 * P0_1_1 + k2a * P0_b1 + k4 * P1_11 - (k1 + k5a +kdr1 ) * P0_11;
    dP1_11dt = i2 * P1_1_1 + k1 * P0_11 + k2a * P1_b1 - (k4 + k5a + i1 + kdr1) * P1_11;
    dP1a_0dt = i1 * P1_a_0 + k5b * P11_0+ k6 * P1a_1 - (k2b +k3 + kdr2) * P1a_0;
    dP1a_1dt = i1 * P1_a_1 + k3 * P1a_0 + k5b * P11_1 - (k2b + k6 + kdr2) * P1a_1;
    dP11_0dt = i1 * P1_1_0 + k2b * P1a_0  + k6 * P11_1 - (k3 + k5b + kdr2) * P11_0;
    dP11_1dt = i1 * P1_1_1 + k2b * P1a_1 + k3 * P11_0 - (k5b + k6 + i2 + kdr2)*P11_1;
    ddegdt = kincl * incl + kskip * skip + kdr1*(P0_b1 + P1_b1 + P0_11 + P1_11) + kdr2*(P1a_0+P1a_1+P11_0+P11_1) + (kdr1+kdr2)*ret;
    dxdt=[dP0_0_0dt;dP1_0_0dt;dP0_a_0dt;dP0_b_0dt;dP0_1_0dt;dP0_0_1dt;dP1_a_0dt;dP1_b_0dt;dP1_1_0dt;dP1_0_1dt;dP0_a_1dt;dP0_b_1dt;dP0_1_1dt;dP1_a_1dt;dP1_b_1dt;dP1_1_1dt;dincldt;dskipdt;dretdt;dP0_b1dt;dP1_b1dt;dP0_11dt;dP1_11dt;dP1a_0dt;dP1a_1dt;dP11_0dt;dP11_1dt;ddegdt];                                                                     
        
end


