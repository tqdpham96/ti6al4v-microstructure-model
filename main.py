import math 	
liquidus = 1650 

# Displacive transformation
def transformation(Tet,Te1t,tt,t1t,fnt,b_ret,a_p0,bet_0,beta_eq, Ms, gamma): 
    if fnt[0] > b_ret: 
        if Te1t < Ms: 
            fnt[2]= a_p0+(bet_0-b_ret)*(1- math.exp(-gamma*(Ms-Te1t)))  
            fnt[0]= 1- fnt[2] - fnt[1] 
    return fnt 
 
# Diffusional transformation
def difussion(Ted,Te1d,td,t1d,k0d,n0d,fb_equi,fnd): 
    if fnd[0]>fb_equi: 
        tfd=((-math.log((fnd[0]-fb_equi)/(1-fb_equi)))/k0d)**(1/n0d) 
        fnd[0]=1-fnd[2]-(1-math.exp(-k0d*((tfd+t1d-td)**n0d)))*(1-fnd[2]-fb_equi)  
        fnd[1]=1-fnd[0]-fnd[2] 
    return fnd 

# Equilibrium state
def equil(T, fa_eq, fb_eq, fa_pr_eq, phase3, phase4, T_beta_tran): 
    if T < T_beta_tran: 
        fa_eq = 0.925*(1 - math.exp(-0.0085*(T_beta_tran - T)))  
    else: 
        fa_eq = 0 
    # beta EVF in the mixture that just contains alpha-beta phase 
    fb_eq0= 1 - fa_eq 
    fb_eq=fb_eq0*(1 - phase3 - phase4) 
 
    # alpha prime EVF 
    fa_pr_eq=(1 + math.tanh((450-T)/80))/2 
     
    return fa_eq, fb_eq, fa_pr_eq, phase3, phase4 
 
# With alpha_prime
def x_alpha_prime(tx,t1x,k0x,n0x,fa_p_equix,fnx): 
    if fnx[2] > fa_p_equix: 
        tfx=((-math.log((fnx[2]-fa_p_equix)/(1-fa_p_equix)))/k0x)**(1/n0x) 
        fa_p0=fnx[2] 
        fnx[2]=1-(1-math.exp(-k0x*((tfx+t1x-tx)**n0x)))*(1-fa_p_equix)  
        fnx[1]=fnx[1]*((1-fnx[2])/(1-fa_p0)) 
        fnx[0]=fnx[0]*((1-fnx[2])/(1-fa_p0)) 
    else: 
        fnx = fnx 
    return fnx 
 
# Without alpha_prime
def no_x_alpha_prime(tn,t1n,k0n,n0n,fb_equi,fa_equi,fnn): 
    if fnn[1] > fa_equi: 
        tf=((-math.log((fnn[1]-fa_equi)/(1-fa_equi)))/k0n)**(1/n0n) 
        fnn[0]=(1-math.exp(-k0n*((tf+t1n-tn)**n0n)))*fb_equi #!Beta value 
        fnn[1]=1-fnn[0]-fnn[2] #!Alpha value          
    else: 
        fnn = fnn 
    return fnn 
         
# Utils function to calculate TDot
def cal_TDot(T1, T2, t1, t2): 
    return (T2-T1)/(t2-t1) 
 
# Linear Interpolation
def interpolate(xData, yData, xVal): 
    for i in range(len(xData)): 
        if xVal < xData[0]: 
            yVal = yData[0] 
        elif xVal > xData[-1]: 
            yVal = yData[-1] 
        elif xData[i] <= xVal <= xData[i+1]: 
            weight = (xVal - xData[i])/(xData[i+1] - xData[i]) 
            yVal = (1-weight)*yData[i] + weight*yData[i+1] 
    return yVal 
 
# Get final microstructure
def calculateFinal(TData, tData, TKN_ref1, TKN_ref2, T_beta_tran, Ms, Bs, gamma): 
    alpha_equili, beta_equili, alpha_prime_equili = 0, 0, 0 
    phase_n = [1, 0, 0, 0] #beta, alpha, alpha_prime 
    beta_phase = [] 
    alpha_phase = [] 
    alpha_prime_phase = [] 
    hb = 300 
    ha = 320 
    ha_prime = 350 
    if phase_n[0] < 0.25: 
        beta_ret = phase_n[0] 
    else: 
        beta_ret = 0.25-0.25*phase_n[0] 
    a_prime0 = phase_n[2] 
    b_0 = phase_n[0] 
 
    for i in range(1, TData.shape[0]-1): 
        beta_phase.append(phase_n[0]) 
        alpha_phase.append(phase_n[1]) 
        alpha_prime_phase.append(phase_n[2]) 
        TDot = cal_TDot(TData[i], TData[i+1], tData[i], tData[i+1]) 

        #interpolate for k and n 
        k1_real = interpolate(TKN_ref1[:,0], TKN_ref1[:,1], TData[i])
        n1_real = interpolate(TKN_ref1[:,0], TKN_ref1[:,2], TData[i]) 
        k2_real = interpolate(TKN_ref2[:,0], TKN_ref2[:,1], TData[i]) 
        n2_real = interpolate(TKN_ref2[:,0], TKN_ref2[:,2], TData[i]) 

        #condition 
        if (TData[i] > liquidus): 
            phase_n = [0, 0, 0, 0] 
            continue 
        elif (TData[i-1] > liquidus) and (TData[i] < liquidus): 
            phase_n = [1, 0, 0, 0] #purely beta phase 
            continue 
        elif TData[i] > T_beta_tran: 
            phase_n = [1, 0, 0, 0] 
            continue 
        else: 
            # Heating 
            if TDot > 0: 
                if phase_n[2] != 0: 
                    if TData[i] > Bs: 
                        alpha_equili, beta_equili, alpha_prime_equili, phase_n[2], phase_n[3] = \ 
                        equil(TData[i], alpha_equili, beta_equili, alpha_prime_equili, phase_n[2], phase_n[3], T_beta_tran) 
                        phase_n = x_alpha_prime(tData[i],tData[i+1],k1_real,n1_real,alpha_prime_equili,phase_n) 
                        continue 
                elif (phase_n[0] < beta_equili): 
                    phase_n = no_x_alpha_prime(tData[i],tData[i+1],k2_real,n2_real,beta_equili,alpha_equili,phase_n)  
                    continue 
                else: 
                    alpha_equili, beta_equili, alpha_prime_equili, phase_n[2], phase_n[3] = \ 
                    equil(TData[i], alpha_equili, beta_equili, alpha_prime_equili, phase_n[2], phase_n[3], T_beta_tran) 
            else: 

            # Cooling 
                if TDot < -20: 
                    if i == 1: 
                        if phase_n[0] < 0.25: 
                            beta_ret = phase_n[0] 
                        else: 
                            beta_ret = 0.25-0.25*phase_n[0] 
                        a_prime0 = phase_n[2] 
                        b_0 = phase_n[0] 
                    elif (TData[i-1] > T_beta_tran) and (TData[i] <= T_beta_tran): 
                        if phase_n[0] < 0.25: 
                            beta_ret = phase_n[0] 
                        else: 
                            beta_ret = 0.25-0.25*phase_n[0] 
                        a_prime0 = phase_n[2] 
                        b_0 = phase_n[0] 
                    elif (TData[i] < TData[i-1]) and (TData[i] < TData[i+1]): 
                        if phase_n[0] < 0.25: 
                            beta_ret = phase_n[0] 
                        else: 
                            beta_ret = 0.25-0.25*phase_n[0] 
                        a_prime0 = phase_n[2] 
                        b_0 = phase_n[0] 
                    alpha_equili, beta_equili, alpha_prime_equili, phase_n[2], phase_n[3] = \ 
                    equil(TData[i], alpha_equili, beta_equili, alpha_prime_equili, phase_n[2], phase_n[3], T_beta_tran) 
                    phase_n = transformation(TData[i],TData[i+1],tData[i],tData[i+1],phase_n,beta_ret,a_prime0,b_0,beta_equili, Ms, gamma) 
                    continue 
                else: 
                    alpha_equili, beta_equili, alpha_prime_equili, phase_n[2], phase_n[3] = \ 
                    equil(TData[i], alpha_equili, beta_equili, alpha_prime_equili, phase_n[2], phase_n[3], T_beta_tran)   
                    phase_n = difussion(TData[i],TData[i+1],tData[i],tData[i+1],k2_real,n2_real,beta_equili,phase_n) 
                    continue 
    
    # Calculate final hardness value
    hard_ness_value = hb*beta_phase[-1] + ha*alpha_phase[-1] + ha_prime*alpha_prime_phase[-1] 
    return alpha_phase,beta_phase,alpha_prime_phase, hard_ness_value 
