from scipy.optimize import brentq
from MortenStorage import get_policy_2_storage
from shortcuts import *

def get_policy_2_storage_modified(mismatch, eta_in=1., eta_out=1., CS=NaN, P_in=None, P_out=None):
    """ This function behaves like Morten's get_policy_2_storage(). However, it allows limiting the charging and discharging capacities of the storage. """
    
    if P_in==None:
        ## No constraints on charging
        P_in = amax(mismatch)
    if P_out==None:
        ## No constraints on discharging
        P_out = -amin(mismatch)
        
    ## The mismatch is cut by the maximum charging and discharging capacities.
    mismatch_ = lambda P_in, P_out: amax([amin([mismatch,P_in*ones_like(mismatch)],axis=0),-P_out*ones_like(mismatch)],axis=0) 
    
    ## The cut mismatch is feed into Mortens code to produce the cut reduced mismatch.     
    mismatch_r_, CS_used = get_policy_2_storage(mismatch_(P_in,P_out),eta_in,eta_out,CS)   
        
    ## Calculate the real reduced mismatch.
    mismatch_r = mismatch_r_ + (mismatch - mismatch_(P_in,P_out))
        
    ## REturns the reduced mismatch and the storage energy capacity that whas actually used (0<=CS_used<=CS).
    return mismatch_r, CS_used
    
def get_min_storage_cap_alt(L, GW, GS, gamma=1, alpha_w=1., CS=None, gain=.9,returnall=False,eta_charge=1.,eta_discharge=1.):
    """Finds upper and lower limits of the storage in and output capacities. acc sets the relative increase in balancing energy that is acceptable.
    WARNING! This function is slow to evaluate."""

    acc = 1-gain

    L, GW, GS, gamma, alpha_w = array(L,ndmin=2), array(GW,ndmin=2), array(GS,ndmin=2), array(gamma,ndmin=1), array(alpha_w,ndmin=1)  #Ensure minimum dimension to 2 to alow the weighed sum to be calculated correctly.
    if len(alpha_w)==1:
        alpha_w = alpha_w*ones_like(gamma)
    
    weighed_sum = lambda x: sum(x,axis=0)/mean(sum(x,axis=0))

    l = weighed_sum(L)
    Gw = weighed_sum(GW)	
    Gs = weighed_sum(GS)

    mismatch = lambda alpha_w, gamma: gamma*(alpha_w*Gw + (1.-alpha_w)*Gs) - l
    mismatch_ = lambda alpha_w, gamma, P_in, P_out: amax([amin([mismatch(alpha_w, gamma),P_in*ones_like(l)],axis=0),-P_out*ones_like(l)],axis=0)

    ## Summed residual load after storage.
    res_load_sum = lambda alpha_w, gamma, P_in, P_out: sum(get_positive(-get_policy_2_storage(mismatch_(alpha_w,gamma,P_in,P_out),storage_capacity = CS)[0]) + get_positive(-mismatch(alpha_w, gamma)-P_out))
    
    ## Difference between residual load before and after storage.
    storage_benefit = lambda alpha_w, gamma, P_in, P_out: sum(get_positive(-mismatch(alpha_w, gamma))) - res_load_sum(alpha_w,gamma,P_in,P_out)
    
    Storage_benefit, E_surplus, E_residual, E_discharge, E_charge, N_cycles = zeros_like(gamma), zeros_like(gamma), zeros_like(gamma), zeros_like(gamma), zeros_like(gamma), zeros_like(gamma)
    P_in_ = zeros((len(gamma),2))
    P_out_ = zeros((len(gamma),2))
    for i in arange(len(gamma)):
        Storage_benefit[i] = storage_benefit(alpha_w[i],gamma[i],1e3,1e3) 

        ## Values to be returned
        E_surplus[i] = mean(get_positive(mismatch(alpha_w, gamma))) ## Mean hourly surplus before storage.
        E_residual[i] = mean(get_positive(-mismatch(alpha_w, gamma))) ## Mean hourly residual load before storage. 
        E_discharge[i] = (1.-acc)*Storage_benefit[i]/len(l) ## Mean hourly output of the storage. (almost same as 'Storage_benefit')
        E_charge[i] = E_discharge/eta_discharge/eta_charge ## Mean hourly input of the storage.
        N_cycles[i] = (E_discharge/eta_discharge)/CS ## Average number of accumulated storage cycles per hour.

        print i,
        sys.stdout.flush()    

        if acc==0:  #Intended for quick calculation where the charging and discharging powers are not relevant.
            P_in_[i] = [NaN,NaN]
            P_out_[i] = [NaN,NaN]
        #elif Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], 1e-8, 1e3) > acc*Storage_benefit[i]:
        elif ((1-acc)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], 1e-8, 10) > 0)*((1-acc)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], 10, 1e-8) > 0):

            P_in_0 = brentq(lambda P_in: (1-acc/2.)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], P_in, 1e3), 1e-8, 10)
            P_out_0 = brentq(lambda P_out: (1-acc)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], P_in_0, P_out), 1e-8, 10)

            P_out_1 = brentq(lambda P_out: (1-acc/2.)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], 1e3, P_out), 1e-8, 10)
            P_in_1 = brentq(lambda P_in: (1-acc)*Storage_benefit[i] - storage_benefit(alpha_w[i], gamma[i], P_in, P_out_1), 1e-8, 10)

            P_in_[i] = [P_in_0, P_in_1]
            P_out_[i] = [P_out_0, P_out_1]
            
        else:
            P_in_[i] = [0,0]
            P_out_[i] = [0,0]

    if returnall==True:
        return Storage_benefit, P_in_, P_out_, E_surplus, E_residual, E_discharge, E_charge, N_cycles
    else:
        return Storage_benefit, P_in_, P_out_