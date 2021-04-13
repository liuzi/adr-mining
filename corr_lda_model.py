import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymc3 as pm
import theano
import theano.tensor as tt
from utils._tools import write2file
from os.path import join

# print(f"Running on PyMC3 v{pm.__version__}")

def get_sample_data(n_patients=10,n_feature_pres=32,n_feature_diag=36):
    # simulated test data
    pres_data = np.random.randint(2,size=(n_patients,n_feature_pres))
    diag_data = np.random.randint(2,size =(n_patients,n_feature_diag))

    return pres_data, diag_data

def p_vector(p,n):
    return p * np.ones((1,n))

### Random sample zdi from zd for each patient
def pymc_corrlda_zdrnd(pres_data, diag_data):
    D,S = pres_data.shape
    D, tS = diag_data.shape 

    K =10
    tK=15

    alpha =0.01
    b_0,b_1=0.01, 0.01
    gamma=0.01
    tb_0,tb_1=0.01,0.01

    with pm.Model() as model:
        
        
    #     p_alpha = alpha * np.ones((1, K))
    #     p_beta = pm.Beta('p_beta', alpha=b_0, beta=b_1, shape=(K, S))
        
        theta = pm.Dirichlet("theta", a=p_vector(alpha,K), shape=(D, K))
    #     phi = pm.distributions.Dirichlet('phi', a=p_beta, shape=(K, S))
        phi = pm.Beta('phi', alpha=b_0, beta=b_1, shape=(K, S))
        psi = pm.Dirichlet("psi",a=p_vector(gamma,tK),shape=(K,tK))
        tphi = pm.Beta("tphi", alpha=tb_0, beta=tb_1, shape=(tK, tS))
        
        # For each DRUG <theta>: D*K
        zs = [pm.Categorical("z_d{}".format(d), p=theta[d], shape=S) for d in range(D)]
        ws = [pm.Bernoulli("w_{}_{}".format(d,i), p=phi[zs[d][i]], observed=pres_data[d][i]) 
            for d in range(D) for i in range(S)]
        
        # For each disease random.choice(zs[d]) <psi>: K*tK
        tzs = [
            pm.Categorical("tz_d{}".format(d), \
                p=psi[zs[d][np.random.randint(S)]], shape=tS) for d in range(D)]
        tws = [pm.Bernoulli("tw_{}_{}".format(d,i), p=tphi[tzs[d][i]], observed=diag_data[d][i]) 
            for d in range(D) for i in range(tS)]

# def p_vector(p,n):
#     return p * np.ones((1,n))


def pij(zsd,phi,s):
    return phi[zsd[s]][s]

## get actual drug topic for each document: mat D*tS*D
def nd3_2(mat,D,tS):
    return tt.reshape([mat[d,:,d] for d in range(D)],(D,tS))

### uniform get zdi from zd for each disease in each drug
def pymc_corrlda_zduniform(pres_data, diag_data):
    D, S = pres_data.shape
    D, tS = diag_data.shape

    K =10
    tK=15

    alpha =0.01
    b_0,b_1=0.01, 0.01
    gamma=0.01
    tb_0,tb_1=0.01,0.01

    with pm.Model() as model:
        
        theta = pm.Dirichlet("theta", a=p_vector(alpha,K), shape=(D, K))
        phi = pm.Beta('phi', alpha=b_0, beta=b_1, shape=(K, S))
        psi = pm.Dirichlet("psi",a=p_vector(gamma,tK),shape=(K,tK))
        tphi = pm.Beta("tphi", alpha=tb_0, beta=tb_1, shape=(tK, tS))
        
        ## For each DRUG <theta>: D*K,  sample(theta_d)->column d
        zs=pm.Categorical("zs",p=theta,shape=(S,D))
        ws=[pm.Bernoulli("ws_{}_{}".format(d,s),p=pij(zs.T[d],phi,s),observed=pres_data[d][s]) for d in range(D) for s in range(S)]
        
        ## For each disease random.choice(zs[d]) <psi>: K*tK
        sp_zs_idx = pm.DiscreteUniform("super_z", lower=np.zeros([D,tS]), upper=np.ones([D,tS])*S,shape=(D,tS))
        ## Map idx to z_s, each z_s is vector with length of D, s is from 1 to S
        super_zs = nd3_2(zs[sp_zs_idx],D,tS) #shape:(D*tS*D) -> (D,tS)
        
        tzs = pm.Categorical("tzs", p=psi[super_zs],shape=(D,tS)) 
        tws = [pm.Bernoulli("tw_{}_{}".format(d,ts), p=pij(tzs[d],tphi,ts), observed=diag_data[d][ts]) 
            for d in range(D) for ts in range(tS)]
    
    n_sample_iter=200
    write_prefix="/data/liu/mimic3/LDA_MODEL/CORR_LDA"
    with model:
    #     step = pm.metropolis.CategoricalGibbsMetropolis(zs+ws+tzs+tws)
        trace=pm.sample(n_sample_iter, init=None, cores=2)
        write2file(pd.DataFrame(trace["phi"][n_sample_iter-1]), join(write_prefix,"corrlda_presphi"))
        write2file(pd.DataFrame(trace['theta'][n_sample_iter-1]),join(write_prefix,"corrlda_prestheta"))
        write2file(pd.DataFrame(trace["psi"][n_sample_iter-1]),join(write_prefix,"corrlda_diagtheta"))
        write2file(pd.DataFrame(trace["tphi"][n_sample_iter-1]),join(write_prefix,"corrlda_diagphi"))

# def run_pymc_corrlda(z_distribution="uniform"):
#     n_sample_iter=200
#     write_prefix="/data/liu/mimic3/LDA_MODEL/CORR_LDA"
#     with model:
#     #     step = pm.metropolis.CategoricalGibbsMetropolis(zs+ws+tzs+tws)
#         trace=pm.sample(n_sample_iter, init=None, cores=2)
#         write2file(pd.DataFrame(trace["phi"][n_sample_iter-1]), join(write_prefix,"corrlda_presphi"))
#         write2file(pd.DataFrame(trace['theta'][n_sample_iter-1]),join(write_prefix,"corrlda_prestheta"))
#         write2file(pd.DataFrame(trace["psi"][n_sample_iter-1]),join(write_prefix,"corrlda_diagtheta"))
#         write2file(pd.DataFrame(trace["tphi"][n_sample_iter-1]),join(write_prefix,"corrlda_diagphi"))

pymc_corrlda_zduniform(*get_sample_data())