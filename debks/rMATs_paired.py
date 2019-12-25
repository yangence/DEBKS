#!/usr/bin/env python
# -*- coding:utf-8 -*-
# import libraries
import subprocess
import re,os,warnings,math,scipy,itertools,logging,time,datetime,pysam,traceback
from scipy import stats
from multiprocessing import Pool
from scipy.optimize import fmin_cobyla
from scipy.optimize import fmin_l_bfgs_b
from math import log
import pandas as pd
import numpy as np
from numpy import *
warnings.filterwarnings('ignore');
np.seterr(divide='ignore',invalid='ignore')
rho=0.9
np.random.seed(1231);
########## functions ############
#binomial MLE optimization functions
def logit(x):
    if x<0.01:
        x=0.01;
    if x>0.99:
        x=0.99;
    return(log(x/(1-x)));
    
def veclogit(x):
    res=[];
    for i in x:
        res.append(logit(i));
    return(res);

#START: log Marginal distribution for the paired analysis
#this is the full log likelihood function within the Laplace approximiation of the marginal distribution of paired replicates
#for one exon only
def myfunc_marginal_full_paired(x, *args):
    beta1=x[0];beta2=x[1];rho=x[2];var1=x[3];var2=x[4];
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    sigma=matrix([[var1,rho*sqrt(var1)*sqrt(var2)],[rho*sqrt(var1)*sqrt(var2),var2]]);
    sigma_inv=sigma.getI();
    vec=matrix([[logit(psi1)-logit(beta1)],[logit(psi2)-logit(beta2)]]);
    new_psi1=inclusion_length*psi1/(inclusion_length*psi1+skipping_length*(1-psi1));
    new_psi2=inclusion_length*psi2/(inclusion_length*psi2+skipping_length*(1-psi2));
    ##print('vec');#print(vec);
    temp1=0;temp2=0;temp3=0;temp4=0;
    temp1+=-0.5*log(abs(linalg.det(sigma)));
    temp2+=-0.5*dot(dot(vec.getT(),sigma_inv),vec);
    temp3+=I1*log(new_psi1)+S1*log(1-new_psi1)+log(psi1)+log(1-psi1);
    temp4+=I2*log(new_psi2)+S2*log(1-new_psi2)+log(psi2)+log(1-psi2);
    return(temp1+temp2+temp3+temp4);

#return the det of sigma1 matrix for one exon
def myfunc_marginal_sigma1_paired(x, *args):
    beta1=x[0];beta2=x[1];rho=x[2];var1=x[3];var2=x[4];
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    sigma=matrix([[var1,rho*sqrt(var1)*sqrt(var2)],[rho*sqrt(var1)*sqrt(var2),var2]]);
    sigma_inv=sigma.getI();
    vec=matrix([[logit(psi1)-logit(beta1)],[logit(psi2)-logit(beta2)]]);
    vec_der=matrix([[1/psi1/(1-psi1)],[1/psi2/(1-psi2)]]);
    vec_2der=matrix([[(-1+2*psi1)/pow(psi1,2)/pow(1-psi1,2)],[(-1+2*psi2)/pow(psi2,2)/pow(1-psi2,2)]]);
    new_psi1=inclusion_length*psi1/(inclusion_length*psi1+skipping_length*(1-psi1));
    new_psi2=inclusion_length*psi2/(inclusion_length*psi2+skipping_length*(1-psi2));
    a_apostrophe=sigma_inv[0,0]*pow(vec_der[0,0],2)+(sigma_inv[0,0]*vec[0,0]+sigma_inv[0,1]*vec[1,0])*vec_2der[0,0];
    a_apostrophe+=-1*pow(psi1,-2)+pow(1-psi1,-2);
    a_apostrophe+=-1*I1*skipping_length*((2*inclusion_length+skipping_length)*psi1+skipping_length*(1-psi1))/pow(psi1,2)/pow(inclusion_length*psi1+skipping_length*(1-psi1),2);
    a_apostrophe+=-1*S1*inclusion_length*((inclusion_length+2*skipping_length)*(1-psi1)+inclusion_length*psi1)/pow(1-psi1,2)/pow(inclusion_length*psi1+skipping_length*(1-psi1),2);
    b_apostrophe=sigma_inv[0,1]*vec_der[0,0]*vec_der[1,0];
    d_apostrophe=sigma_inv[1,1]*pow(vec_der[1,0],2)+(sigma_inv[1,1]*vec[1,0]+sigma_inv[0,1]*vec[0,0])*vec_2der[1,0];
    d_apostrophe+=-1*pow(psi2,-2)+pow(1-psi2,-2);
    d_apostrophe+=-1*I2*skipping_length*((2*inclusion_length+skipping_length)*psi2+skipping_length*(1-psi2))/pow(psi2,2)/pow(inclusion_length*psi2+skipping_length*(1-psi2),2);
    d_apostrophe+=-1*S2*inclusion_length*((inclusion_length+2*skipping_length)*(1-psi2)+inclusion_length*psi2)/pow(1-psi2,2)/pow(inclusion_length*psi2+skipping_length*(1-psi2),2);
    sigma_1=matrix([[a_apostrophe,b_apostrophe],[b_apostrophe,d_apostrophe]]);
    return(-1*sigma_1);

def myfunc_marginal_paired_rho0(x, *args):
    beta1=x[0];beta2=x[1];rho=0;
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    var1=args[8];var2=args[9];
    sum=0;temp1=matrix([[1]]);temp2=0;
    sum_temp1=0;sum_temp2=0;
    for i in range(len(psi1)):
        temp1=myfunc_marginal_sigma1_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        temp2=myfunc_marginal_full_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        sum+=(-0.5*log(abs(linalg.det(temp1))+1)+temp2);
        sum_temp1+=temp1;sum_temp2+=temp2;
    return(-1*sum);

def myfunc_marginal_der_paired_rho0(x, *args):
    beta1=x[0];beta2=x[1];rho=0;
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    var1=args[8];var2=args[9];
    sum=array([0.0,0.0]);
    sum_temp1=array([0.0,0.0]);sum_temp2=array([0.0,0.0]);
    temp1=0;temp2=0;
    for i in range(len(psi1)):
        temp1=myfunc_marginal_sigma1_der_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length,1);
        temp2=myfunc_marginal_full_der_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        sum_temp1+=temp1[0:2];sum_temp2+=temp2[0:2];
        sum=sum+temp1[0:2]+temp2[0:2];
        sigma1=myfunc_marginal_sigma1_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        #if len(args)>10:
        #    ofile=open('test_output.txt','a');ofile.write(str(temp1[0])+'\t'+str(temp1[1])+'\t'+str(psi1[i])+'\t'+str(psi2[i])+'\t'+str(linalg.det(sigma1))+'\t'+str(-0.5/(linalg.det(sigma1)+1))+'\t'+str(list(sigma1))+'\n');
    return(-1*array(sum));    
    
def myfunc_marginal_paired_rho(x, *args):
    beta1=x[0];beta2=x[1];rho=x[2];
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    var1=args[8];var2=args[9];
    sum=0;temp1=matrix([[1]]);temp2=0;
    sum_temp1=0;sum_temp2=0;
    for i in range(len(psi1)):
        temp1=myfunc_marginal_sigma1_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        temp2=myfunc_marginal_full_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        sum+=(-0.5*log(abs(linalg.det(temp1))+1)+temp2);
    return(-1*sum);

def myfunc_marginal_der_paired_rho(x, *args):
    beta1=x[0];beta2=x[1];rho=x[2];
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    var1=args[8];var2=args[9];
    sum=array([0.0,0.0,0.0]);
    sum_temp1=array([0.0,0.0,0.0,0.0,0.0]);sum_temp2=array([0.0,0.0,0.0,0.0,0.0]);
    temp1=0;temp2=0;
    for i in range(len(psi1)):
        temp1=myfunc_marginal_sigma1_der_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        temp2=myfunc_marginal_full_der_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        sum=sum+(temp1[0:3]+temp2[0:3]);
    return(-1*array(sum));
    
def myfunc_marginal_paired_fixvar(x, *args):
    beta1=x[0];beta2=x[1];
    #beta1=x[0];beta2=args[11];
    #beta2=x[0];beta1=args[11];
    rho=args[8];var1=args[9];var2=args[10];
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    sum=0;temp1=matrix([[1]]);temp2=0;
    sum_temp1=0;sum_temp2=0;
    ##print('fixvar_beta1');#print(beta1);#print('fixvar_beta2');#print(beta2);
    for i in range(len(psi1)):
        temp1=myfunc_marginal_sigma1_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        temp2=myfunc_marginal_full_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        sum+=(0.5*log(abs(linalg.det(temp1))+1)-temp2);
        sum_temp1+=temp1;sum_temp2+=temp2;
    ##print('sum_temp1');#print(sum_temp1);#print('sum_temp2');#print(sum_temp2);#print('sum');#print(sum);
    return(sum);

def myfunc_marginal_der_paired_fixvar(x, *args):
    beta1=x[0];beta2=x[1];
    #beta1=x[0];beta2=args[11];
    #beta2=x[0];beta1=args[11];
    rho=args[8];var1=args[9];var2=args[10];
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    sum=array([0.0,0.0]);
    #sum_temp1=array([0.0,0.0,0.0,0.0,0.0]);sum_temp2=array([0.0,0.0,0.0,0.0,0.0]);
    temp1=0;temp2=0;
    for i in range(len(psi1)):
        temp1=myfunc_marginal_sigma1_der_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        temp2=myfunc_marginal_full_der_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        #sum_temp1=sum_temp1-1*temp1;sum_temp2=sum_temp2-1*temp2;
        sum=sum+(-1*temp1[0:2]-1*temp2[0:2]);
    ##print('sum_temp1_der');#print(sum_temp1);#print('sum_temp2_der');#print(sum_temp2);
    return(array(sum));

def myfunc_marginal_paired(x, *args):
    beta1=x[0];beta2=x[1];rho=x[2];var1=x[3];var2=x[4];
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    sum=0;temp1=matrix([[1]]);temp2=0;
    for i in range(len(psi1)):
        temp1=myfunc_marginal_sigma1_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        temp2=myfunc_marginal_full_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        sum+=(-0.5*log(abs(linalg.det(temp1))+1)+temp2);
    return(-1*sum);
#END: log Marginal distribution for the paired analysis

#START: Derivative of the log Marginal distribution for the paired analysis
def myfunc_marginal_full_der_paired(x, *args):
    beta1=x[0];beta2=x[1];rho=x[2];var1=x[3];var2=x[4];
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    sigma=matrix([[var1,rho*sqrt(var1)*sqrt(var2)],[rho*sqrt(var1)*sqrt(var2),var2]]);
    sigma_inv=sigma.getI();
    vec=matrix([[logit(psi1)-logit(beta1)],[logit(psi2)-logit(beta2)]]);
    temp=-1*dot(transpose(vec),dot(sigma_inv,matrix([[-1/beta1/(1-beta1),0],[0,-1/beta2/(1-beta2)]])));
    temp1=temp[0,0];temp2=temp[0,1];
    diag=(-2*rho-1+pow(rho,2))/sqrt(var1)/sqrt(var2);
    temp3=-0.5*dot(transpose(vec),dot(matrix([[2*rho/var1,diag],[diag,2*rho/var2]]),vec))/pow(1-rho*rho,2)-0.5/linalg.det(sigma)*(-1*var1*var2);
    #temp4=(-1/pow(var1,2)*pow(vec[0,0],2)+rho/pow(var1,1.5)/pow(var2,0.5)*vec[0,0]*vec[1,0])/(1-rho*rho)+(1-rho*rho)*var2/linalg.det(sigma);
    #temp5=(-1/pow(var2,2)*pow(vec[1,0],2)+rho/pow(var2,1.5)/pow(var1,0.5)*vec[0,0]*vec[1,0])/(1-rho*rho)+(1-rho*rho)*var1/linalg.det(sigma);
    return(array([float(temp1),float(temp2),float(temp3)]));

def myfunc_marginal_sigma1_der_paired(x, *args):
    beta1=x[0];beta2=x[1];rho=x[2];var1=x[3];var2=x[4];
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    sigma=matrix([[var1,rho*sqrt(var1)*sqrt(var2)],[rho*sqrt(var1)*sqrt(var2),var2]]);
    sigma_inv=sigma.getI();
    ##print('sigma_inv');#print(sigma_inv);
    vec=matrix([[logit(psi1)-logit(beta1)],[logit(psi2)-logit(beta2)]]);
    vec_der=matrix([[1/psi1/(1-psi1)],[1/psi2/(1-psi2)]]);
    vec_2der=matrix([[(-1+2*psi1)/pow(psi1,2)/pow(1-psi1,2)],[(-1+2*psi2)/pow(psi2,2)/pow(1-psi2,2)]]);
    sigma1=myfunc_marginal_sigma1_paired([beta1,beta2,rho,var1,var2],I1,S1,psi1,I2,S2,psi2,inclusion_length,skipping_length);
    temp1=sigma_inv[0,0]*(-1/beta1/(1-beta1))*vec_2der[0,0]*sigma1[1,1]*(-1)+sigma_inv[0,1]*(-1/beta1/(1-beta1))*vec_2der[1,0]*sigma1[0,0]*(-1);
    temp2=sigma_inv[1,1]*(-1/beta2/(1-beta2))*vec_2der[1,0]*sigma1[0,0]*(-1)+sigma_inv[0,1]*(-1/beta2/(1-beta2))*vec_2der[0,0]*sigma1[1,1]*(-1);
    #temp1=(-1/var1/beta1/(1-beta1)*vec_2der[0,0])*sigma1[1,1]/(1-rho*rho)*(-1)-(rho/sqrt(var1)/sqrt(var2)/beta1/(1-beta1)*vec_2der[1,0])*sigma1[0,0]/(1-rho*rho);
    #temp2=(-1/var2/beta2/(1-beta2)*vec_2der[1,0])*sigma1[0,0]/(1-rho*rho)*(-1)-(rho/sqrt(var1)/sqrt(var2)/beta2/(1-beta2)*vec_2der[0,0])*sigma1[1,1]/(1-rho*rho);
    
    temp3=(-1/sqrt(var1)/sqrt(var2)*vec[1,0]*vec_2der[0,0])*sigma1[1,1]/(1-rho*rho)*(-1);
    temp3+=(-1/sqrt(var1)/sqrt(var2)*vec[0,0]*vec_2der[1,0])*sigma1[0,0]/(1-rho*rho)*(-1);
    temp3+=-2*(-1/sqrt(var1)/sqrt(var2)*vec_der[0,0]*vec_der[1,0])*sigma1[0,1]/(1-rho*rho)*(-1);
    temp3+=4*rho/(1-rho*rho)*linalg.det(sigma1);
    ##print('test_der');#print(linalg.det(sigma1));#print(x);#print(-0.5/(linalg.det(sigma1))*array([float(temp1),float(temp2),float(temp3)]));
    #if len(args)>8:
    #    ofile=open('test_output_2.txt','a');ofile.write(str([float(temp1),float(temp2),float(temp3)])+'\n');
    if linalg.det(sigma1)>0:
        return(-0.5/(linalg.det(sigma1)+1)*array([float(temp1),float(temp2),float(temp3)]));
    else:
        return(-0.5/(linalg.det(sigma1)-1)*array([float(temp1),float(temp2),float(temp3)]));
        
def myfunc_marginal_der_paired(x, *args):
    beta1=x[0];beta2=x[1];rho=x[2];var1=x[3];var2=x[4];
    I1=args[0];S1=args[1];psi1=args[2];
    I2=args[3];S2=args[4];psi2=args[5];
    inclusion_length=args[6];skipping_length=args[7];
    sum=array([0.0,0.0,0.0,0.0,0.0]);
    sum_temp1=array([0.0,0.0,0.0,0.0,0.0]);sum_temp2=array([0.0,0.0,0.0,0.0,0.0]);
    temp1=0;temp2=0;
    for i in range(len(psi1)):
        temp1=myfunc_marginal_sigma1_der_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        temp2=myfunc_marginal_full_der_paired([beta1,beta2,rho,var1,var2],I1[i],S1[i],psi1[i],I2[i],S2[i],psi2[i],inclusion_length,skipping_length);
        sum_temp1=sum_temp1-1*temp1;sum_temp2=sum_temp2-1*temp2;
        sum=sum+(-1*temp1-1*temp2);
    return(array(sum[0:5]));
#END: Derivative of the log Marginal distribution for the paired analysis
            

#START: log-likelihood function of each individual replicate
def myfunc_individual_paired(x, *args):
    psi1=x[0];psi2=x[1];
    I1=args[0];S1=args[1];beta1=args[2];var1=args[3];
    I2=args[4];S2=args[5];beta2=args[6];var2=args[7];
    inclusion_length=args[8];skipping_length=args[9];
    rho=args[10];
    temp=myfunc_marginal_full_paired([beta1,beta2,rho,var1,var2],I1,S1,psi1,I2,S2,psi2,inclusion_length,skipping_length);
    return(-1*temp);  

def myfunc_individual_der_paired(x,*args):
    psi1=x[0];psi2=x[1];
    I1=args[0];S1=args[1];beta1=args[2];var1=args[3];
    I2=args[4];S2=args[5];beta2=args[6];var2=args[7];
    inclusion_length=args[8];skipping_length=args[9];
    rho=args[10];
    sigma=matrix([[var1,rho*sqrt(var1)*sqrt(var2)],[rho*sqrt(var1)*sqrt(var2),var2]]);
    sigma_inv=sigma.getI();
    vec=matrix([[logit(psi1)-logit(beta1)],[logit(psi2)-logit(beta2)]]);
    vec_der=matrix([[1/psi1/(1-psi1)],[1/psi2/(1-psi2)]]);
    new_psi1=inclusion_length*psi1/(inclusion_length*psi1+skipping_length*(1-psi1));
    new_psi1_der=inclusion_length*skipping_length/pow(inclusion_length*psi1+skipping_length*(1-psi1),2);
    new_psi2=inclusion_length*psi2/(inclusion_length*psi2+skipping_length*(1-psi2));
    new_psi2_der=inclusion_length*skipping_length/pow(inclusion_length*psi2+skipping_length*(1-psi2),2);
    temp1=0;temp2=0;
    temp1+=-1*(sigma_inv[0,0]*vec[0,0]*vec_der[0,0]+sigma_inv[0,1]*vec_der[0,0]*vec[1,0]);
    temp1+=1*(I1/new_psi1*new_psi1_der-S1/(1-new_psi1)*new_psi1_der+1/psi1-1/(1-psi1));
    temp2+=-1*(sigma_inv[1,1]*vec[1,0]*vec_der[1,0]+sigma_inv[0,1]*vec_der[1,0]*vec[0,0]);
    temp2+=1*(I2/new_psi2*new_psi2_der-S2/(1-new_psi2)*new_psi2_der+1/psi2-1/(1-psi2));
    return(-1*array([temp1,temp2]));
#END: log-likelihood function of each individual replicate
    
def myfunc_multivar(x,*args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
    ##print('psi1');#print(psi1);#print('psi2');#print(psi2);
    sum1=0;sum2=0;
    for i in psi1:
        sum1+=pow(logit(i)-logit(x[0]),2);
    sum1=sum1/var1/2;
    for i in psi2:
        sum2+=pow(logit(i)-logit(x[1]),2);
    sum2=sum2/var2/2;
    return(sum1+sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x[0]),2)+pow(stats.norm.ppf(x[1]),2)-2*rho*stats.norm.ppf(x[0])*stats.norm.ppf(x[1])));

def myfunc_multivar_der(x,*args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
    sum1=0;sum2=0;
    for i in psi1:
        sum1+=-2*(logit(i)-logit(x[0]))/x[0]/(1-x[0]);
    sum1=sum1/var1/2;
    res1=sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x[0])-2*rho*stats.norm.ppf(x[1]))/stats.norm.pdf(stats.norm.ppf(x[0]));
    ##print('1');#print(x[1]);#print(res1);#print(stats.norm.pdf(stats.norm.ppf(x[0])));
    for i in psi2:
        sum2+=-2*(logit(i)-logit(x[1]))/x[1]/(1-x[1]);
    sum2=sum2/var2/2;
    res2=sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x[1])-2*rho*stats.norm.ppf(x[0]))/stats.norm.pdf(stats.norm.ppf(x[1]));
    return(np.array([res1,res2]));

def myfunc_1(x, *args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
    sum1=0;sum2=0;
    for i in psi1:
        sum1+=pow(logit(i)-logit(x+cutoff),2);
    sum1=sum1/var1/2;
    for i in psi2:
        sum2+=pow(logit(i)-logit(x),2);
    sum2=sum2/var2/2;
    return(sum1+sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x+cutoff),2)+pow(stats.norm.ppf(x),2)-2*rho*stats.norm.ppf(x+cutoff)*stats.norm.ppf(x)));

def myfunc_der_1(x, *args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
    sum1=0;sum2=0;
    for i in psi1:
        sum1+=-2*(logit(i)-logit(x+cutoff))/(x+cutoff)/(1-x-cutoff);
    sum1=sum1/var1/2;
    res1=sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x+cutoff)-2*rho*stats.norm.ppf(x))/stats.norm.pdf(stats.norm.ppf(x+cutoff));
    for i in psi2:
        sum2+=-2*(logit(i)-logit(x))/x/(1-x);
    sum2=sum2/var2/2;
    res2=sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x)-2*rho*stats.norm.ppf(x+cutoff))/stats.norm.pdf(stats.norm.ppf(x));
    return(res1+res2);

def myfunc_2(x, *args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
    sum1=0;sum2=0;
    for i in psi1:
        sum1+=pow(logit(i)-logit(x),2);
    sum1=sum1/var1/2;
    for i in psi2:
        sum2+=pow(logit(i)-logit(x+cutoff),2);
    sum2=sum2/var2/2;
    return(sum1+sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(pow(stats.norm.ppf(x+cutoff),2)+pow(stats.norm.ppf(x),2)-2*rho*stats.norm.ppf(x+cutoff)*stats.norm.ppf(x)));

def myfunc_der_2(x, *args):
    psi1=args[0];psi2=args[1];var1=args[2];var2=args[3];
    sum1=0;sum2=0;
    for i in psi1:
        sum1+=-2*(logit(i)-logit(x))/(x)/(1-x);
    sum1=sum1/var1/2;
    res1=sum1+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x+cutoff)-2*rho*stats.norm.ppf(x))/stats.norm.pdf(stats.norm.ppf(x+cutoff));
    for i in psi2:
        sum2+=-2*(logit(i)-logit(x+cutoff))/(x+cutoff)/(1-x-cutoff);
    sum2=sum2/var2/2;
    res2=sum2+0.1*0.5*(pow(rho,2))/(1-pow(rho,2))*(2*stats.norm.ppf(x)-2*rho*stats.norm.ppf(x+cutoff))/stats.norm.pdf(stats.norm.ppf(x));
    return(res1+res2);

def myfunc_marginal_2der(x, args):
    I=args[0];S=args[1];beta=args[2];var=args[3];
    inclusion_length=args[4];
    skipping_length=args[5]
    temp1=1/pow(x,2)/pow(1-x,2)*(((2*x-1)*(logit(beta)-logit(x))-1)/var-1);#temp1=(2*x-1)/pow(x,2)/pow(1-x,2)*((logit(beta)-logit(x)-1/(2*x-1))/var+1);
    temp2=I*skipping_length*((2*inclusion_length+skipping_length)*x+skipping_length*(1-x))/pow(x,2)/pow(inclusion_length*x+skipping_length*(1-x),2);
    temp3=S*inclusion_length*((inclusion_length+2*skipping_length)*(1-x)+inclusion_length*x)/pow(1-x,2)/pow(inclusion_length*x+skipping_length*(1-x),2);
    ##print('test');#print(beta);#print(x);#print(var);#print(temp1);#print(temp2);#print(temp3);
    return(temp1-temp2-temp3);
    
def myfunc_marginal(x, *args):
    beta=x;
    I=args[0];S=args[1];psi=args[2];var=args[3];
    inclusion_length=args[4];skipping_length=args[5];
    sum=0;temp1=0;temp2=0;
    for i in range(len(psi)):
        new_psi=inclusion_length*psi[i]/(inclusion_length*psi[i]+skipping_length*(1-psi[i]));
        f1=I[i]*log(new_psi)+S[i]*log(1-new_psi)-pow(logit(psi[i])-logit(beta),2)/2/var-log(psi[i])-log(1-psi[i]);
        temp2+=-1*f1;
        f1_2der=abs(myfunc_marginal_2der(psi[i],[I[i],S[i],beta,var,inclusion_length,skipping_length]));
        temp1+=0.5*log(abs(f1_2der));
        sum+=(0.5*log(abs(f1_2der)+0.00001)-f1);
    ##print('test_temp1');#print(temp1);#print(temp2);
    return(sum);

def myfunc_marginal_der(x, *args):
    beta=x;
    I=args[0];S=args[1];psi=args[2];var=args[3];
    inclusion_length=args[4];skipping_length=args[5];
    sum=0;
    for i in range(len(psi)):
        new_psi=inclusion_length*psi[i]/(inclusion_length*psi[i]+skipping_length*(1-psi[i]));
        f1_3der=1*(2*psi[i]-1)/pow(psi[i],2)/pow(1-psi[i],2)/beta/(1-beta)/var;
        f1_2der=myfunc_marginal_2der(psi[i],[I[i],S[i],beta,var,inclusion_length,skipping_length]);
        ##print('test2');#print(f1_2der);
        f1_der=(logit(psi[i])-logit(beta))/beta/(1-beta)/var;
        f1=I[i]*log(new_psi)+S[i]*log(1-new_psi)-pow(logit(psi[i])-logit(beta),2)/2/var+log(psi[i])+log(1-psi[i]);
        sum+=(0.5*f1_3der/(f1_2der)-f1_der);
    return(sum);

def myfunc_marginal_1(x, *args):
    beta2=x;beta1=x+cutoff;
    I1=args[0];S1=args[1];psi1=args[2];var1=args[3];
    I2=args[4];S2=args[5];psi2=args[6];var2=args[7];
    inclusion_length=args[8];skipping_length=args[9];
    return(myfunc_marginal(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length));

def myfunc_marginal_1_der(x, *args):
    beta2=x;beta1=x+cutoff;
    I1=args[0];S1=args[1];psi1=args[2];var1=args[3];
    I2=args[4];S2=args[5];psi2=args[6];var2=args[7];
    inclusion_length=args[8];skipping_length=args[9];
    return(myfunc_marginal_der(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal_der(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length));

def myfunc_marginal_2(x, *args):
    beta2=x+cutoff;beta1=x;
    I1=args[0];S1=args[1];psi1=args[2];var1=args[3];
    I2=args[4];S2=args[5];psi2=args[6];var2=args[7];
    inclusion_length=args[8];skipping_length=args[9];
    return(myfunc_marginal(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length));

def myfunc_marginal_2_der(x, *args):
    beta2=x+cutoff;beta1=x;
    I1=args[0];S1=args[1];psi1=args[2];var1=args[3];
    I2=args[4];S2=args[5];psi2=args[6];var2=args[7];
    inclusion_length=args[8];skipping_length=args[9];
    return(myfunc_marginal_der(beta1,I1,S1,psi1,var1,inclusion_length,skipping_length)+myfunc_marginal_der(beta2,I2,S2,psi2,var2,inclusion_length,skipping_length));
    
def myfunc_individual(x,*args):
    I=args[0];S=args[1];beta=args[2];var=args[3];
    inclusion_length=args[4];
    skipping_length=args[5]
    new_psi=inclusion_length*x/(inclusion_length*x+skipping_length*(1-x));
    return(-1*(I*log(new_psi)+S*log(1-new_psi)-(logit(x)-logit(beta))*(logit(x)-logit(beta))/2/var-log(x)-log(1-x)));

def myfunc_individual_der(x,*args):
    I=args[0];S=args[1];beta=args[2];var=args[3];
    inclusion_length=args[4];
    skipping_length=args[5];
    new_psi=inclusion_length*x/(inclusion_length*x+skipping_length*(1-x));
    new_psi_der=inclusion_length*skipping_length/pow(inclusion_length*x+skipping_length*(1-x),2);
    return(-1*(I/new_psi*new_psi_der-S/(1-new_psi)*new_psi_der-(logit(x)-logit(beta))/var/x/(1-x)-1/x+1/(1-x) ));

def myfunc_likelihood(x, args):
    I=args[0];S=args[1];beta=args[2];var=args[3];sum=0;N=I+S;
    #return(-1*(-log(sqrt((I+S)*x*(1-x)))-(I-(I+S)*x)*(I-(I+S)*x)/2/((I+S)*x*(1-x))-log(sqrt(var))-(x-beta)*(x-beta)/2/var));
    ##print('debug');#print(N);#print(var);#print(x);#print(beta);
    if N==0:
        return(0);
    return(-0.5*((I-N*x)*(I-N*x)/(N*x)+(S-N*(1-x))*(S-N*(1-x))/(N*(1-x)))-log(sqrt(var))-(x-beta)*(x-beta)/2/var);

def MLE_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
    psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
    beta_0=sum(psi1)/len(psi1);
    beta_1=sum(psi2)/len(psi2);
    var1=10*scipy.var(np.array(psi1)-beta_0);
    var2=10*scipy.var(np.array(psi2)-beta_1);
    if var1<=0.01:
        var1=0.01;
    if var2<=0.01:
        var2=0.01;
    #print('var1');#print(var1);#print('var2');#print(var2);
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1;
        #iteration of beta
        beta_0=sum(psi1)/len(psi1);
        beta_1=sum(psi2)/len(psi2);
        #print('var1');#print(var1);#print('var2');#print(var2);
        #if abs(sum(psi1)/len(psi1)-sum(psi2)/len(psi2))>cutoff:
        if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_1,[sum(psi2)/len(psi2)],myfunc_der_1,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            theta2 = max(min(float(xopt[0]),1-cutoff),0);theta1=theta2+cutoff;
        else:#minize psi1 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_2,[sum(psi1)/len(psi1)],myfunc_der_2,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            theta1 = max(min(float(xopt[0]),1-cutoff),0);theta2=theta1+cutoff;
        #print('constrain_1xopt');#print('theta');#print(theta1);#print(theta2);#print(xopt);
        #else:
        #    theta1=sum(psi1)/len(psi1);theta2=sum(psi2)/len(psi2);
        beta_0=theta1;beta_1=theta2;
        #iteration of psi
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
        #print('constrain_2xopt');
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt);
            #likelihood_sum+=myfunc_marginal(new_psi1[i],[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length]);
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt);
            #likelihood_sum+=myfunc_marginal(new_psi2[i],[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length]);
        #print('new_psi[0]');#print(new_psi1[0]);#print(new_psi2[0]);
        psi1=new_psi1;psi2=new_psi2;
        #print('count');#print(count);#print('previous_sum');#print(previous_sum);#print('current_sum');#print(current_sum);
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum);
        previous_sum=current_sum;
    ##print('constrain');#print(theta1);#print(theta2);#print(psi1);#print(psi2);#print(current_sum);#print(likelihood_sum);
    ##print(xopt);
    return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);
    #return([likelihood_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);

def MLE_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
    psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
    beta_0=sum(psi1)/len(psi1);
    beta_1=sum(psi2)/len(psi2);
    var1=10*scipy.var(np.array(psi1)-beta_0);
    var2=10*scipy.var(np.array(psi2)-beta_1);
    if var1<=0.01:
        var1=0.01;
    if var2<=0.01:
        var2=0.01;
    #print('var1');#print(var1);#print('var2');#print(var2);
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1;
        #iteration of beta
        beta_0=sum(psi1)/len(psi1);
        beta_1=sum(psi2)/len(psi2);
        xopt=fmin_l_bfgs_b(myfunc_multivar,[beta_0,beta_1],myfunc_multivar_der,args=[psi1,psi2,var1,var2],bounds=[[0.01,0.99],[0.01,0.99]],iprint=-1);
        beta_0=float(xopt[0][0]);
        beta_1=float(xopt[0][1]);
        #print('unconstrain_1xopt');#print(xopt);
        #print('theta');#print(beta_0);#print(beta_1);#print('theta_end');
        #iteration of psi
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt);
            #likelihood_sum+=myfunc_marginal(new_psi1[i],[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length]);
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt);
            #likelihood_sum+=myfunc_marginal(new_psi2[i],[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length]);
        #print('new_psi[0]');#print(new_psi1[0]);#print(new_psi2[0]);
        psi1=new_psi1;psi2=new_psi2;print
        #print('count');#print(count);('previous_sum');#print(previous_sum);#print('current_sum');#print(current_sum);
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum);
        previous_sum=current_sum;
    ##print('unconstrain');#print(beta_0);#print(beta_0+beta_1);#print(psi1);#print(psi2);#print(current_sum);#print(likelihood_sum);
    ##print(xopt);
    if count>iter_maxrun:
        return([current_sum,[psi1,psi2,0,0,var1,var2]]);
    return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);
    #return([likelihood_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);

def MLE_marginal_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
    #initial value
    psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
    beta_0=sum(psi1)/len(psi1);
    beta_1=sum(psi2)/len(psi2);
    #var1=scipy.var(np.array(veclogit(psi1)));
    #var2=scipy.var(np.array(veclogit(psi2)));
    var1=scipy.var(np.array(psi1)-beta_0);
    var2=scipy.var(np.array(psi2)-beta_1);
    if var1<=0.01:
        var1=0.01;
    if var2<=0.01:
        var2=0.01;
    #MLE of the full likelihood
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1;
        beta_0=sum(psi1)/len(psi1);
        beta_1=sum(psi2)/len(psi2);
        xopt=fmin_l_bfgs_b(myfunc_multivar,[beta_0,beta_1],myfunc_multivar_der,args=[psi1,psi2,var1,var2],bounds=[[0.01,0.99],[0.01,0.99]],iprint=-1);
        beta_0=float(xopt[0][0]);
        beta_1=float(xopt[0][1]);
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);
        psi1=new_psi1;psi2=new_psi2;
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum);
        previous_sum=current_sum;
    if count>iter_maxrun:
        psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
        iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
        beta_0=sum(psi1)/len(psi1);
        beta_1=sum(psi2)/len(psi2);
        #var1=scipy.var(np.array(veclogit(psi1)));
        #var2=scipy.var(np.array(veclogit(psi2)));
        var1=10*scipy.var(np.array(psi1)-beta_0);
        var2=10*scipy.var(np.array(psi2)-beta_1);
        if var1<=0.01:
            var1=0.01;
        if var2<=0.01:
            var2=0.01;
        while((iter_cutoff>0.01)&(count<=iter_maxrun)):
            count+=1;
            beta_0=sum(psi1)/len(psi1);
            beta_1=sum(psi2)/len(psi2);
            xopt=fmin_l_bfgs_b(myfunc_multivar,[beta_0,beta_1],myfunc_multivar_der,args=[psi1,psi2,var1,var2],bounds=[[0.01,0.99],[0.01,0.99]],iprint=-1);
            beta_0=float(xopt[0][0]);
            beta_1=float(xopt[0][1]);
            new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
            for i in range(len(psi1)):
                xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
                new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);
            for i in range(len(psi2)):
                xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
                new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);
            psi1=new_psi1;psi2=new_psi2;
            if count>1:
                iter_cutoff=abs(previous_sum-current_sum);
            previous_sum=current_sum;
            if count>iter_maxrun:
                return([current_sum,[psi1,psi2,0,0,var1,var2,0]]);
    
    #MLE of the marginal likelihood
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1;
        ##print('unconstrain_MLE_marginal_value_1');#print(myfunc_marginal(0.51,i1,s1,psi1,var1,effective_inclusion_length,effective_skipping_length));
        ##print('unconstrain_MLE_marginal_value_2');#print(myfunc_marginal(0.56,i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length));
        #exit();
        xopt1=fmin_l_bfgs_b(myfunc_marginal,[beta_0],myfunc_marginal_der,args=[i1,s1,psi1,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
        xopt2=fmin_l_bfgs_b(myfunc_marginal,[beta_1],myfunc_marginal_der,args=[i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
        beta_0=float(xopt1[0][0]);
        beta_1=float(xopt2[0][0]);
        #print(xopt1)
        ##if re.findall('\d+',str(xopt1[2]))[0]!='0':
        #    return([current_sum,[psi1,psi2,0,0,var1,var2,0]]);
        ##if re.findall('\d+',str(xopt2[2]))[0]!='0':
        #    return([current_sum,[psi1,psi2,0,0,var1,var2,0]]);
        #print('unconstrain_1xopt');print(xopt1);
        #print('unconstrain_2xopt');print(xopt2);
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);
        ##if re.findall('\d+',str(xopt[2]))[0]!='0':
            #return([current_sum,[psi1,psi2,0,0,var1,var2,0]]);
        psi1=new_psi1;psi2=new_psi2;
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum);
        previous_sum=current_sum;
    cor=max(min(stats.pearsonr(psi1,psi2)[0],0.9),-0.9);
    #print(cor)
    if math.isnan(cor):
        cor=0;
    if count>iter_maxrun:
        return([current_sum,[psi1,psi2,0,0,var1,var2,cor]]);
    return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2,cor]]);

def MLE_marginal_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length):
    #initial value
    psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
    beta_0=sum(psi1)/len(psi1);
    beta_1=sum(psi2)/len(psi2);
    #var1=scipy.var(np.array(veclogit(psi1)));
    #var2=scipy.var(np.array(veclogit(psi2)));
    var1=scipy.var(np.array(psi1)-beta_0);
    var2=scipy.var(np.array(psi2)-beta_1);
    if var1<=0.01:
        var1=0.01;
    if var2<=0.01:
        var2=0.01;
    #MLE of the full likelihood
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1;
        beta_0=sum(psi1)/len(psi1);
        beta_1=sum(psi2)/len(psi2);
        if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_1,[sum(psi2)/len(psi2)],myfunc_der_1,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            theta2 = max(min(float(xopt[0]),1-cutoff),0);theta1=theta2+cutoff;
        else:#minize psi1 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_2,[sum(psi1)/len(psi1)],myfunc_der_2,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            theta1 = max(min(float(xopt[0]),1-cutoff),0);theta2=theta1+cutoff;
        beta_0=theta1;beta_1=theta2;
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt);
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt);
        psi1=new_psi1;psi2=new_psi2;
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum);
        previous_sum=current_sum;
    if count>iter_maxrun:
        psi1=vec2psi(i1,s1,effective_inclusion_length,effective_skipping_length);psi2=vec2psi(i2,s2,effective_inclusion_length,effective_skipping_length);
        iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
        beta_0=sum(psi1)/len(psi1);
        beta_1=sum(psi2)/len(psi2);
        #var1=scipy.var(np.array(veclogit(psi1)));
        #var2=scipy.var(np.array(veclogit(psi2)));
        var1=10*scipy.var(np.array(psi1)-beta_0);
        var2=10*scipy.var(np.array(psi2)-beta_1);
        if var1<=0.01:
            var1=0.01;
        if var2<=0.01:
            var2=0.01;
        #MLE of the full likelihood
        while((iter_cutoff>0.01)&(count<=iter_maxrun)):
            count+=1;
            beta_0=sum(psi1)/len(psi1);
            beta_1=sum(psi2)/len(psi2);
            if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
                xopt = fmin_l_bfgs_b(myfunc_1,[sum(psi2)/len(psi2)],myfunc_der_1,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
                theta2 = max(min(float(xopt[0]),1-cutoff),0);theta1=theta2+cutoff;
            else:#minize psi1 if this is the case
                xopt = fmin_l_bfgs_b(myfunc_2,[sum(psi1)/len(psi1)],myfunc_der_2,args=[psi1,psi2,var1,var2],bounds=[[0.001,0.999-cutoff]],iprint=-1)
                theta1 = max(min(float(xopt[0]),1-cutoff),0);theta2=theta1+cutoff;
            beta_0=theta1;beta_1=theta2;
            new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
            for i in range(len(psi1)):
                xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
                new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt);
            for i in range(len(psi2)):
                xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
                new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt);
            psi1=new_psi1;psi2=new_psi2;
            if count>1:
                iter_cutoff=abs(previous_sum-current_sum);
            previous_sum=current_sum;
        if count>iter_maxrun:
            return([current_sum,[psi1,psi2,0,0,var1,var2]]);

    
    #MLE of the marginal likelihood
    iter_cutoff=1;iter_maxrun=100;count=0;previous_sum=0;
    while((iter_cutoff>0.01)&(count<=iter_maxrun)):
        count+=1;
        if (sum(psi1)/len(psi1))>(sum(psi2)/len(psi2)):#minize psi2 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_marginal_1,[beta_1],myfunc_marginal_1_der,args=[i1,s1,psi1,var1,i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            beta_1 = max(min(float(xopt[0]),1-cutoff),0);beta_0=beta_1+cutoff;
        else:#minize psi1 if this is the case
            xopt = fmin_l_bfgs_b(myfunc_marginal_2,[beta_0],myfunc_marginal_2_der,args=[i1,s1,psi1,var1,i2,s2,psi2,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.001,0.999-cutoff]],iprint=-1)
            beta_0 = max(min(float(xopt[0]),1-cutoff),0);beta_1=beta_0+cutoff;
        #if re.findall('\d+',str(xopt[2]))[0]!='0':
            #return([current_sum,[psi1,psi2,0,0,var1,var2]]);
        #print('constrain_xopt');#print(xopt);
        new_psi1=[];new_psi2=[];current_sum=0;likelihood_sum=0;
        for i in range(len(psi1)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi1[i]],myfunc_individual_der,args=[i1[i],s1[i],beta_0,var1,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi1.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt);
        for i in range(len(psi2)):
            xopt = fmin_l_bfgs_b(myfunc_individual,[psi2[i]],myfunc_individual_der,args=[i2[i],s2[i],beta_1,var2,effective_inclusion_length,effective_skipping_length],bounds=[[0.01,0.99]],iprint=-1);
            new_psi2.append(float(xopt[0]));current_sum+=float(xopt[1]);#print(xopt);
        #if re.findall('\d+',str(xopt[2]))[0]!='0':
            #return([current_sum,[psi1,psi2,0,0,var1,var2]]);
        psi1=new_psi1;psi2=new_psi2;
        if count>1:
            iter_cutoff=abs(previous_sum-current_sum);
        previous_sum=current_sum;
    if count>iter_maxrun:
        return([current_sum,[psi1,psi2,0,0,var1,var2]]);

    return([current_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);
    #return([likelihood_sum,[psi1,psi2,beta_0,beta_1,var1,var2]]);
    


    
#Random Sampling Function
def likelihood_test(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length,flag,id):
    #print('testing'+str(id));
    if flag==0:
        return([1,0]);
    else:
        res=MLE_marginal_iteration(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length);
        #print(res)
        if abs(res[1][2]-res[1][3])<=cutoff:
            #print('1<=cutoff');#print(res);#print((res[1][2]-res[1][3]));
            #return([1-scipy.stats.norm.cdf(abs(res[1][2]-res[1][3]),cutoff,sqrt(res[1][4]+res[1][5])),res[0]]);
            return([1,0]);
        else:
            res_constrain=MLE_marginal_iteration_constrain(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length);
            #print('2>cutoff');#print(res);#print(res_constrain);
            #print(abs(res_constrain[0]-res[0]));#print(pow(1-pow(res[1][6],2),2));#print('2end');
            ##print(abs(res_constrain[2]-res[2]));#print('2end_marginal');            
            return([scipy.stats.chi2.sf(2*(abs(res_constrain[0]-res[0]))/pow(1-pow(res[1][6],2),2),1),res[1][-1]]);

#MultiProcessorFunction
def MultiProcessorPool(n_original_diff):
    i1=n_original_diff[0];i2=n_original_diff[1];s1=n_original_diff[2];s2=n_original_diff[3];
    effective_inclusion_length=n_original_diff[4];effective_skipping_length=n_original_diff[5];
    flag=n_original_diff[6];id=n_original_diff[7];
    P=likelihood_test(i1,i2,s1,s2,effective_inclusion_length,effective_skipping_length,flag,id);
    return(P);

#Function for vector handling
def vec2float(vec):
    res=[];
    for i in vec:
        res.append(float(i));
    return(res);

#add 1 in both inclusion and skipping counts for robustness in small sample size
def vecAddOne(vec):
    res=[];
    for i in vec:
        res.append(i+1);
    return(res);

def vecprod(vec):
    res=1;
    for i in vec:
        res=res*i;
    return(res);

def vecadd(vec1,vec2):
    res=[];
    for i in range(len(vec1)):
        res.append(vec1[i]+vec2[i]);
    return(res);
    
def vec2remove0psi(inc,skp):
    res1=[];res2=[];
    for i in range(len(inc)):
        if (inc[i]!=0) | (skp[i]!=0):
            res1.append(inc[i]);res2.append(skp[i]);
    return([res1,res2]);

def vec2psi(inc,skp,effective_inclusion_length,effective_skipping_length):
    psi=[];
    inclusion_length=effective_inclusion_length;
    skipping_length=effective_skipping_length;
    for i in range(len(inc)):
        if (float(inc[i])+float(skp[i]))==0:
            psi.append(0.5);
        else:
            psi.append(float(inc[i])/inclusion_length/(float(inc[i])/inclusion_length+float(skp[i])/skipping_length));
    return(psi);

def vec210(vec):
    res=[];
    for i in vec:
        if i>0:
            res.append(1);
        else:
            res.append(-1);
    return(res);

#this function removes the paired samples with inc1+skp1=0 or inc2+skp2=0
def vec_remove_na_pair(inc1,skp1,inc2,skp2):
    res_inc1=[];res_skp1=[];res_inc2=[];res_skp2=[];
    for i in range(len(inc1)):
        if ((float(inc1[i])+float(skp1[i]))>0) & ((float(inc2[i])+float(skp2[i]))>0):
            res_inc1.append(inc1[i]);
            res_skp1.append(skp1[i]);
            res_inc2.append(inc2[i]);
            res_skp2.append(skp2[i]);
    return([res_inc1,res_skp1,res_inc2,res_skp2]);
    
def vec_remove_na(inc,skp):
    res_inc=[];res_skp=[];
    for i in range(len(inc)):
        if (float(inc[i])+float(skp[i]))>0:
            res_inc.append(inc[i]);
            res_skp.append(skp[i]);
    return([res_inc,res_skp]);
def myorder(p,reverse):
    p_dict={}
    ##print('len(p)');#print(len(p))
    for index in range(len(p)):
        if not(p[index] in p_dict):
            p_dict[p[index]]=[]
        p_dict[p[index]].append(index+1)
    ##print('p_dict');#print(p_dict)
    o_index=sorted(p_dict,reverse=reverse);o=[]
    ##print('o_index');#print(o_index)
    for index in o_index:
        for this_p in p_dict[index]:
            o.append(this_p)
    return(o)
def mycummin(p):
        res=[]
        for i in range(len(p)):
            res.append(min(p[:(i+1)]))
        return(res)

def myFDR(p):
    lp=len(p)
    i=range(lp,0,-1)
    o=myorder(p,True)
    ro=myorder(o,False)
    p_new=[]
    for index in range(len(o)):
        p_new.append(p[o[index]-1]*lp/i[index])
    p_mycummin=mycummin(p_new)
    p_mycummin_new=[]
    for index in range(len(p_mycummin)):
        p_mycummin_new.append(min(p_mycummin[index],1))
    res=[]
    for index in range(len(p_mycummin_new)):
        res.append(p_mycummin_new[ro[index]-1])
    return(res)

def adjustPsi(inc,skp,effective_inclusion_length,effective_skipping_length):
    #if effective_inclusion_length == 0:
        #print([inc,skp,effective_inclusion_length,effective_skipping_length])
    inc_adj=np.array(inc)/effective_inclusion_length
    skp_adj=np.array(skp)/effective_skipping_length
    return(','.join([str(i) for i in (inc_adj/(inc_adj+skp_adj)).tolist()]))

def myttest(vec1,vec2):
    if (len(vec1)==1) & (len(vec2)==1):
        res=stats.ttest_ind([vec1[0],vec1[0]],[vec2[0],vec2[0]]);
    else:
        res=stats.ttest_ind(vec1,vec2);
    return(res);

def myand(l1,l2):
    return([i for i in range(len(l1)) if l1[i] & l2[i] ==True])

    ######## end of functions ##############
def getrMATS_P(list_n_original_diff,mycutoff,threads):
    global cutoff
    cutoff=mycutoff
    probability=[]
    if threads>1:
        pool=Pool(processes=threads)
        probability=pool.map(MultiProcessorPool,list_n_original_diff)
        pool.close()
        pool.join()
    else:
        for i in range(len(list_n_original_diff)):
            probability.append(MultiProcessorPool(list_n_original_diff[i]))
    P=pd.DataFrame(probability).iloc[:,0].tolist()
    return(P)
