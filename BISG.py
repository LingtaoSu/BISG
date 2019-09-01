#################################################################################################
"""
Python module for BISG.Version 1.

Copyright: ((This work was supported by the National Nature Science Foundation of China (No. 61373051and NO. 61772226) and US National Institutes of Health (R35-GM126985)).

Lingtao Su(major author).sulin@missouri.edu.
Mingrui Ma(assistant).280406100@qq.com

Statement: this module uses librfn, downloaded from https://github.com/bioinf-jku/librfn-python,
    Copyright © 2015-2017 Thomas Unterthiner.

An example for using this module:
    you can set values to the parameters (learnrateW, learnratePsi, dropout_rate, minP, thr_p, thr_rp, thr_t, thr_t2, zero)
		BISG(GSMT, sample, character, GeneSnpMTGene, GeneSnpMTSnp, data, ResultFolder,
             TMup_1, n_hidden, n_iter,ITt=80,
             learnrateW=0.1, learnratePsi=0.1, dropout_rate=0.1, minP=1e-1,
             thr_p=0.05, thr_rp=0.5, thr_t=1.3, thr_t2=3, zero=0.8,
             GPUID=0)
	or you can use the default parameters:
        BISG(GSMT, sample, character, GeneSnpMT, GeneSnpMTGene, GeneSnpMTSnp, data, ResultFolder,
			 TMup_1, n_hidden, n_iter,ITt)

Parameters
----------
GSMT: Gene SNP map matrix, if a SNP located on a gene then the value is 1, else the value is 0.
sample: a sample list file contains all sample ids of gene expression matrix.
character: a gene/SNP list files contains all genes/SNPs of gene expression matrix.
GeneSnpMTGene: all genes in the GSMT matrix.
data: patients survival data.
ResultFolder: the folder path for the results.
TMup_1: binary matrix, 1 means up regulation,0 means no change.
n_hidden: the bicluster number, which need to set before running RFN.
n_iter: the iteration times of RFN.
learnrateW: learning rate of RFN
learnratePsi: parameters in RFN.
dropout_rate: parameters when do RFN modle training.
minP: RFN paramter.
thr_p: p-value use to filter biclusters.
thr_rp: p-value use to filter biclusters.
thr_t: threshold value use to filter bicluster samples.
thr_t2: threshold value use to filter bicluster genes.
zero: zero ratio in a bicluster, for bicluster quality measure.
GPUID: RFN parameter, to signify wheter use CPU or GPU
ITt: iteration times when do detecting survival related biclusers.


"""

### load package #######################
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import pandas as pd
from pkg_resources import resource_filename
kmf = KaplanMeierFitter()
import math
import pandas as pd
import numpy as np
from scipy import stats
from sklearn import preprocessing
from rfn import *
import random
from scipy.spatial import distance

################# methods defined##################
def _load_dataset(filename, **kwargs):
    """
    Load a dataset from lifelines.datasets

    Parameters
    ----------
    filename : string
        for example "larynx.csv"
    usecols : list
        list of columns in file to use

    Returns
    -------
        output: DataFrame
    """
    return pd.read_csv(resource_filename("lifelines", "datasets/" + filename), engine="python", **kwargs)




def multilogrank_test(Samindx,newsample,newdata, ST,thr_p): #ST sample times
    """ Test bicluster genes in survival classification during ST samplings """
    count=0;# count the number of significant test during ST sampling.
    for itr in range(ST):
        T = newdata["OS_MONTHS"]
        E = [0] * len(newsample)
        for n in range(0, len(Samindx)):
            loc = Samindx[n]
            E[loc] = 1
        Samindex2 = []  # Samples not included in bicluster
        for m in range(0, len(newsample)):
            if (m not in Samindx):
                Samindex2.append(m);
        saminx = random.sample(Samindex2, len(Samindx)) #error control by bicluster size.
        if(len(saminx)>2):# at least include two samples
            for k in range(0, len(saminx)):
                loc2 = saminx[k]
                E[loc2] = 2
            newdata["AGE"] = E
            f = newdata.AGE == 2
            T = newdata[f]['OS_MONTHS']
            C = newdata[f]['OS_STATUS']
            f2 = newdata.AGE == 1
            T2 = newdata[f2]['OS_MONTHS']
            C2 = newdata[f2]['OS_STATUS']
            results = logrank_test(T, T2, C, C2, alpha=0.99)
            p = results.p_value
            if (p < thr_p):
                count = count + 1;
        else:
            itr=itr-1;
    return float(count)/ST;



def multilogrank_test4(samp,MT,data,temchar,char,thr_p): # use to validate genes in a bicluster
    """ test whether genes can seperate patient groups by their significant different survival cures """
    rdsampleindex = random.sample(range(MT.shape[1]), 100)  # random select 100 samples no rep
    TestMT = MT[:, rdsampleindex]# Test set expression matrix
    Testsample = samp[rdsampleindex]
    Testdata = data[0:0]
    for id in range(0, len(Testsample)):
        tempsample = Testsample[id]
        lation = [x for x in range(len(data["PATIENT_ID"])) if data["PATIENT_ID"][x] == tempsample]
        Testdata = Testdata.append(data[lation[0]:(lation[0] + 1)])
    T3 = Testdata["OS_MONTHS"]
    E3 = [0] * len(Testsample);
    Samindex3 = []  # Samples not included in group 1
    for inx in range(0, len(Testsample)):
        ix = [i for i, v in enumerate(TestMT[:, inx] > 0) if v]
        tempChar = char[ix]  # Genes in sample inx
        temcount = len(np.intersect1d(tempChar, temchar))
        if (temcount > (len(temchar) * 0.8)):# sample include over 80% bicluster genes
            E3[inx] = 1;  # seperate all samples into two groups, 1 if the sample include 80% genes in the bicluster and 0 otherwise
        else:
            Samindex3.append(inx);
    saminx2 = random.sample(Samindex3, E3.count(1))
    for k in range(0, len(saminx2)):
        loc2 = saminx2[k]
        E3[loc2] = 2
    Testdata["AGE"] = E3
    f3 = Testdata.AGE == 2
    T3 = Testdata[f3]['OS_MONTHS']
    C3 = Testdata[f3]['OS_STATUS']
    f4 = Testdata.AGE == 1
    T4 = Testdata[f4]['OS_MONTHS']
    C4 = Testdata[f4]['OS_STATUS']
    results2 = logrank_test(T3, T4, C3, C4, alpha=0.99)
    p2 = results2.p_value
    groups = Testdata["AGE"];
    T6 = Testdata["OS_MONTHS"]
    E6 = Testdata['AGE'];
    ix = (groups == 1);
    kmf.fit(T6[~ix], E6[~ix], label='Without up-regulated genes');
    ax = kmf.plot();
    kmf.fit(T6[ix], E6[ix], label='With up-regulated genes');
    ax = kmf.plot(ax=ax);
    return p2, ax;

#### test gene sets classification significane  ###########

def multilogrank_test5(ExpGene,ExpSam,Bicgene,ExpM,data): #ST sample times# use to validate genes in a bicluster
    """ return the union of two lists """
    #rdsampleindex = random.sample(range(MT.shape[1]), 100)  # random select 100 sample no reptive
    #TestMT = MT[:, rdsampleindex]# Test set expression matrix
    TestMT =ExpM  # Expression matrix
    Testsample =ExpSam # Samples
    Testdata = data
    T3 = Testdata["OS_MONTHS"]
    E3 = [0] * len(Testsample);
    Samindex3 = []  # Samples not included in group 1
    for inx in range(0, len(Testsample)):
        # for inx in range(0, 100):
        ix = [i for i, v in enumerate(TestMT[:, inx] > 0) if v]
        # ix = [i for i, v in enumerate(MTout.transpose()[:, inx] >= 0.5) if v]
        #print(len(ix))
        temSNP = ExpGene[ix]  # SNPs in sample inx
        #print(temSNP)
        temcount = len(np.intersect1d(temSNP, Bicgene))
        #print(temcount)
        if (temcount > 2):
        #if (temcount > (len(Bicgene) * 0.5)):
            # if (temcount >= len(temchar)*0.5):  # this sample include all the SNPs in the bicluster
            E3[inx] = 1;  # seperate all samples into two groups 1 if the sample include all SNPs in the bicluster and 0 otherwise
        else:
            Samindex3.append(inx);
    #print(len(Samindex3))
    #print(E3);
    #print("exe here")
    rnloc=random.randint(0,len(E3)-1)# random selection a location in E3 and set its value to 1.
    E3[rnloc]=1;
    #print(E3)
    #print(E3.count(1))
    #if (E3.count(1) > 2):
    if(E3.count(1)>2 and E3.count(1)<(len(E3)/2)): #make sure there are enough sample for random sampling., However, when E3.count(1)==0, error will occur: TypeError: 'NoneType' object is not iterable. to resolve this, we simplely random set a location as 1. other solutions??
        #print("E3's 1 number : ", E3.count(1))
        saminx2 = random.sample(Samindex3, E3.count(1))
        #print(len(saminx2))
        for k in range(0, len(saminx2)):
            loc2 = saminx2[k]
            #print(loc2)
            E3[loc2] = 2
        #print(E3.count(2))
        Testdata["AGE"] = E3
        f3 = Testdata.AGE == 2
        T3 = Testdata[f3]['OS_MONTHS']
        C3 = Testdata[f3]['OS_STATUS']
        f4 = Testdata.AGE == 1
        T4 = Testdata[f4]['OS_MONTHS']
        C4 = Testdata[f4]['OS_STATUS']
        results2 = logrank_test(T3, T4, C3, C4, alpha=0.99)
        p2 = results2.p_value
        #print(p2)
        groups = Testdata["AGE"];
        T6 = Testdata["OS_MONTHS"]
        E6 = Testdata['AGE'];
        ix = (groups == 1);
        kmf.fit(T6[~ix], E6[~ix], label='Without the gene set');
        ax = kmf.plot();
        kmf.fit(T6[ix], E6[ix], label='With the gene set');
        ax = kmf.plot(ax=ax);
        plt.title("p-value: %s "% p2);
        #print(ax)
        return p2, ax;
    else:
        return 0,0;


###################


def data_Process(rawpath,TMuppath,TMdnpath): # generate discritized datasets
    rawM = np.loadtxt(rawpath, np.float32);
    TMup = np.zeros(rawM.shape, dtype=np.float32)  # up-regulated
    TMdn = np.zeros(rawM.shape, dtype=np.float32)  # down-regulated
    for r in range(rawM.shape[0]):
        TMr1 = np.quantile(rawM[r, :], 0.75)
        TMr2 = np.quantile(rawM[r, :], 0.25)
        for j in range(rawM.shape[1]):
            if (rawM[r, j] > TMr1 and rawM[r, j] > 0):  # only keep real high-up regulated genes
                TMup[r, j] = rawM[r, j];
            elif (rawM[r, j] < TMr2 and rawM[r, j] < 0):
                TMdn[r, j] = rawM[r, j];
    # transform Matrix to [0-1] values.
    TMup_1 = np.zeros(rawM.shape, dtype=np.float32)
    TMdn_1 = np.zeros(rawM.shape, dtype=np.float32)
    for r in range(rawM.shape[0]):
        TMr1 = np.max(TMup[r, :])  # get max value of each row
        TMr2 = np.min(TMdn[r, :])  # get min value of each row.
        for j in range(rawM.shape[1]):
            if (TMr1 > 0 and TMup[r, j] > 0):
                TMup_1[r, j] = TMup[r, j] / TMr1;
            if (TMr2 < 0 and TMdn[r, j] < 0):
                TMdn_1[r, j] = TMdn[r, j] / TMr2;
    with open(TMuppath, 'wb') as f:
        np.savetxt(f, TMup_1, fmt='%1f', delimiter='	')
    with open(TMdnpath, 'wb') as f:
        np.savetxt(f, TMdn_1, fmt='%1f', delimiter='	')
    f.close();
    TM=TMup_1+TMdn_1;
    return TMup_1, TMdn_1,TM;





def BISG(GSMT,sample,character,GeneSnpMTGene,GeneSnpMTSnp,data,TM,n_hidden,n_iter,learnrateW,learnratePsi,dropout_rate,minP,thr_p,thr_rp,thr_t,thr_t2,zero,ResultFolder): # main code of BISG
    char = np.loadtxt(character, dtype='str')
    samp = np.loadtxt(sample, dtype='str')
    MGene = np.loadtxt(GeneSnpMTGene, dtype='str')
    MSnp = np.loadtxt(GeneSnpMTSnp, dtype='str')
    MT = TM
    Bicp = [];  # use to storage p value
    Bicp2 = [];  # use to storage p2 value
    num = 0
    ST = 100  # sampling times
    E=7; # E=-ln(10-3), for bicluster statistical evaluation
    for y in range(ITt):
        print("Iter times:", y);
        newsampleindex = random.sample(range(MT.shape[1]), 100)  # random select 100 sample no reptive
        newMT = MT[:, newsampleindex]
        newsample = samp[newsampleindex]
        newdata = data[0:0]
        for id in range(0, len(newsample)):
            tempsample = newsample[id]
            lation = [x for x in range(len(data["PATIENT_ID"])) if data["PATIENT_ID"][x] == tempsample]
            newdata = newdata.append(data[lation[0]:(lation[0] + 1)])
        W, P, Wout = train_rfn(newMT.transpose(), n_hidden, n_iter, learnrateW, learnratePsi, minP, dropout_rate,
                               activation="relu", gpu_id=0)
        H = np.maximum(0, np.dot(Wout, newMT))
        for index in range(0, len(W)):
            temmax = np.max(W[index])
            if (temmax > 0.06):  # if max value in W <0 then no bicluster elements are found
                Snpindx = []  # storege SNP location in a bicluster
                for ix3 in range(0, len(W[index])):
                    temperct = temmax / np.abs(W[index][ix3])
                    if (temperct < thr_t and W[index][ix3] > 0.06):
                        Snpindx.append(ix3)
                temchar = char[Snpindx]
                if (len(Snpindx) > 2):  # at least two chars included in a bicluster
                    Samindx = []  # storege SNP location in a bicluster
                    temmax2 = np.max(H[index])  # max value in H
                    for ix2 in range(0, len(H[index])):
                        if (H[index][ix2] > 1.5):
                            temperct2 = temmax2 / H[index][ix2]
                            if (temperct2 < thr_t2):
                                Samindx.append(ix2)
                    temsam = newsample[Samindx]
                    if (len(Samindx) > 2):  # at least two samples included
                        BicM1 = newMT[Snpindx, :]
                        BicM = BicM1[:, Samindx]
                        Rato = (float(np.count_nonzero(BicM)) / (len(Snpindx) * len(Samindx)))  # none zero ratio in bicluster detected, quality of bicluster
                        #print(Rato)
                        Ratoall = (float(np.count_nonzero(newMT)) / (
                                    len(newMT) * len(newMT[0])));  # none-zero ratio of the input matrix
                        Cmnk = float(np.count_nonzero(BicM)) - (len(Snpindx) * len(Samindx) * Ratoall) - math.sqrt(
                            (3 * E) * (len(Snpindx) * len(Samindx) * Ratoall)) # statistical significance
                        if (Rato > zero and Cmnk>0):  # bicluster fullfill threshold
                            print(Rato, Ratoall, Cmnk)
                            rp = multilogrank_test(Samindx, newsample, newdata, ST, thr_p);
                            if (rp > thr_rp):
                                rp2, ax = multilogrank_test4(samp, MT, data, temchar, char, thr_p);
                                if (rp2 < thr_p):
                                    print("the bicluster fulfill threshold:")
                                    print(rp, rp2);
                                    Bicp.append(rp);
                                    Bicp2.append(rp2);
                                    num = num + 1
                                    temGene = [];  # use to storag Genes
                                    for inx1 in range(0, len(temchar)):
                                        Sloc = MGene.tolist().index(
                                            temchar[inx1]);  # get the location of Genes in temchar
                                        Gloc = GSMT[:, Sloc].tolist().index(1);  # get the location of SNP related gene
                                        TGene = MSnp[Gloc];
                                        temGene.append(TGene);
                                    with open(ResultFolder+'%s_SNP.txt' % num,'w') as filehandle:
                                        filehandle.writelines("%s\n" % place for place in temGene)
                                    with open(ResultFolder+'%s_Gene.txt' % num,'w') as filehandle:
                                        filehandle.writelines("%s\n" % place for place in temchar)
                                    with open(ResultFolder+'%s_Sample.txt' % num,'w') as filehandle:
                                        filehandle.writelines("%s\n" % place for place in temsam);
                                    plt.savefig(ResultFolder+'Survival_%s.pdf' % num)
                                    plt.close();
    with open(ResultFolder+'/Bicp.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % place for place in Bicp);
    with open(ResultFolder+'/Bicp2.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % place for place in Bicp2)
###################################





