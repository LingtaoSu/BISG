############### RUN CODE###############
import os
os.chdir('/home/dlg/sulin/librfn-python-master/cBioportal/')

cancerID=1;
rawpath = "%s_median_Zscores_1.txt" % cancerID
TMuppath = '%s_TM_up.txt' % cancerID;
TMdnpath = '%s_TM_dn.txt' % cancerID;
TM, TMup_1, TMdn_1 = data_Process(rawpath, TMuppath, TMdnpath)
for i in range (0,3):
    if (i==0):
        ITt=1000;
    else:
        ITt=1000*(5*i);
    sample = "%s_Gene_Exp_Sample.txt" % cancerID  #
    character = "%s_Gene_Exp_Gene.txt" % cancerID
    GeneSnpMT = "SNP_Gene_Matrix_%s.txt" % cancerID
    GeneSnpMTGene = "SNP_Sample_Matrix_%s_Gene.txt" % cancerID
    GeneSnpMTSnp = "SNP_Sample_Matrix_%s_SNP.txt" % cancerID
    data = _load_dataset('Survival_%s_3.csv' % cancerID)  # load survival data /home/dlg/sulin/librfn-python-master
    GSMT = np.loadtxt(GeneSnpMT, np.int)
    ResultFolder = '/cBioportalResults2/%s/%s/' % (cancerID, ITt);
    n_hidden = 50
    n_iter = 100
    learnrateW = 0.1
    learnratePsi = 0.10
    dropout_rate = 0.1
    minP = 1e-1
    # thr_s=0 # threshold to filter samples
    # thr_c=0.03 # threshold to filter characters
    thr_p = 0.05
    thr_rp = 0.5  # number of times significant during ST sampling
    # thr_d=15 # distance between matrix input and matrix output through RFN decompose
    thr_t = 1.3  # the times between max and each element
    thr_t2 = 3  # the threshold used to select element in H
    zero = 0.8  # zero ratio in a bicluster.
    BISG2(GSMT, sample, character, GeneSnpMTGene, GeneSnpMTSnp, data, TMup_1, n_hidden, n_iter,
              learnrateW, learnratePsi, dropout_rate, minP, thr_p, thr_rp, thr_t, thr_t2, zero, ResultFolder)

