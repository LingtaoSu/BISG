############################
import os
os.chdir('/home/dlg/sulin/librfn-python-master/cBioportal/')
rawpath = "GSE32062_zscore_2.txt"
TMuppath = 'GSE32062_TM_up.txt'
TMdnpath = 'GSE32062_TM_dn.txt'
TMup_1, TMdn_1 = data_Process(rawpath, TMuppath, TMdnpath)


############
ExpGenepath="GSE3494_Gene.txt";
ExpSampath="GSE3494_Sample.txt";
ExpMpath="GSE3494_TM_up.txt";
ExpGene=np.loadtxt(ExpGenepath, dtype='str')
ExpSam=np.loadtxt(ExpSampath, dtype='str')
ExpM = np.loadtxt(ExpMpath, np.float32);
data = _load_dataset('GSE3494_Survival.csv')

for i in range(1,16):
    #print("iter: ",i)
    Bicgenepath="/%s_Gene.txt" % i;
    Bicgene = np.loadtxt(Bicgenepath, dtype='str')
    #print(len(Bicgene))
    p2, ax = multilogrank_test5(ExpGene, ExpSam, Bicgene, ExpM, data);
    if(p2<0.05):
        print(i," :  ",p2)
        plt.savefig("Survival_%s.pdf" % i)
        plt.close();
    #print (p2)

###################### BISG Compare analysis #############

# for breast cancer, we use two GEO datasets GSE1456 from their paper and GSE3494, test the survival curve of 68 core genesets and 4-gene
# from our results.
ExpGenepath="GSE32062_Gene.txt";
ExpSampath="GSE32062_Sample.txt";
ExpMpath="GSE32062_TM_up.txt";
ExpGene=np.loadtxt(ExpGenepath, dtype='str')
ExpSam=np.loadtxt(ExpSampath, dtype='str')
ExpM = np.loadtxt(ExpMpath, np.float32);
data = _load_dataset('GSE32062_Survival.csv')

for i in range(1,9):
    #print("iter: ",i)
    Bicgenepath="4_%s_Gene.txt" % i;
    Bicgene = np.loadtxt(Bicgenepath, dtype='str')
    #print(len(Bicgene))
    p2, ax = multilogrank_test5(ExpGene, ExpSam, Bicgene, ExpM, data);
    if(p2<1.05):
        print(i," :  ",p2)
        plt.savefig("Survival_4_%s.pdf" % i)
        plt.close();
    print (p2)
#######################

