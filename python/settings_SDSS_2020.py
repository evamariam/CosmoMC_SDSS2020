# sample settings for a particular grid run
import copy
from paramgrid import batchjob

# Directory to find .ini files
ini_dir = 'batch3/'

# directory to look for existing covariance matrices
cov_dir = 'DR16covmat/'

# ini files you want to base each set of runs on
defaults = ['common.ini']
getdist_options = {'ignore_rows': 0.3, 'marker[nrun]': 0,
                   'marker[w]': -1,'marker[wa]': 0,'marker[omegak]': 0}
importanceDefaults = ['importance_sampling.ini']



#dataset names
CMB = 'CMB'
SN = 'SN'
BAO = 'BAO'
RSD = 'RSD'
BAORSD = 'BAORSD'
SDSS = 'SDSS'
allGrowth = 'allGrowth'
WL = 'WL'
CMBLens = 'CMBLens'
JLA = 'JLA'
DR7 = 'DR7'
DES = 'DES'
WMAP = 'WMAP'
lens = 'lens'
BAOnoLya = 'BAOnoLya'

BAOdata = 'SDSS_final_BAO.ini'
RSDdata = 'SDSS_final_RSD.ini'
BAORSDdata = 'SDSS_final_BAORSD.ini'
SNdata = 'Pantheon18.ini'
CMBData =  ['plik_rd12_HM_v22_TTTEEE.ini', 'lowl.ini', 'simall_EE.ini']
CMBLensData =  ['plik_rd12_HM_v22_TTTEEE.ini', 'lowl.ini', 'simall_EE.ini','lensing.ini']
GrowthData = ['SDSS_final_RSD.ini','DES_lensing.ini','lensing.ini']
GrowthOnlyData = ['SDSS_final_RSD.ini','DES_lensing.ini','lensonly.ini']
WLdata = 'DES_lensing.ini'
lensData = 'lensonly.ini'
DESdata = 'DES.ini'
JLAdata = 'JLA.ini'
DR7data = 'SDSS_DR7_BAO.ini'
WMAPdata = 'WMAP.ini'
BAOnoLyaData = 'SDSS_final_BAO_noLya.ini'

noCMB = {'lmin_store_all_cmb': 0, 'lmin_computed_cl': 0,
         'get_sigma8': False}

# set up list of groups of parameters and data sets
groups = []

# make first group of runs (all parameter variations with all data combinations)
g = batchjob.jobGroup('main')
g.params = [[],['omegak'],['w'],['mnu']]
g.datasets = []
g.datasets.append(batchjob.dataSet([CMB], [CMBData]))
g.datasets.append(batchjob.dataSet([CMB,SN], [CMBData,SNdata]))
g.datasets.append(batchjob.dataSet([CMB,BAO], [CMBData, BAOdata]))
g.datasets.append(batchjob.dataSet([CMB,BAOnoLya], [CMBData, BAOnoLyaData]))
g.datasets.append(batchjob.dataSet([CMB,BAO,SN], [CMBData,BAOdata,SNdata]))
g.datasets.append(batchjob.dataSet([CMBLens,BAORSD], [CMBLensData,BAORSDdata]))
g.datasets.append(batchjob.dataSet([CMBLens,BAORSD,SN], [CMBLensData,BAORSDdata,SNdata]))
g.datasets.append(batchjob.dataSet([CMBLens,BAORSD,DES,SN], [CMBLensData,BAORSDdata,DESdata,SNdata]))
groups.append(g)

gG = batchjob.jobGroup('growth')
gG.params = [[],['omegak'],['w']]
gG.datasets = []
gG.datasets.append(batchjob.dataSet([CMB,RSD], [CMBData,RSDdata]))
gG.datasets.append(batchjob.dataSet([CMB,WL], [CMBData,WLdata]))
gG.datasets.append(batchjob.dataSet([CMB,allGrowth], [CMBData,GrowthData]))
groups.append(gG)

gDE = batchjob.jobGroup('extraDE')
gDE.params = [['w','wa'],['omegak','w'],['omegak','w','wa']]
gDE.datasets = []
gDE.datasets.append(batchjob.dataSet([CMB,BAO,SN], [CMBData,BAOdata,SNdata]))
gDE.datasets.append(batchjob.dataSet([CMB,BAOnoLya,SN], [CMBData,BAOnoLyaData,SNdata]))
gDE.datasets.append(batchjob.dataSet([CMBLens,BAORSD,DES,SN], [CMBLensData,BAORSDdata,DESdata,SNdata]))
gDE.datasets.append(batchjob.dataSet([CMBLens,BAORSD], [CMBLensData,BAORSDdata]))
gDE.datasets.append(batchjob.dataSet([CMBLens,SN], [CMBLensData,SNdata]))
groups.append(gDE)

gN = batchjob.jobGroup('extraNeutrino')
gN.params = [['mnu']]
gN.datasets = []
gN.datasets.append(batchjob.dataSet([CMBLens,BAORSD,DES], [CMBLensData,BAORSDdata,DESdata]))
gN.datasets.append(batchjob.dataSet([CMBLens], [CMBLensData]))
#gN.datasets.append(batchjob.dataSet([CMBLens,RSD], [CMBLensData,RSDdata]))
#gN.datasets.append(batchjob.dataSet([CMBLens,DES], [CMBLensData,DESdata]))
gN.datasets.append(batchjob.dataSet([CMBLens,SN], [CMBLensData,SNdata]))
gN.datasets.append(batchjob.dataSet([CMBLens,BAO], [CMBLensData,BAOdata]))
gN.datasets.append(batchjob.dataSet([CMBLens,BAO,SN], [CMBLensData,BAOdata,SNdata]))
groups.append(gN)

gNDE = batchjob.jobGroup('NeutrinoDE')
gNDE.params = [['mnu','w']]
gNDE.datasets = []
gNDE.datasets.append(batchjob.dataSet([CMBLens,BAORSD], [CMBLensData,BAORSDdata]))
gNDE.datasets.append(batchjob.dataSet([CMBLens,BAORSD,SN], [CMBLensData,BAORSDdata,SNdata]))
gNDE.datasets.append(batchjob.dataSet([CMBLens,BAORSD,DES], [CMBLensData,BAORSDdata,DESdata]))
gNDE.datasets.append(batchjob.dataSet([CMBLens,BAORSD,DES,SN], [CMBLensData,BAORSDdata,DESdata,SNdata]))
#gNDE.datasets.append(batchjob.dataSet([CMBLens], [CMBLensData]))
gNDE.datasets.append(batchjob.dataSet([CMBLens,SN], [CMBLensData,SNdata]))
gNDE.datasets.append(batchjob.dataSet([CMBLens,BAO], [CMBLensData,BAOdata]))
gNDE.datasets.append(batchjob.dataSet([CMBLens,BAO,SN], [CMBLensData,BAOdata,SNdata]))
groups.append(gNDE)


gDEC = batchjob.jobGroup('Decadal')
gDEC.params = [['omegak','w','mnu']]
gDEC.datasets = []
gDEC.datasets.append(batchjob.dataSet([WMAP,JLA,DR7], [WMAPdata,JLAdata,DR7data]))
gDEC.datasets.append(batchjob.dataSet([CMBLens,JLA,DR7], [CMBLensData,JLAdata,DR7data]))
gDEC.datasets.append(batchjob.dataSet([WMAP,SN,DR7], [WMAPdata,SNdata,DR7data]))
gDEC.datasets.append(batchjob.dataSet([WMAP,JLA,DR7,DES], [WMAPdata,JLAdata,DR7data,DESdata]))
gDEC.datasets.append(batchjob.dataSet([WMAP,JLA,BAORSD], [WMAPdata,JLAdata, BAORSDdata]))
gDEC.datasets.append(batchjob.dataSet([CMBLens,BAORSD,DES,SN], [CMBLensData,BAORSDdata,DESdata,SNdata]))
#gDEC.datasets.append(batchjob.dataSet([CMBLens,DES,SN], [CMBLensData,DESdata,SNdata]))
gDEC.datasets.append(batchjob.dataSet([WMAP,JLA,DR7,'BAOLya'], [WMAPdata,JLAdata,DR7data,'SDSS_final_BAO_Lya.ini']))
groups.append(gDEC)

BAOonlydata = batchjob.dataSet(['BAOonly'], [BAOdata, 'DR16_BAOonly_priors.ini', noCMB])

gBAO = batchjob.jobGroup('BAOonly')
gBAO.params = [[],['omegak'],['w']]
gBAO.datasets = [BAOonlydata]
gBAO.datasets.append(batchjob.dataSet(['BAOonlyNoLya'], [BAOnoLyaData, 'DR16_BAOonly_priors.ini', noCMB]))
gBAO.extra_opts ={'MPI_Converge_Stop': 0.002}
groups.append(gBAO)

gT = batchjob.jobGroup('BAOtracers')
gT.params = [[]]
gT.datasets.append(batchjob.dataSet(['BAOonlyLRG'], ['SDSS_final_BAO_LRG', 'DR16_BAOonly_priors.ini', noCMB]))
gT.datasets.append(batchjob.dataSet(['BAOonlyELG'], ['SDSS_final_BAO_ELG', 'DR16_BAOonly_priors.ini', noCMB]))
gT.datasets.append(batchjob.dataSet(['BAOonlyQSO'], ['SDSS_final_BAO_QSO', 'DR16_BAOonly_priors.ini', noCMB]))
gT.datasets.append(batchjob.dataSet(['BAOonlyLya'], ['SDSS_final_BAO_Lya', 'DR16_BAOonly_priors.ini', noCMB]))
gT.datasets.append(batchjob.dataSet(['BAOonlyDR12LRG'], ['SDSS_DR12_BAO_LRG', 'DR16_BAOonly_priors.ini', noCMB]))
gT.datasets.append(batchjob.dataSet(['BAOonlyMGS'], ['SDSS_MGS_BAO.ini', 'DR16_BAOonly_priors.ini', noCMB]))
gT.extra_opts ={'MPI_Converge_Stop': 0.002}
groups.append(gT)

gBBN = batchjob.jobGroup('BBN')
gBBN.params = [[]]
gBBN.datasets = [copy.deepcopy(BAOonlydata)]
gBBN.datasets.append(batchjob.dataSet(['BAOonlyLowz'], ['SDSS_final_BAO_lowz', 'DR16_BAOonly_priors.ini', noCMB]))
gBBN.datasets.append(batchjob.dataSet(['BAOonlyHighz'], ['SDSS_final_BAO_highz', 'DR16_BAOonly_priors.ini', noCMB]))
gBBN.datasets.append(batchjob.dataSet(['BAOonlyNoLya'], [BAOnoLyaData, 'DR16_BAOonly_priors.ini', noCMB]))
for d in gBBN.datasets:
    d.add('BBN','Cooke17BBN.ini')
groups.append(gBBN)

gH0 = batchjob.jobGroup('H0')
gH0.params = [[]]
gH0.datasets = [copy.deepcopy(BAOonlydata).add('H0','HST_Riess2019.ini')]
groups.append(gH0)

gH0lad = batchjob.jobGroup('H0ladder')
gH0lad.params = [[],['omegak','w','wa']]
gH0lad.datasets = [copy.deepcopy(BAOonlydata).add('H0','HST_Riess2019.ini').add(SN,SNdata)]
groups.append(gH0lad)


RSDonlyData = batchjob.dataSet(['RSD'], [RSDdata, noCMB])
BAORSDonlyData = batchjob.dataSet(['BAORSD'], [BAORSDdata, noCMB])

glp = batchjob.jobGroup('lenspriors')
glp.params = [[]]
glp.datasets = [RSDonlyData]
glp.datasets += [BAORSDonlyData]
glp.datasets.append(batchjob.dataSet([DES],[DESdata,noCMB]))
glp.datasets.append(batchjob.dataSet([WL],[WLdata,noCMB]))
glp.datasets.append(batchjob.dataSet([lens],[lensData,noCMB]))
glp.datasets.append(batchjob.dataSet([allGrowth],[GrowthOnlyData,noCMB]))
for d in glp.datasets:
    d.add('lenspriors','lensonly_priors.ini')
glp.extra_opts ={'get_sigma8': True}
groups.append(glp)

gSN = batchjob.jobGroup('SN')
gSN.params = [['w'],['omegak']]
gSN.datasets.append(batchjob.dataSet(['SNonly'], [SNdata,'background.ini',noCMB]))
gBAO.extra_opts ={'MPI_Converge_Stop': 0.002,'accuracy_level': 1.35,'high_accuracy_default':False}
groups.append(gSN)

# ranges for parameters when they are varied (can delete params if you just want to use defaults)
params = dict()
params['w'] = '-0.99 -3. 1 0.02 0.04'
params['wa'] = '0 -3 0.75 0.05 0.08'
params['mnu'] = '0.06 0.00 5 0.1 0.03'
params['omegak'] = '-0.0008 -0.8 0.8 0.001 0.003'  # starting exactly on flat seems to confuse minimizer


# extra parameters that are set only when specific parameters are varied. Can deleted to get defaults.
param_extra_opts = {
    'mnu': {'num_massive_neutrinos': 3,'neutrino_hierarchy': 'degenerate','MPI_converge_Stop':0.005}}



# try to match each new run name to exisitng covmat (covariance matrix for efficient exploration)
# e.g. try to map to get name without particular data combinations
covWithoutNameOrder = ['lensing', 'lowl','lowE','RSD']
# or try replacing various names (these are standard for the provided planck_covmats)
covNameMappings = {'BAORSD':SDSS,'SN':'Pantheon18','CMB':'Planck','DESlens':'WL','BAOnoLya':'BAO'}
# for mapping to nominal mission names try
# covNameMappings = {'plikHM':'planck','TT':'', 'lowTEB':'lowLike'}

covrenames = []
covrenames.append(['RSD','RSDonly'])
covrenames.append(['base_mnu_w_CMBLens_BAO','base_mnu_w_SDSS_Planck'])
covrenames.append(['base_mnu_w_CMBLens_BAO_SN','base_mnu_w_SDSS_Planck_SN'])
covrenames.append(['CMB_WL','Planck_DESlens'])
covrenames.append(['CMB_allGrowth','Planck_RSD_lensing_DESlens'])
covrenames.append(['CMBLens_SN','Planck_SN'])
covrenames.append(['CMBLens_BAO_SN','Planck_BAO_SN'])
covrenames.append(['CMBLens_BAO','Planck_BAO'])
covrenames.append(['CMB_SN','Planck_SN'])
covrenames.append(['CMB_BAO_SN','Planck_BAO_SN'])
covrenames.append(['CMB_BAOnoLya_SN','Planck_BAO_SN'])
covrenames.append(['CMBLens_BAORSD_DES','SDSS_Planck_DES'])
covrenames.append(['CMBLens_BAORSD','SDSS_Planck'])
covrenames.append(['CMBLens_BAORSD_DES_SN', 'Planck_SDSS_DES_lensing_SN'])
covrenames.append(['CMBLens_JLA_DR7','Planck_SDSS_DES_lensing_SN'])
covrenames.append(['WMAP_SN_DR7','WMAP_JLA_DR7'])
covrenames.append(['WMAP_JLA_DR7_DES','WMAP_JLA_DR7'])
covrenames.append(['CMBLens_DES_SN', 'Planck_DES_lensing_SN'])
covrenames.append(['WMAP_SN_DR7_BAOLya','WMAP_JLA_DR7'])
