"""
This module is to set common variables
Left to update.
"""
import os 
class config(object):
    def __init__(self, root_path, log_dir="", data_type="f"):
        # path
        self.log_dir = log_dir
        self.root_path = root_path
        self.data_path = os.path.join(self.root_path, "data/" + data_type + "-data/")
        self.assay_data_path = os.path.join(self.data_path, "assay-data/bedGraph/")
        self.pilot_region_path = os.path.join(self.data_path, "pilot-data")
        self.pilot_region_file = os.path.join(self.data_path, "pilot-data/encodePilotRegions.hg19.bed")
        self.processed_pilot_data_path = os.path.join(self.data_path, "processed-pilot-data/")
        self.blacklist_region_file = os.path.join(self.data_path, "blacklist-data/ENCFF419RSJ.bed")
        self.blacklist_rm_data_path = os.path.join(self.data_path, "blacklist-rm-data/")
        self.experiment_path = os.path.join(self.root_path, "experiment-2")
        self.experiment1_path = os.path.join(self.root_path, "experiment-1")
        self.assay_region_file = os.path.join(self.data_path, "assay-data/assay-region.bed")
    
        # others
        self.chromosome_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 
                                'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                                'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 
                                'chr19', 'chr20', 'chr21', 'chr22','chrX']
        self.chromosome_list_test = ['chr21', 'chr22']
        self.roadmap_assay_names_file = os.path.join(self.data_path, "assay-data", "roadmap_assays.txt")
        self.e003_assays = ['E003-DNase', 
                'E003-H2A.Z', 'E003-H2AK5ac', 
                'E003-H2BK5ac', 'E003-H2BK12ac', 'E003-H2BK15ac', 'E003-H2BK20ac', 'E003-H2BK120ac', 
                'E003-H3K4ac', 'E003-H3K4me1', 'E003-H3K4me2', 'E003-H3K4me3', 
                'E003-H3K9ac', 'E003-H3K9me3', 'E003-H3K14ac', 'E003-H3K18ac', 'E003-H3K23ac', 
                'E003-H3K23me2', 'E003-H3K27ac', 'E003-H3K27me3', 'E003-H3K36me3', 'E003-H3K56ac', 
                'E003-H3K79me1', 'E003-H3K79me2', 'E003-H4K5ac', 'E003-H4K8ac', 'E003-H4K20me1', 
                'E003-H4K91ac']
        self.e017_assays = ['E017-DNase','E017-H2A.Z','E017-H2AK5ac','E017-H2AK9ac','E017-H2BK5ac','E017-H2BK12ac','E017-H2BK15ac','E017-H2BK20ac',
                'E017-H2BK120ac','E017-H3K4ac','E017-H3K4me1','E017-H3K4me2','E017-H3K4me3','E017-H3K9ac','E017-H3K9me1','E017-H3K9me3',
                'E017-H3K14ac','E017-H3K18ac','E017-H3K23ac','E017-H3K27ac','E017-H3K27me3','E017-H3K36me3','E017-H3K56ac','E017-H3K79me1',
                'E017-H3K79me2','E017-H4K5ac','E017-H4K8ac','E017-H4K20me1','E017-H4K91ac']
        self.e114_assays = ['E114-H4K20me1','E114-H3K79me2','E114-H3K36me3','E114-H3K27me3','E114-H3K27ac','E114-H3K9me3','E114-H3K9ac','E114-H3K4me3',
                'E114-H3K4me2','E114-H3K4me1','E114-H2A.Z','E114-DNase']
        self.e116_assays = ['E116-H4K20me1','E116-H3K79me2','E116-H3K36me3','E116-H3K27me3','E116-H3K27ac','E116-H3K9me3','E116-H3K9ac','E116-H3K4me3',
                'E116-H3K4me2','E116-H3K4me1','E116-H2A.Z','E116-DNase']
        self.e122_assays = ['E122-H4K20me1','E122-H3K79me2','E122-H3K36me3','E122-H3K27me3','E122-H3K27ac','E122-H3K9me3','E122-H3K9me1','E122-H3K9ac',
                'E122-H3K4me3','E122-H3K4me2','E122-H3K4me1','E122-H2A.Z','E122-DNase']
        self.e123_assays = ['E123-DNase','E123-H2A.Z','E123-H3K4me1','E123-H3K4me2','E123-H3K4me3','E123-H3K9ac','E123-H3K9me1','E123-H3K9me3',
                'E123-H3K27ac','E123-H3K27me3','E123-H3K36me3','E123-H3K79me2','E123-H4K20me1']
        self.e127_assays = ['E127-DNase','E127-H2A.Z','E127-H3K4me1','E127-H3K4me2','E127-H3K4me3','E127-H3K9ac','E127-H3K9me1','E127-H3K9me3',
                'E127-H3K27ac','E127-H3K27me3','E127-H3K36me3','E127-H3K79me2','E127-H4K20me1']
        self.e128_assays = ['E128-DNase','E128-H2A.Z','E128-H3K4me1','E128-H3K4me2','E128-H3K4me3','E128-H3K9ac','E128-H3K9me3','E128-H3K27ac',
                'E128-H3K27me3','E128-H3K36me3','E128-H3K79me2','E128-H4K20me1']
        self.e116_chromHMM_assays = ['E116-H3K4me3', 'E116-H3K4me1', 'E116-H3K36me3', 'E116-H3K27me3', 'E116-H3K9me3']
        self.e003_chromHMM_assays = ['E003-H3K4me3', 'E003-H3K4me1', 'E003-H3K36me3', 'E003-H3K27me3', 'E003-H3K9me3']
        self.gexp_signal_file = self.data_path + "tss-data/Ensembl_v65.Gencode_v10.ENSG.gene_info"
        self.gexp_region_file = self.data_path + "tss-data/57epigenomes.RPKM.pc"
        self.gexp_file = self.data_path + "tss-data/E116-gexp.bed"
        self.annotation_result_dir = self.data_path + "annotation-result"
        self.region_result_dir = self.data_path + "region-result"
        self.enhancer_file = self.data_path + "enhancer-data/enhancer.bed"
        self.filtered_gexp_file = os.path.join(self.experiment1_path, "35-0728-18-gexpEnhancerFilter-p", "exp1", "filtered_gexp.bed")
        self.filtered_enhancer_file = os.path.join(self.experiment1_path, "35-0728-18-gexpEnhancerFilter-p", "exp2", "filtered_enhancer.bed")
        self.genome_browser_dir = self.data_path + "genome-browser"
        self.chromhmm_data_dir = os.path.join(self.data_path, 'assay-data', 'chromhmm')
        self.segway_data_dir = os.path.join(self.data_path, 'assay-data', 'segway')
        self.models_dir = os.path.join(self.root_path, 'src', 'models')

        self.exp_ids = [
            'exp6-1', 'exp6-2', 'exp6-3', 'exp6-4', 'exp6-5', 'exp6-6',
            'exp8-1', 'exp8-2', 'exp8-3', 'exp8-4', 'exp8-5', 'exp8-6',
            'exp9-1-con', 'exp9-2-con', 'exp9-3-con', 'exp9-4-con', 'exp9-5-con', 'exp9-6-con',
            'exp9-1-dis', 'exp9-2-dis', 'exp9-3-dis', 'exp9-4-dis', 'exp9-5-dis', 'exp9-6-dis',
            'exp11-1-con', 'exp11-2-con', 'exp11-3-con', 'exp11-4-con', 'exp11-5-con', 'exp11-6-con',
            'exp11-1-dis', 'exp11-2-dis', 'exp11-3-dis', 'exp11-4-dis', 'exp11-5-dis', 'exp11-6-dis',
            'exp4-1_1', 'exp4-1_2', 'exp4-1_3', 'exp4-1_4', 'exp4-1_5', 'exp4-1_6',
            'exp0-1', 'exp0-2', 'exp0-3', 'exp0-4', 'exp0-5', 'exp0-6',
            'exp17-1', 'exp17-2', 'exp17-3', 'exp17-4', 'exp17-5', 'exp17-6', 'exp17-7', 'exp17-8', 'exp17-9', 'exp17-10', 'exp17-11', 'exp17-12',
            'exp47-1', 'exp47-2', 'exp47-3', 'exp47-4', 'exp47-5', 'exp47-6',
            'exp48-1', 'exp48-2', 'exp48-3', 'exp48-4', 'exp48-5', 'exp48-6',
            'exp49-1', 'exp49-2', 'exp49-3', 'exp49-4', 'exp49-5', 'exp49-6',
            'exp50-1', 'exp50-2', 'exp50-3', 'exp50-4', 'exp50-5', 'exp50-6',
            'exp51-1', 'exp51-2', 'exp51-3', 'exp51-4', 'exp51-5', 'exp51-6',
            'exp52-1', 'exp52-2', 'exp52-3', 'exp52-4', 'exp52-5', 'exp52-6',
            'exp53-1', 'exp53-2', 'exp53-3', 'exp53-4', 'exp53-5', 'exp53-6',
            'exp54-1', 'exp54-2', 'exp54-3', 'exp54-4', 'exp54-5', 'exp54-6',
            'exp55-1', 'exp55-2', 'exp55-3', 'exp55-4', 'exp55-5', 'exp55-6',
            'exp56-1', 'exp56-2', 'exp56-3', 'exp56-4', 'exp56-5', 'exp56-6',
            'exp57-1', 'exp57-2', 'exp57-3', 'exp57-4', 'exp57-5', 'exp57-6',
            'exp58-1', 'exp58-2', 'exp58-3', 'exp58-4', 'exp58-5', 'exp58-6',
            'exp59-1', 'exp59-2', 'exp59-3', 'exp59-4', 'exp59-5', 'exp59-6',
            'exp60-1', 'exp60-2', 'exp60-3', 'exp60-4', 'exp60-5', 'exp60-6',
            'exp61-1', 'exp61-2', 'exp61-3', 'exp61-4', 'exp61-5', 'exp61-6',
            'exp62-1', 'exp62-2', 'exp62-3', 'exp62-4', 'exp62-5', 'exp62-6',
            'exp63-1', 'exp63-2', 'exp63-3', 'exp63-4', 'exp63-5', 'exp63-6',
            'exp64-1', 'exp64-2', 'exp64-3', 'exp64-4', 'exp64-5', 'exp64-6',
            'exp65-1', 'exp65-2', 'exp65-3', 'exp65-4', 'exp65-5', 'exp65-6',
            'exp66-1', 'exp66-2', 'exp66-3', 'exp66-4', 'exp66-5', 'exp66-6',
            'exp68-1', 'exp68-2', 'exp68-3', 'exp68-4', 'exp68-5', 'exp68-6',
            'exp68-7', 'exp68-8', 'exp68-9', 'exp68-10', 'exp68-11', 'exp68-12',
            'exp85-1-con', 'exp85-2-con', 'exp85-3-con', 'exp85-4-con', 'exp85-5-con', 'exp85-6-con',
            'exp85-1-dis', 'exp85-2-dis', 'exp85-3-dis', 'exp85-4-dis', 'exp85-5-dis', 'exp85-6-dis',
            'exp86-1-con', 'exp86-2-con', 'exp86-3-con', 'exp86-4-con', 'exp86-5-con', 'exp86-6-con',
            'exp86-1-dis', 'exp86-2-dis', 'exp86-3-dis', 'exp86-4-dis', 'exp86-5-dis', 'exp86-6-dis',
            'exp87-1', 'exp87-2', 'exp87-3', 'exp87-4', 'exp87-5', 'exp87-6',
            'exp88-1', 'exp88-2', 'exp88-3',
            'exp89-1', 'exp89-2', 'exp89-3',
            'exp90-1', 'exp90-2', 'exp90-3', 'exp90-4', 'exp90-5', 'exp90-6',
            'exp91-1', 'exp91-2', 'exp91-3', 'exp91-4', 'exp91-5', 'exp91-6']

        self.state_dic = {
            'exp6-1':5, 'exp6-2':3, 'exp6-3':8, 'exp6-4':10, 'exp6-5':12, 'exp6-6':15,
            'exp8-1':5, 'exp8-2':3, 'exp8-3':8, 'exp8-4':10, 'exp8-5':12, 'exp8-6':15,
            'exp9-1-con':3, 'exp9-2-con':5, 'exp9-3-con':8, 'exp9-4-con':10, 'exp9-5-con':12, 'exp9-6-con':15,
            'exp9-1-dis':3, 'exp9-2-dis':5, 'exp9-3-dis':8, 'exp9-4-dis':10, 'exp9-5-dis':12, 'exp9-6-dis':15,
            'exp11-1-con':3, 'exp11-2-con':5, 'exp11-3-con':8, 'exp11-4-con':10, 'exp11-5-con':12, 'exp11-6-con':15,
            'exp11-1-dis':3, 'exp11-2-dis':5, 'exp11-3-dis':8, 'exp11-4-dis':10, 'exp11-5-dis':12, 'exp11-6-dis':15,
            'exp4-1_1':3, 'exp4-1_2':5, 'exp4-1_3':8, 'exp4-1_4':10, 'exp4-1_5':12, 'exp4-1_6':15,
            'exp17-1':3, 'exp17-2':5, 'exp17-3':8, 'exp17-4':10, 'exp17-5':12, 'exp17-6':15, 'exp17-7':3, 'exp17-8':5, 'exp17-9':8, 'exp17-10':10, 'exp17-11':12, 'exp17-12':15,
            'exp47-1':3, 'exp47-2':5, 'exp47-3':8, 'exp47-4':10, 'exp47-5':12, 'exp47-6':15,
            'exp48-1':3, 'exp48-2':5, 'exp48-3':8, 'exp48-4':10, 'exp48-5':12, 'exp48-6':15,
            'exp49-1':3, 'exp49-2':5, 'exp49-3':8, 'exp49-4':10, 'exp49-5':12, 'exp49-6':15,
            'exp50-1':3, 'exp50-2':5, 'exp50-3':8, 'exp50-4':10, 'exp50-5':12, 'exp50-6':15,
            'exp51-1':3, 'exp51-2':5, 'exp51-3':8, 'exp51-4':10, 'exp51-5':12, 'exp51-6':15,
            'exp52-1':3, 'exp52-2':5, 'exp52-3':8, 'exp52-4':10, 'exp52-5':12, 'exp52-6':15,
            'exp53-1':3, 'exp53-2':5, 'exp53-3':8, 'exp53-4':10, 'exp53-5':12, 'exp53-6':15,
            'exp54-1':3, 'exp54-2':5, 'exp54-3':8, 'exp54-4':10, 'exp54-5':12, 'exp54-6':15,
            'exp55-1':3, 'exp55-2':5, 'exp55-3':8, 'exp55-4':10, 'exp55-5':12, 'exp55-6':15,
            'exp56-1':3, 'exp56-2':5, 'exp56-3':8, 'exp56-4':10, 'exp56-5':12, 'exp56-6':15,
            'exp57-1':3, 'exp57-2':5, 'exp57-3':8, 'exp57-4':10, 'exp57-5':12, 'exp57-6':15,
            'exp58-1':3, 'exp58-2':5, 'exp58-3':8, 'exp58-4':10, 'exp58-5':12, 'exp58-6':15,
            'exp59-1':3, 'exp59-2':5, 'exp59-3':8, 'exp59-4':10, 'exp59-5':12, 'exp59-6':15,
            'exp60-1':3, 'exp60-2':5, 'exp60-3':8, 'exp60-4':10, 'exp60-5':12, 'exp60-6':15,
            'exp61-1':3, 'exp61-2':5, 'exp61-3':8, 'exp61-4':10, 'exp61-5':12, 'exp61-6':15,
            'exp62-1':3, 'exp62-2':5, 'exp62-3':8, 'exp62-4':10, 'exp62-5':12, 'exp62-6':15,
            'exp63-1':3, 'exp63-2':5, 'exp63-3':8, 'exp63-4':10, 'exp63-5':12, 'exp63-6':15,
            'exp64-1':3, 'exp64-2':5, 'exp64-3':8, 'exp64-4':10, 'exp64-5':12, 'exp64-6':15,
            'exp65-1':3, 'exp65-2':5, 'exp65-3':8, 'exp65-4':10, 'exp65-5':12, 'exp65-6':15,
            'exp66-1':3, 'exp66-2':5, 'exp66-3':8, 'exp66-4':10, 'exp66-5':12, 'exp66-6':15,
            'exp68-1':3, 'exp68-2':5, 'exp68-3':8, 'exp68-4':10, 'exp68-5':12, 'exp68-6':15,
            'exp68-7':3, 'exp68-8':5, 'exp68-9':8, 'exp68-10':10, 'exp68-11':12, 'exp68-12':15,
            'exp85-1-con':3, 'exp85-2-con':5, 'exp85-3-con':8, 'exp85-4-con':10, 'exp85-5-con':12, 'exp85-6-con':15,
            'exp85-1-dis':3, 'exp85-2-dis':5, 'exp85-3-dis':8, 'exp85-4-dis':10, 'exp85-5-dis':12, 'exp85-6-dis':15,
            'exp86-1-con':3, 'exp86-2-con':5, 'exp86-3-con':8, 'exp86-4-con':10, 'exp86-5-con':12, 'exp86-6-con':15,
            'exp86-1-dis':3, 'exp86-2-dis':5, 'exp86-3-dis':8, 'exp86-4-dis':10, 'exp86-5-dis':12, 'exp86-6-dis':15,
            'exp87-1':3, 'exp87-2':5, 'exp87-3':8, 'exp87-4':10, 'exp87-5':12, 'exp87-6':15,
            'exp88-1':1, 'exp88-2':3, 'exp88-3':5,
            'exp89-1':1, 'exp89-2':3, 'exp89-3':5,
            'exp90-1-dis':3, 'exp90-2-dis':5, 'exp90-3-dis':8, 'exp90-4-dis':10, 'exp90-5-dis':12, 'exp90-6-dis':15,
            'exp90-1-con':3, 'exp90-2-con':5, 'exp90-3-con':8, 'exp90-4-con':10, 'exp90-5-con':12, 'exp90-6-con':15,
            'exp91-1':3, 'exp91-2':5, 'exp91-3':8, 'exp91-4':10, 'exp91-5':12, 'exp91-6':15
            }

        self.name_dic = {
            'exp6-1':'NMF', 'exp6-2':'NMF', 'exp6-3':'NMF', 'exp6-4':'NMF', 'exp6-5':'NMF', 'exp6-6':'NMF',
            'exp8-1':'PCA', 'exp8-2':'PCA', 'exp8-3':'PCA', 'exp8-4':'PCA', 'exp8-5':'PCA', 'exp8-6':'PCA',
            'exp9-1-con':'HMMber-con', 'exp9-2-con':'HMMber-con', 'exp9-3-con':'HMMber-con', 'exp9-4-con':'HMMber-con', 'exp9-5-con':'HMMber-con', 'exp9-6-con':'HMMber-con',
            'exp9-1-dis':'HMMber-dis', 'exp9-2-dis':'HMMber-dis', 'exp9-3-dis':'HMMber-dis', 'exp9-4-dis':'HMMber-dis', 'exp9-5-dis':'HMMber-dis', 'exp9-6-dis':'HMMber-dis',
            'exp11-1-con':'HMMgaus-con', 'exp11-2-con':'HMMgaus-con', 'exp11-3-con':'HMMgaus-con', 'exp11-4-con':'HMMgaus-con', 'exp11-5-con':'HMMgaus-con', 'exp11-6-con':'HMMgaus-con',
            'exp11-1-dis':'HMMgaus-dis', 'exp11-2-dis':'HMMgaus-dis', 'exp11-3-dis':'HMMgaus-dis', 'exp11-4-dis':'HMMgaus-dis', 'exp11-5-dis':'HMMgaus-dis', 'exp11-6-dis':'HMMgaus-dis',
            'exp4-1_1':'Raw', 'exp4-1_2':'Raw', 'exp4-1_3':'Raw', 'exp4-1_4':'Raw', 'exp4-1_5':'Raw', 'exp4-1_6':'Raw',
            'exp17-1':'SSM', 'exp17-2':'SSM', 'exp17-3':'SSM', 'exp17-4':'SSM', 'exp17-5':'SSM', 'exp17-6':'SSM', 'exp17-7':'SSMnn', 'exp17-8':'SSMnn', 'exp17-9':'SSMnn', 'exp17-10':'SSMnn', 'exp17-11':'SSMnn', 'exp17-12':'SSMnn',
            'exp47-1':'SSM', 'exp47-2':'SSM', 'exp47-3':'SSM', 'exp47-4':'SSM', 'exp47-5':'SSM', 'exp47-6':'SSM',
            'exp48-1':'SSMnn', 'exp48-2':'SSMnn', 'exp48-3':'SSMnn', 'exp48-4':'SSMnn', 'exp48-5':'SSMnn', 'exp48-6':'SSMnn',
            'exp49-1':'SSM', 'exp49-2':'SSM', 'exp49-3':'SSM', 'exp49-4':'SSM', 'exp49-5':'SSM', 'exp49-6':'SSM',
            'exp50-1':'SSMnn', 'exp50-2':'SSMnn', 'exp50-3':'SSMnn', 'exp50-4':'SSMnn', 'exp50-5':'SSMnn', 'exp50-6':'SSMnn',
            'exp51-1':'SSM', 'exp51-2':'SSM', 'exp51-3':'SSM', 'exp51-4':'SSM', 'exp51-5':'SSM', 'exp51-6':'SSM',
            'exp52-1':'SSMnn', 'exp52-2':'SSMnn', 'exp52-3':'SSMnn', 'exp52-4':'SSMnn', 'exp52-5':'SSMnn', 'exp52-6':'SSMnn',
            'exp53-1':'SSM', 'exp53-2':'SSM', 'exp53-3':'SSM', 'exp53-4':'SSM', 'exp53-5':'SSM', 'exp53-6':'SSM',
            'exp54-1':'SSMnn', 'exp54-2':'SSMnn', 'exp54-3':'SSMnn', 'exp54-4':'SSMnn', 'exp54-5':'SSMnn', 'exp54-6':'SSMnn',
            'exp55-1':'SSM', 'exp55-2':'SSM', 'exp55-3':'SSM', 'exp55-4':'SSM', 'exp55-5':'SSM', 'exp55-6':'SSM',
            'exp56-1':'SSMnn', 'exp56-2':'SSMnn', 'exp56-3':'SSMnn', 'exp56-4':'SSMnn', 'exp56-5':'SSMnn', 'exp56-6':'SSMnn',
            'exp57-1':'SSM', 'exp57-2':'SSM', 'exp57-3':'SSM', 'exp57-4':'SSM', 'exp57-5':'SSM', 'exp57-6':'SSM',
            'exp58-1':'SSMnn', 'exp58-2':'SSMnn', 'exp58-3':'SSMnn', 'exp58-4':'SSMnn', 'exp58-5':'SSMnn', 'exp58-6':'SSMnn',
            'exp59-1':'SSM', 'exp59-2':'SSM', 'exp59-3':'SSM', 'exp59-4':'SSM', 'exp59-5':'SSM', 'exp59-6':'SSM',
            'exp60-1':'SSMnn', 'exp60-2':'SSMnn', 'exp60-3':'SSMnn', 'exp60-4':'SSMnn', 'exp60-5':'SSMnn', 'exp60-6':'SSMnn',
            'exp61-1':'SSM', 'exp61-2':'SSM', 'exp61-3':'SSM', 'exp61-4':'SSM', 'exp61-5':'SSM', 'exp61-6':'SSM',
            'exp62-1':'SSMnn', 'exp62-2':'SSMnn', 'exp62-3':'SSMnn', 'exp62-4':'SSMnn', 'exp62-5':'SSMnn', 'exp62-6':'SSMnn',
            'exp63-1':'SSM', 'exp63-2':'SSM', 'exp63-3':'SSM', 'exp63-4':'SSM', 'exp63-5':'SSM', 'exp63-6':'SSM',
            'exp64-1':'SSMnn', 'exp64-2':'SSMnn', 'exp64-3':'SSMnn', 'exp64-4':'SSMnn', 'exp64-5':'SSMnn', 'exp64-6':'SSMnn',
            'exp65-1':'SSM', 'exp65-2':'SSM', 'exp65-3':'SSM', 'exp65-4':'SSM', 'exp65-5':'SSM', 'exp65-6':'SSM',
            'exp66-1':'SSMnn', 'exp66-2':'SSMnn', 'exp66-3':'SSMnn', 'exp66-4':'SSMnn', 'exp66-5':'SSMnn', 'exp66-6':'SSMnn',
            'exp68-1':'epigenome-ssm', 'exp68-2':'epigenome-ssm', 'exp68-3':'epigenome-ssm', 'exp68-4':'epigenome-ssm', 'exp68-5':'epigenome-ssm', 'exp68-6':'epigenome-ssm',
            'exp68-7':'epigenome-ssm(non-neg)', 'exp68-8':'epigenome-ssm(non-neg)', 'exp68-9':'epigenome-ssm(non-neg)', 'exp68-10':'epigenome-ssm(non-neg)', 'exp68-11':'epigenome-ssm(non-neg)', 'exp68-12':'epigenome-ssm(non-neg)',
            'exp85-1-con':'HMMber-con', 'exp85-2-con':'HMMber-con', 'exp85-3-con':'HMMber-con', 'exp85-4-con':'HMMber-con', 'exp85-5-con':'HMMber-con', 'exp85-6-con':'HMMber-con',
            'exp85-1-dis':'HMMber-dis', 'exp85-2-dis':'HMMber-dis', 'exp85-3-dis':'HMMber-dis', 'exp85-4-dis':'HMMber-dis', 'exp85-5-dis':'HMMber-dis', 'exp85-6-dis':'HMMber-dis',
            'exp86-1-con':'HMMgaus-con', 'exp86-2-con':'HMMgaus-con', 'exp86-3-con':'HMMgaus-con', 'exp86-4-con':'HMMgaus-con', 'exp86-5-con':'HMMgaus-con', 'exp86-6-con':'HMMgaus-con',
            'exp86-1-dis':'HMMgaus-dis', 'exp86-2-dis':'HMMgaus-dis', 'exp86-3-dis':'HMMgaus-dis', 'exp86-4-dis':'HMMgaus-dis', 'exp86-5-dis':'HMMgaus-dis', 'exp86-6-dis':'HMMgaus-dis',
            'exp87-1':'NMF', 'exp87-2':'NMF', 'exp87-3':'NMF', 'exp87-4':'NMF', 'exp87-5':'NMF', 'exp87-6':'NMF',
            'exp88-1':'PCA', 'exp88-2':'PCA', 'exp88-3':'PCA',
            'exp89-1':'Raw', 'exp89-2':'Raw', 'exp89-3':'Raw',
            'exp90-1-dis':'chromHMM-dis', 'exp90-2-dis':'chromHMM-dis', 'exp90-3-dis':'chromHMM-dis', 'exp90-4-dis':'chromHMM-dis', 'exp90-5-dis':'chromHMM-dis', 'exp90-6-dis':'chromHMM-dis',
            'exp90-1-con':'chromHMM-con', 'exp90-2-con':'chromHMM-con', 'exp90-3-con':'chromHMM-con', 'exp90-4-con':'chromHMM-con', 'exp90-5-con':'chromHMM-con', 'exp90-6-con':'chromHMM-con',
            'exp91-1':'Segway', 'exp91-2':'Segway', 'exp91-3':'Segway', 'exp91-4':'Segway', 'exp91-5':'Segway', 'exp91-6':'Segway'
        }    
        
    
    def pilot_index_region(self, resolution_size):
        return os.path.join(self.pilot_region_path, "index_resol" + str(resolution_size) + "bp.bed")

    def openLog(self, log_name = "log"):
        self.log_file = os.path.join(self.log_dir, log_name + ".txt")
        self.log_f = open(self.log_file, 'w')

    def toLog(self, msg):
        # using print just for convenience, not a good implementation anyway
        print(msg, file=self.log_f)

    def closeLog(self):
        self.log_f.close()

