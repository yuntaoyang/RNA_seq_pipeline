#!/usr/bin/env python
# coding: utf-8

# # run RNA_seq pipeline using python script (Salmon for single-end data)

# set up path

# In[1]:


# path to fasta file (raw data)
path_fasta = '/data/yyang18/test/'
# file_paired (filename and basename)
file_paired = '/home/yyang18/Project/SBMI_Zheng_NGS_Python_Script/RNA_seq/file_paired.csv'
# path to the output of trim_galore
path_trim_galore = '/home/yyang18/Project/SBMI_Zheng_NGS_Python_Script/RNA_seq/trim_galore/'
# path to the output of salmon
path_salmon = '/home/yyang18/Project/SBMI_Zheng_NGS_Python_Script/RNA_seq/salmon/'
# path to salmon index
path_index = '/data/yyang18/ref_geonome/salmon_index_human/salmon_index/'


# set up parameters

# In[2]:


# number of threads for salmon
thread = 32


# In[3]:


import subprocess
import os
import logging
import pandas as pd


# step1: read file_single

# In[4]:


file = pd.read_csv(file_paired)


# step2: create the output directory of trim_galore

# In[5]:


try:
    os.mkdir(path_trim_galore)
    logging.basicConfig(level=logging.DEBUG, 
                        filename="logfile", 
                        filemode="a",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logger = logging.getLogger(__name__)
    logger.info("create the output directory of trim_galore: "+path_trim_galore)
except:
    logging.basicConfig(level=logging.DEBUG, 
                        filename="logfile", 
                        filemode="a",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logger = logging.getLogger(__name__)
    logging.info("File exists: "+path_trim_galore)


# step3: create the output directory of salmon

# In[6]:


try:
    os.mkdir(path_salmon)
    logging.basicConfig(level=logging.DEBUG, 
                        filename="logfile", 
                        filemode="a",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logger = logging.getLogger(__name__)
    logger.info("create the output directory of salmon: "+path_salmon)
except:
    logging.basicConfig(level=logging.DEBUG, 
                        filename="logfile", 
                        filemode="a",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logger = logging.getLogger(__name__)
    logging.info("File exists: "+path_salmon)


# step4: trim_galore

# In[7]:


n = 0
for filename_1,filename_2,basename in zip(file['filename_1'],file['filename_2'],file['basename']):
    n = n + 1
    trim_galore = 'trim_galore -q 20'+' '                  '--phred33'+' '                  '--fastqc'+' '                  '--gzip'+' '                  '--length 36'+' '                  '--trim-n'+' '                  '--paired'+' '                  '-o'+' '+path_trim_galore+' '                  '--basename'+' '+basename+' '+                  path_fasta+filename_1+' '+                  path_fasta+filename_2
    if n == len(file['filename_1']):
        process = subprocess.Popen(trim_galore.split(), stdout=subprocess.PIPE)
        process.communicate()
    else:
        subprocess.Popen(trim_galore.split(), stdout=subprocess.PIPE)
logging.basicConfig(level=logging.DEBUG, 
                        filename="logfile", 
                        filemode="a",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
logger = logging.getLogger(__name__)
logging.info('trim_galore is done!')


# step5: salmon

# In[10]:


n = 0
for basename in file['basename']:
    n = n + 1
    salmon = 'salmon quant -i'+' '+path_index+' '+             '-p'+' '+str(thread)+' '+             '-l A'+' '+             '-1'+' '+path_trim_galore+basename+'_val_1.fq.gz'+' '+             '-2'+' '+path_trim_galore+basename+'_val_2.fq.gz'+' '+             '--validateMappings'+' '+             '-o'+' '+path_salmon+basename
    if n == len(file['basename']):
        process = subprocess.Popen(salmon.split(), stdout=subprocess.PIPE)
        process.communicate()
    else:
        subprocess.Popen(salmon.split(), stdout=subprocess.PIPE)
logging.basicConfig(level=logging.DEBUG, 
                        filename="logfile", 
                        filemode="a",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
logger = logging.getLogger(__name__)
logging.info('salmon is done!')    


# In[ ]:





# In[ ]:



