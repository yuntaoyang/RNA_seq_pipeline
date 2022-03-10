#!/usr/bin/env python
# coding: utf-8

# # run RNA_seq pipeline using python script (star for single-end data)

# set up path

# In[1]:


# path to fasta file (raw data)
path_fasta = '/data/yyang18/RNA_seq_single/'
# file_single (filename and basename)
file_single = '/home/yyang18/Project/SBMI_Zheng_NGS_Python_Script/RNA_seq/file_single.csv'
# path to the output of trim_galore
path_trim_galore = '/data/yyang18/trim_galore/'
# path to the output of star
path_star = '/data/yyang18/star/'
# path to star index
path_index = '/data/yyang18/ref_geonome/star_index_human/index/'


# set up parameters

# In[2]:


# number of threads for star
thread = 8
# number of threads for bam sorting
thread_sort = 4
# name of logfile
log = 'logfile_star'


# In[3]:


import subprocess
import os
import logging
import pandas as pd


# step1: read file_single

# In[4]:


file = pd.read_csv(file_single)


# In[5]:


logging.basicConfig(level=logging.DEBUG, 
                        filename=log, 
                        filemode="a",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
logger = logging.getLogger(__name__)


# step2: create the output directory of trim_galore

# In[6]:


if os.path.isdir(path_trim_galore):
    logger.info("File exists: "+path_trim_galore)
else:
    os.mkdir(path_trim_galore)
    logger.info("create the output directory of trim_galore: "+path_trim_galore)


# step3: create the output directory of star

# In[7]:


if os.path.isdir(path_star):
    logger.info("File exists: "+path_star)
else:
    os.mkdir(path_star)
    logger.info("create the output directory of star: "+path_star)


# step4: trim_galore

# In[8]:


for filename,basename in zip(file['filename'],file['basename']):
    trim_galore = 'trim_galore -q 20'+' '                  '--phred33'+' '                  '--fastqc'+' '                  '--gzip'+' '                  '--length 36'+' '                  '--trim-n'+' '                  '-o'+' '+path_trim_galore+' '                  '--basename'+' '+basename+' '+                  path_fasta+filename+' '+                  '>'+' '+path_trim_galore+basename+'_trim.log'+' '+'2>&1'
    process = subprocess.Popen(trim_galore,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    process.communicate()
    logger.info(basename+' trim_galore is done!')


# step5: star

# In[9]:


for basename in file['basename']:
    star = 'STAR --runThreadN'+' '+str(thread)+' '+            '--outBAMsortingThreadN'+' '+str(thread_sort)+' '+            '--genomeDir'+' '+path_index+' '+            '--readFilesIn'+' '+path_trim_galore+basename+'_trimmed.fq.gz'+' '+            '--readFilesCommand gunzip -c'+' '+            '--outSAMtype BAM SortedByCoordinate'+' '+            '--twopassMode Basic'+' '+            '--quantMode GeneCounts'+' '+            '--outFileNamePrefix'+' '+path_star+basename+' '+            '>'+' '+path_star+basename+'_star.log'+' '+'2>&1'
    process = subprocess.Popen(star,shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    process.communicate()
    logger.info(basename+' star is done!')     


# In[ ]:




