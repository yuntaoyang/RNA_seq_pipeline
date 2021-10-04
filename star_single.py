#!/usr/bin/env python
# coding: utf-8

# # run RNA_seq pipeline using python script (star for single-end data)

# set up path

# In[4]:


# path to fasta file (raw data)
path_fasta = '/data/yyang18/fastq/'
# file_single (filename and basename)
file_single = '/home/yyang18/Project/SBMI_Zheng_NGS_Python_Script/RNA_seq/file_single.csv'
# path to the output of trim_galore
path_trim_galore = '/home/yyang18/Project/SBMI_Zheng_NGS_Python_Script/RNA_seq/trim_galore/'
# path to the output of star
path_star = '/home/yyang18/Project/SBMI_Zheng_NGS_Python_Script/RNA_seq/star/'
# path to star index
path_index = '/data/yyang18/ref_geonome/star_index_human/index/'


# set up parameters

# In[5]:


# number of threads for star
thread = 64
# number of threads for bam sorting
thread_sort = 4


# In[6]:


import subprocess
import os
import logging
import pandas as pd


# step1: read file_single

# In[7]:


file = pd.read_csv(file_single)


# step2: create the output directory of trim_galore

# In[8]:


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


# step3: create the output directory of star

# In[9]:


try:
    os.mkdir(path_star)
    logging.basicConfig(level=logging.DEBUG, 
                        filename="logfile", 
                        filemode="a",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logger = logging.getLogger(__name__)
    logger.info("create the output directory of star: "+path_star)
except:
    logging.basicConfig(level=logging.DEBUG, 
                        filename="logfile", 
                        filemode="a",
                        format="%(asctime)-15s %(levelname)-8s %(message)s")
    logger = logging.getLogger(__name__)
    logging.info("File exists: "+path_star)


# step4: trim_galore

# In[11]:


n = 0
for filename,basename in zip(file['filename'],file['basename']):
    n = n + 1
    trim_galore = 'trim_galore -q 20'+' '+                  '--phred33'+' '+\ 
                  '--fastqc'+' '+                  '--gzip'+' '+\ 
                  '--length 36'+' '+\ 
                  '--trim-n'+' '+                  '-o'+' '+path_trim_galore+' '+                  '--basename'+' '+basename+' '+                  path_fasta+filename
    if n == len(file['filename']):
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


# step5: star

# In[12]:


for basename in file['basename']:
    star = 'STAR --runThreadN'+' '+str(thread)+' '+            '--outBAMsortingThreadN'+' '+str(thread_sort)+' '+            '--genomeDir'+' '+path_index+' '+            '--readFilesIn'+' '+path_trim_galore+basename+'_trimmed.fq.gz'+' '+            '--readFilesCommand gunzip -c'+' '+            '--outSAMtype BAM SortedByCoordinate'+' '+            '--twopassMode Basic'+' '+            '--quantMode GeneCounts'+' '+            '--outFileNamePrefix'+' '+path_star+basename
    process = subprocess.Popen(star.split(), stdout=subprocess.PIPE)
    process.communicate()
logging.basicConfig(level=logging.DEBUG, 
                    filename="logfile", 
                    filemode="a",
                    format="%(asctime)-15s %(levelname)-8s %(message)s")
logger = logging.getLogger(__name__)
logging.info('star is done!')    


# In[ ]:




