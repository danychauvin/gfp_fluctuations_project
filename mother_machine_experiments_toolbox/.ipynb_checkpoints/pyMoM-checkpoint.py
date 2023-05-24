import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
import shutil
import subprocess

def extractAndTidyCurated(experimentFolderPath,preprocFolder,subset,rx0,gls_dic,analysis_id="/"):
    
    proms=list(subset['promoter'])
    poss=list(subset['pos'])
    vecs=list(subset['vector'])
    unique_proms=list(set(proms))
    strainToPosDic={s:[poss[i] for i in range(0,len(proms)) if proms[i]==s ] for s in unique_proms}
    posToStrainDic={'Pos'+str(e):key for key in strainToPosDic.keys() for e in strainToPosDic[key]}

    for i in range(0,len(proms)):
        prom=proms[i]
        vec=vecs[i]
        pos=poss[i]
        for gl in gls_dic[pos]:
            for root, dirs, files in os.walk("%s/Pos%s/Pos%s_GL%i"%(experimentFolderPath+"/"+preprocFolder,pos,pos,gl)):
                #walking in the preprocFolder/Posi
                for fileName in files:
                    if re.search(rx0,fileName): # if a relevant CSV file is found, and it is not in a backed-up folder
                        fullPath=os.path.join(root,fileName) #full path to the csv file, should NOT contain `curated`
                        #print(fullPath)

                        if analysis_id!="/":
                            folder_name="/"+analysis_id+"/"
                        else:
                            folder_name=analysis_id
                        #print(folder_name)

                        if 'BKP' not in fullPath and folder_name in fullPath: #not a backup and found in the desired batch folder
                            #print(fullPath)
                            outputDir="%s/curated_data/%s_%s_curated"%(experimentFolderPath,prom,vec) #date_bottom|top_strain_curated
                            outputPath="%s/curated_data/%s_%s_curated/%s"%(experimentFolderPath,prom,vec,fileName) #date_bottom|top_strain_curated
                            #print(outputPath)
                            if "curated" not in fullPath: #avoid copying already curated files
                                try:
                                    os.mkdir(outputDir)
                                except:
                                    print("")
                                try: 
                                    shutil.copy(fullPath,outputPath)
                                    print("Copying %s to %s."%(fullPath,outputPath))
                                except:
                                    print("Error for %s: cannot copy curated folder. Curated CSV file certainly already exists."%(outputFolderPath))

def walklevel(some_dir, level=1):
    some_dir = some_dir.rstrip(os.path.sep)
    #try:
    assert os.path.isdir(some_dir)
    #except:
    #print("The path to preprocessing data does not exist. They are certainly wrong.")
    #return()
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]
            
            
def makeDataFrame(a_,i_):
    b_=a_[1]
    b_["peak_locations"]=a_[0]
    b_["index"]=np.array([i_]*len(a_[0]))
    final_data_frame = pd.DataFrame(data=b_)
    return(final_data_frame)
    
def submit_preprocessing(experimentFolderPath,preprocFolder,template,NORMALIZATION_REGION_OFFSET,pos_var,RAW_PATH,FLATFIELD_PATH,ROTATION,TMAX):    
    preproc_dir="%s/%s"%(experimentFolderPath,preprocFolder)
    PREPROC_DIR_TPL=preproc_dir+"/mmpreproc_%s/"
    GL_DETECTION_TEMPLATE_PATH="%s/%s/template_config.json"%(experimentFolderPath,template)
    IMAGE_REGISTRATION_METHOD="1"
    #FRAMES_TO_IGNORE="()"
    NORMALIZATION_CONFIG_PATH="true"
    
    decision=input("Do you really want to submit preprocessing tasks to the scheduler? If the preprocessing folder target path has not been changed, it might result in overwritting your data. Enter 'Yes', to confirm.")
    
    if decision=="Yes":
        
        if os.path.exists("./pos_names.sh"):
            os.remove("./pos_names.sh")

        f=open("./pos_names.sh","w")
        f.write(pos_var)
        f.close()

        subprocess.run(["%s/submit_slurm.sh"%(os.getcwd()),PREPROC_DIR_TPL, RAW_PATH, FLATFIELD_PATH,"pos_names.sh",GL_DETECTION_TEMPLATE_PATH,ROTATION,TMAX,IMAGE_REGISTRATION_METHOD,NORMALIZATION_CONFIG_PATH,NORMALIZATION_REGION_OFFSET])
        
    else:
        print("Aborting scheduler preprocessing submission")
        

def print_curated_data_summary(experimental_details,experimentFolderPath):
    #Do not edit
    #curated_data folder should exist
    experimentDate=experimental_details[0]["date"]
    experiment_summary_table = pd.read_csv('%s/curated_data/positions.csv'%(experimentFolderPath))
    return(experiment_summary_table)
    
def create_moma_batch_directory(experimentFolderPath,batch_folder):
    if not os.path.exists("%s/%s"%(experimentFolderPath,batch_folder)):
        os.mkdir("%s/%s"%(experimentFolderPath,batch_folder))
        print("Making batch directory %s."%(batch_folder))
    else:
        print("Batch folder already exist.")

    if not os.path.exists("%s/%s/mm.properties"%(experimentFolderPath,batch_folder)):
        shutil.copy("./mm.properties","%s/%s"%(experimentFolderPath,batch_folder))
        print("Copying standard mm.properties file to %s."%(batch_folder))
    else:
        print("mm.properties already exists, copying cancelled.")
        
def create_yaml_file(batch_file,curation_id,experimentFolderPath,preprocFolder,batch_folder,subset):
    preproc_dir="%s/%s"%(experimentFolderPath,preprocFolder)
    print(preproc_dir)
    if os.path.exists("%s/%s/%s"%(experimentFolderPath,batch_folder,batch_file)):
        #os.remove("%s/%s/%s"%(experimentFolderPath,batch_folder,batch_file))
        print("Batch YAML file already exist. Aborting")
    
    else:
        print("Batch file does not exist. It is created for the following positions, for all detected growth lanes:")

        f=open("%s/%s/%s"%(experimentFolderPath,batch_folder,batch_file),"w")
        f.write("file_format: 0.1.0\n")
        f.write("preprocessing_path: \'%s\'\n"%(preproc_dir))
        f.write("default_moma_arg:\n")
        f.write("  p: \'%s/%s/mm.properties\'\n"%(experimentFolderPath,batch_folder))
        f.write("  analysis: \'%s\'\n"%(curation_id))
        f.write("pos:\n")

        for pos in list(subset["pos"]):
            print("Pos%i"%(pos))
            f.write("  %i:\n    gl:\n"%(pos))
            for root, dirs, files in walklevel("%s/%s"%(preproc_dir,"Pos%i"%(pos)),0):
                gls=sorted([int(e[6+len(str(pos)):]) for e in dirs])
                #print([e for e in dirs])
                for i in gls:
                    f.write("      %i: {}\n"%(i))
        f.close()
        
def write_export_yaml_file(batch_file,experimentFolderPath,preprocFolder,batch_folder,gls_dic):
    preproc_dir="%s/%s"%(experimentFolderPath,preprocFolder)
    if os.path.exists("%s/%s/%s_%s"%(experimentFolderPath,batch_folder,"export_curated",batch_file)):
        #os.remove("%s/%s/%s"%(experimentFolderPath,batch_folder,batch_file))
        print("Batch YAML file already exist. Aborting")
    
    else:
        print("Batch file does not exist. It is created for the following positions, for all curated growth lanes in the export dictionnary:")

        f=open("%s/%s/%s_%s"%(experimentFolderPath,batch_folder,"export_curated",batch_file),"w")
        f.write("file_format: 0.1.0\n")
        f.write("preprocessing_path: \'%s\'\n"%(preproc_dir))
        f.write("default_moma_arg:\n")
        f.write("  p: \'%s/%s/mm.properties\'\n"%(experimentFolderPath,batch_folder))
        f.write("  analysis: \'%s\'\n"%(batch_file[:-5]))
        f.write("pos:\n")

        for pos in list(gls_dic.keys()):
            print("Pos%i"%(pos))
            f.write("  %i:\n    gl:\n"%(pos))
            for i in gls_dic[pos]:
                f.write("      %i: {}\n"%(i))
        f.close()
        
def generate_import_script_input(experimentFolderPath,subset,experiment_summary_table,experimental_details):
    proms=list(subset['promoter'])
    poss=list(subset['pos'])
    unique_proms=list(set(proms))
    strainToPosDic={s:[poss[i] for i in range(0,len(proms)) if proms[i]==s ] for s in unique_proms}

    #print("date,description,f_start,f_end,condition,t_interval,data_path,promoter,vector,cell_detection_offset\n")

    experiment_summary_table_reduced=experiment_summary_table.drop_duplicates(subset=['promoter','vector'])

    for entry in experimental_details:
        #print(entry)
        proms=list(experiment_summary_table_reduced['promoter'])
        vecs=list(experiment_summary_table_reduced['vector'])
        for k in range(0,len(proms)):
            print(entry["date"],entry["description"],entry["f_start"],entry["f_end"],entry["condition"],entry["t_interval"],experimentFolderPath+"/curated_data/%s_%s_curated"%(proms[k],vecs[k]),proms[k],vecs[k],entry["cell_detection_offset"],sep=",")
                   
  
