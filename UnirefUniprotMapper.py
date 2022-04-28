"""
Created on Wed Apr 27 16:27:34 2022

@author: hbckleikamp
"""
#%% parameters

accs=[] #put here your Uniref accessions

input_data_type='N50' #Uniref50, please look at abbreviations from https://www.uniprot.org/help/api_idmapping

columns=[ #these are the columns that you would like to have in your output data, please look at: https://www.uniprot.org/help/uniprotkb_column_names
    "id",
    "entry%20name",
    "protein%20names",
    "organism-id",
    "lineage(SUPERKINGDOM)",
    "lineage(PHYLUM)",
    "lineage(CLASS)",
    "lineage(ORDER)",
    "lineage(FAMILY)",
    "lineage(GENUS)",
    "lineage(SPECIES)",
    "ec",
    "go"
    
    
    ]



#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())

basedir=os.getcwd()


import pandas as pd
import time
import threading
import urllib

#%% functions

#Uniprot mapping multithreaded
def uniprot_mapping_mt(fromtype, totype, identifier,rs):
    """Takes an identifier, and types of identifier 
    (to and from), and calls the UniProt mapping service"""
    base = 'http://www.uniprot.org'
    tool = 'mapping'
    params = {'from':fromtype,
                'to':totype,
                'format':'tab',
                'query':identifier,
    }
    #urllib turns the dictionary params into an encoded url suffix
    data = urllib.parse.urlencode(params)
    #construct the UniProt URL
    url = base+'/'+tool+'?'+data
    url=url.replace("%2B","+")
    #and grab the mapping
    
    while True:
    
        try:
            r = urllib.request.urlopen(url).read()
            rs.append(pd.DataFrame([str(i).split("\\t") for i in str(r).split("\\n")],
                                   columns=["From","To"]))
        
            break
        except:
            "sleeping"
            time.sleep(2)



#uniprot mapping with collumn argument
def uniprot_mapping_mtc(fromtype, totype, identifier,columns,rs):
    """Takes an identifier, and types of identifier 
    (to and from), and calls the UniProt mapping service"""
    base = 'http://www.uniprot.org'
    tool = 'mapping'
    params = {'from':fromtype,
                'to':totype,
                'format':'tab',
                'query':identifier,
                'columns':columns,
    }
    #urllib turns the dictionary params into an encoded url suffix
    data = urllib.parse.urlencode(params)
    #construct the UniProt URL
    url = base+'/'+tool+'?'+data
    url=url.replace("%2B","+")
    #and grab the mapping
    # while True:
    
    #     try:
    r = urllib.request.urlopen(url).read()
    yourlist=str(r).split("\\n")[0].split("b'")[1] #get the batch submission webid
    
    
    #urllib turns the dictionary params into an encoded url suffix
    data = urllib.parse.urlencode(params)
    url = "https://www.uniprot.org/uniprot/?query="+yourlist+"&format=tab&columns="+",".join(columns)
    
    
    r = urllib.request.urlopen(url).read()
    
    
    rs.append(pd.DataFrame([str(i).split("\\t") for i in str(r).split("\\n")],
                            columns=columns))
            
        #     break
        # except:
        #     "sleeping"
        #     time.sleep(2)



# chunker
def chunker(lst,n):
    for i in range(0,len(lst),n):
        yield lst[i:i+n]
        
        
def lca(x,rank_names):
        ix=0
        for r in rank_names:
            if x[r].nunique()!=1:
                break
            else:
                ix+=1
        blca=x.iloc[0,0:ix].tolist()
        return blca+[""]*(len(rank_names)-len(blca)) #pad






#%% Map Uniref to UniprotKB
s=time.time()
rs=list()
threads=[]
counter=0
base_thread=threading.active_count()

batchsize=200
chunks=chunker(accs,batchsize)
for chunk in chunks: #multithread this!
    time.sleep(0.1)
    counter+=1
    print(counter)
    t=threading.Thread(target=uniprot_mapping_mt, args=[input_data_type,
                                                    'ACC',
                                                    "+OR+".join(chunk.tolist()),
                                                    rs])
    t.start()
    threads.append(t)
    
    #unwind in case of thread overload, manage server traffic
    if counter%25==0:
        print("unwinding, query at: "+str(counter/(len(accs)//batchsize))+" elapsed time="+str(time.time()-s))
        for thread in threads:
            thread.join()
        threads=[] #this seems to act different on windows?
         
for thread in threads:
    thread.join()
    

rdf=pd.concat(rs)
rdf=rdf[~rdf.applymap(lambda x: x is None).any(axis=1)]
rdf=rdf[~rdf.applymap(lambda x: x.startswith("b'")).any(axis=1)]


            

#%% Map UniprotKB accessions to lineages and functions

s=time.time()
rs=list()
threads=[]
counter=0

accs=rdf["To"].drop_duplicates()
batchsize=200
chunks=chunker(accs,batchsize)
for chunk in chunks: #multithread this!
    
    time.sleep(0.2)
    counter+=1
    print(counter)
    

  
    t=threading.Thread(target=uniprot_mapping_mtc, args=['ACC',
                                                    'ACC',
                                                    "+OR+".join(chunk.tolist()),
                                                    columns,
                                                    rs])
    t.start()
    threads.append(t)
    
    #unwind in case of thread overload, manage server traffic
    if counter%50==0:
        print("unwinding, query at: "+str(counter/(len(accs)//batchsize))+" elapsed time="+str(time.time()-s))
        for thread in threads:
            thread.join()
        threads=[] #this seems to act different on windows?
      
for thread in threads:
    thread.join()
#%%
tdf=pd.concat(rs)
tdf=tdf[~tdf.applymap(lambda x: x is None).any(axis=1)]
tdf=tdf[~tdf.applymap(lambda x: x.startswith("b'")).any(axis=1)]
mdf=rdf.merge(tdf,how="left",left_on="To",right_on="id") #result is stored in mdf

mdf.to_csv("mapped_result.csv")



