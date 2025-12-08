#!/usr/bin/python3
from collections import Counter,defaultdict
import os,sys,csv
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

def draw_PID_cutoff(PID,name,figurename,colors):
    ns,bins,patches = plt.hist(np.log10(PID['count']),bins=15,edgecolor='black',log=True,color=colors[0])
    for i,n in enumerate(ns[5:],5):
        if i==len(ns)-1:
            thre=10**bins[4]
            break
        if n < ns[i+1]:
            thre = 10**bins[i-1]
            if thre > 100:
                thre = 100
            break
    for i,n in enumerate(ns):
        if 10**bins[i] > thre:
            patches[i].set_facecolor(colors[1])
            
    plt.yscale('log')
    plt.ylabel('number of PID')
    plt.xlabel('log10 PID read depth')
        #plt.ylim([0.5,2*max(np.log10(PID['count']))])
    plt.title(name)
    
    goodPID = PID[PID['count']>thre]['PID'].values
    plt.tight_layout()
    plt.savefig(figurename)
    
    return (thre,goodPID)

def draw_PID_confidence(tdata,name,figurename,colors):
    pidict = defaultdict(dict)
    bcdict = defaultdict(int)
    for i in range(len(tdata)):
        bc = tdata.iloc[i]['barcode']
        pid = tdata.iloc[i]['PID']
        c = tdata.iloc[i]['count']
        pidict[pid][bc] = pidict[pid].get(bc,0)+c
    xlist = []
    freqlist = []

    for pid in pidict:
        count = sum(pidict[pid].values()) #total count
        xlist.append(count)
        if len(pidict[pid] )== 1:
            bc = list(pidict[pid].keys())[0]
            freqlist.append(1)
            bcdict[bc] += 1
            continue
        top=sorted((v,k) for k,v in pidict[pid].items())[-1]
        freqlist.append(float(top[0]/count))
        bc=top[1]
        bcdict[bc]+=1
    
    if len(xlist) < 1000: 
        size = 3; alpha = 1
    else:
        size = 1; alpha = 0.2
    plt.figure(figsize=[5,5])
    plt.scatter(xlist,freqlist,alpha=alpha,s=size,color=colors[0])
    plt.xlim(1,max(xlist)*3)
    plt.ylim([0,1.05])
    plt.xlabel('read depth')
    plt.ylabel('PID correction confidence')
    plt.xscale('log')
    plt.title(name)
    plt.tight_layout()
    plt.savefig(figurename)

    bcdict= {k:v for (v,k) in sorted(((v,k) for k,v in bcdict.items()),reverse=True)}    
    return bcdict

def cluster_bc(index,bcdict,path):
    clusterdict = {}
    bcindex = 0
    tmpindexdict = {}
    
    tempinfile = os.path.join(path,'temp',index)
    tempoutfile = os.path.join(path,'temp','{}_out'.format(index))
    clusterfile = os.path.join(path,'temp','{}_out.clstr'.format(index))
    
    with open(tempinfile,'w') as f:
        for bc in bcdict:
            bcindex += 1
            f.write('>bc'+str(bcindex)+'\n'+bc+'\n')
            tmpindexdict[bcindex] = bc
   
    os.system('cd-hit-est -i {} -o {} -c 0.8'.format(tempinfile,tempoutfile))
    
    with open(clusterfile) as f:
        Cluster = {}
        for line in f:
            if line[0] == '>': continue
            if line[0] == '0':
                key = line.rsplit('>bc')[-1].rsplit('...')[0]
                Cluster[key] = [key]
            else:
                kid = line.rsplit('>bc')[-1].rsplit('...')[0]
                Cluster[key].append(kid)

    for key in Cluster:
        bc = tmpindexdict[int(key)]
        clusterdict[bc] = 0
        for kid in Cluster[key]:
            clusterdict[bc] += bcdict[tmpindexdict[int(kid)]]
    clusterdict = dict(sorted(clusterdict.items(), key=lambda kv: -kv[1]))
    
    return clusterdict

def save_parameter(index,path,folder,depth,goodPID,thre,bcdict,clusterdict):
    info={'depth':int(depth),'copy num':len(goodPID),'PID threshold':float(thre),'good PID':list(goodPID),'bc dict':bcdict,'bc count':len(bcdict),'bc cluster count':len(clusterdict),'cluster dict':clusterdict}
    with open(os.path.join(path,folder,'analysis','info_{}.json'.format(index)),'w') as f:
        json.dump(info,f)
    #save barcode files
    with open(os.path.join(path,folder,'analysis','barcode_{}.csv'.format(index)),'w') as csvfile:
        fieldnames = ['barcode', 'count']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        for barcode, count in bcdict.items():
            writer.writerow({'barcode': barcode, 'count': count})
    #save cluster files:
    with open(os.path.join(path,folder,'analysis','cluster_{}.csv'.format(index)),'w') as csvfile:
        fieldnames = ['barcode', 'count']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        for barcode, count in clusterdict.items():
            writer.writerow({'barcode': barcode, 'count': count})
    return (1)

def main():
    index=sys.argv[2]
    folder='NovaSeq053023'
    reffile=sys.argv[1]

    path='/u/scratch/y/yuanshi'
    subfolders=['analysis','temp','figures']
    for subfolder in subfolders:
        Path(os.path.join(path,folder,subfolder)).mkdir(parents=True, exist_ok=True)
    reftb=pd.read_csv(reffile).set_index('index')
    name=reftb.loc[index]['name']
    #get count file
    tdata=pd.read_csv(os.path.join(path,folder,'output','{}_count.txt'.format(index)),sep='\t',names=['barcode','PID','count'])

    #correct bc
    def correct_bc(bc): #if 3rd position is N, change to C
        pos1=bc[::3]
        pos2=bc[1::3]
        pos3=bc[2::3].replace('N','C')
        for i,j,k in zip(pos1,pos2,pos3):
            yield(i+j+k)

    tdata['barcode_correct'] = tdata['barcode'].apply(lambda x: ''.join(list(correct_bc(x))))
    #filter out 'N'
    mask1 = tdata['barcode_correct'].apply(lambda x: 'N' not in x)
    mask2 = tdata['PID'].apply(lambda x: 'N' not in x)
    tdata = tdata[mask1 & mask2][['barcode_correct','PID','count']]
    tdata.rename(columns={'barcode_correct': 'barcode'}, inplace=True)
    tdata = tdata.groupby(['barcode', 'PID']).sum().reset_index()
    depth = np.sum(tdata['count'])
    PID = tdata.groupby('PID').sum().reset_index()

    tdata.to_csv(os.path.join(path,folder,'analysis','data_{}.csv'.format(index)),index=False)
    PID.to_csv(os.path.join(path,folder,'analysis','UMI_{}.csv'.format(index)),index=False)
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    # PID cutoff
    figurename=os.path.join(os.path.join(path,folder,'figures','PID_cutoff_{}'.format(index)))
    thre,goodPID=draw_PID_cutoff(PID,name,figurename,colors)

    # PID confidence
    tdata = tdata[tdata['PID'].isin(goodPID)]
    tdata.to_csv(os.path.join(path,folder,'analysis','filtered_data_{}.csv'.format(index)),index=False)
    figurename=os.path.join(path,folder,'figures','PID_confidence_{}'.format(index))
    bcdict=draw_PID_confidence(tdata,name,figurename,colors)

    # barcode cluster
    clusterdict=cluster_bc(index,bcdict,os.path.join(path,folder))
    
    # save all parameters
    save_parameter(index,path,folder,depth,goodPID,thre,bcdict,clusterdict)


if __name__ == '__main__':
    main()
Text to speech button
