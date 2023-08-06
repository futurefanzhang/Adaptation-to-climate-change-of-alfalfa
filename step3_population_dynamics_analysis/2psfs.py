import pandas as pd
import os
import sys
import matplotlib.pyplot as plt 


def file_sfs(i,o,gp1,gp2):
    dic_g={'M.truncatula':'d0_','M.falcata.diploid':'d1_','M.caerulea':'d2_'} 
    pf  =pd.read_csv(i,sep='\t')
    pf_c=pf.groupby(gp1)
    # c=pd.DataFrame()
    # print(list(range(0,max(pf[gp2])+1)))
    # input()
    # input()
    c=pd.DataFrame(index=list(range(0,max(pf[gp2]+1))))
    # input()
    for x,y in pf_c:
        pf_2=y.groupby(gp2).count()
        c[dic_g[gp1]+str(x)]=pf_2[gp1]
               
    new_index={}
    for row in c.index.to_list():
        new_index[row]=dic_g[gp2]+str(row)
    c.rename(index=new_index, inplace=True)
    # c.index.names = ['']
    # print(c)
    # c.index.name=None
    c.fillna(value='0',inplace=True)
    # print(c)
    c.to_csv(o+gp1+'_'+gp2+'_cache',sep='\t')
    
    with open(o+gp1+'_'+gp2+'_cache', 'r') as f:
        with open(o+'T1_jointDAFpop'+dic_g[gp2][-2]+'_'+dic_g[gp1][-2]+'.obs', 'w') as f1:
            f1.write('1 observation'+'\n')
        for l in f:
            with open(o+'T1_jointDAFpop'+dic_g[gp2][-2]+'_'+dic_g[gp1][-2]+'.obs', 'a') as f1:
                f1.write(l)
    os.remove(o+gp1+'_'+gp2+'_cache')
        


i=r'/public/agis/zhouyongfeng_group/zhangfan02/vcf_file/count'
o=r'/public/agis/zhouyongfeng_group/zhangfan02/vcf_file/Q'
file_sfs(i,o,'M.truncatula','M.falcata.diploid')
file_sfs(i,o,'M.truncatula','M.caerulea')
file_sfs(i,o,'M.falcata.diploid','M.caerulea')
