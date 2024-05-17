import matplotlib.patches as mpatches
from collections import Counter
import logomaker
import ast
import numpy as np
from scipy import stats
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

def fig_1():
    fig, axs = plt.subplots(2,2, sharey='row')
    
    X, Y = dfg['Observed mutation frequency'],dfg['DeepSHM predicted mutation frequency']
    slope, intercept, r, p, se = stats.linregress(X,Y)

    axs[0,0].plot(np.linspace(0,.6,5),np.linspace(0,.6,5),c='k',linewidth=1)
    axs[0,0].plot(X,[i*slope+intercept for i in X],c='r',linewidth=1)
    axs[0,0].set_title('A'+r'$\bf{G}$'+'CT central G, r = %.2f, $\it{P}$<10$^{%s}$' % (r,str(p).split('e')[1]))
    sns.scatterplot(ax=axs[0,0],data=dfg,x='Observed mutation frequency',y='DeepSHM predicted mutation frequency')

    X,Y =  dfc['Observed mutation frequency'],dfc['DeepSHM predicted mutation frequency']
    slope, intercept, r, p, se = stats.linregress(X,Y)

    axs[0,1].plot(np.linspace(0,.6,5),np.linspace(0,.6,5),c='k',linewidth=1)
    axs[0,1].plot(X,[i*slope+intercept for i in X],c='r',linewidth=1)
    axs[0,1].set_title('AG'+r'$\bf{C}$'+'T central C, r = %.2f, $\it{P}$<10$^{%s}$' % (r,str(p).split('e')[1]))
    sns.scatterplot(ax=axs[0,1],data=dfc,x='Observed mutation frequency',y='DeepSHM predicted mutation frequency')
    
    sns.violinplot(ax=axs[1,0],data=df,x='IMGT Subregion',y='Observed mutation frequency',order=['FW1','CDR1','FW2','CDR2','FW3'],color='lightblue',linewidth=1,s=5)
    axs[1,0].set_title('Observed mutation frequency')
    axs[1,0].set_xlabel('')

    sns.violinplot(ax=axs[1,1],data=df,x='IMGT Subregion',y='DeepSHM predicted mutation frequency',order=['FW1','CDR1','FW2','CDR2','FW3'],color='lightblue',linewidth=1,s=5)
    axs[1,1].set_title('DeepSHM predicted mutation frequency')
    axs[1,1].set_xlabel('')

    plt.tight_layout()
    plt.savefig('figures/Figure_1.png',dpi=300)

    plt.show()

def fig_2():
    fig, axs = plt.subplots(15,2,gridspec_kw={'height_ratios':[3,1,.25]*5},figsize=(7,6))
    subregs = ['FW1', 'CDR1', 'FW2', 'CDR2', 'FW3']

    for i, subreg in enumerate(subregs):
        seqs = dfg[dfg['IMGT Subregion']==subreg]['sequence'].tolist()
        seqs = logomaker.alignment_to_matrix(seqs)

        sns.boxplot(data=np.array([ast.literal_eval(i) for i in dfg[dfg['IMGT Subregion']==subreg]['explanation'].tolist()]),ax=axs[i*3,0],linewidth=1,color='lightblue',flierprops={"marker": "."},fliersize=2)
        logomaker.Logo(seqs,ax=axs[i*3+1,0])
        axs[i*3,0].set_ylabel(subreg)
        axs[i*3+1,0].set_yticklabels([])

        seqs = dfc[dfc['IMGT Subregion']==subreg]['sequence'].tolist()
        seqs = logomaker.alignment_to_matrix(seqs)
        sns.boxplot(data=np.array([ast.literal_eval(i) for i in dfc[dfc['IMGT Subregion']==subreg]['explanation'].tolist()]),ax=axs[i*3,1],linewidth=1,color='lightblue',flierprops={"marker": "."},fliersize=2)
        logomaker.Logo(seqs,ax=axs[i*3+1,1])

        axs[i*3,0].sharey(axs[i*3,1])
        axs[i*3,1].tick_params(labelleft=False)
        axs[i*3+1,0].set_yticklabels([])
        axs[i*3+1,1].set_yticklabels([])
        axs[i*3+2,0].remove()
        axs[i*3+2,1].remove()

    axs[0,0].set_title('A'+r'$\bf{G}$'+'CT central G mutated')
    axs[0,1].set_title('AG'+r'$\bf{C}$'+'T central C mutated')
    axs[13,0].set_xticklabels([])
    axs[13,1].set_xticklabels([])

    fig.supylabel('Integrated gradient score')
    plt.subplots_adjust(wspace=0.05, hspace=0)
    plt.savefig('figures/Figure_2.png',dpi=300)
    plt.show()

def fig_3():
    df = pd.concat([dfg,dfc])
    sns.scatterplot(data=df,x="5' Integrated gradients score",y="3' Integrated gradients score",style='Motif',style_order=['CAGCTG','non-CAGCTG'],color='black')
    plt.axvline(0,color='b',linestyle='dashed',lw=1)
    plt.axhline(0,color='b',linestyle='dashed',lw=1)
    plt.title('Integrated gradients scores for nucleotides flanking AGCT across IGHVs')
    plt.tight_layout()
    plt.savefig('figures/Figure_3.png',dpi=300)
    plt.show()

def fig_4():
    df = pd.concat([dfg,dfc])
    colors = {"FW1":"tab:blue",
              "CDR1":"tab:orange", 
              "FW2":"tab:green",
              "CDR2":"tab:brown",
              "FW3":"tab:red"}
    
    flank_mfs = {}
    for i in set(df['Flanking Bases'].tolist()):
        flank_mfs[i] = df[df['Flanking Bases']==i]['Observed mutation frequency'].mean()
    flank_mfs = dict(sorted(flank_mfs.items(), key=lambda x: x[1]))
    g = sns.swarmplot(data=df,x='Flanking Bases',y='Observed mutation frequency',hue='IMGT Subregion',order=flank_mfs.keys(),palette=colors,hue_order=['FW1','CDR1','FW2','CDR2','FW3'])
    g.set_xticklabels(["5'-"+i[0]+", 3'-"+i[1] for i in flank_mfs.keys()],rotation=45)
    g.set_xlabel("5' and 3' nucleotides flanking AGCT motif")
    g.set_title('A'+r'$\bf{GC}$'+'T mutation frequency by flanking nucleotide identity')
    plt.tight_layout()
    plt.savefig('figures/Figure_4.png',dpi=300)
    plt.show()

def fig_5(supp=False):
    fig, axs = plt.subplots(1,2,sharey=True)
    dfc = pd.read_csv('data/moods_c.csv')
    dfg = pd.read_csv('data/moods_g.csv')

    if supp==True:
        dfc = dfc[dfc['matrix']=='TFAP4.pfm']
        dfg = dfg[dfg['matrix']=='TFAP4.pfm']
    else:
        dfc = dfc[dfc['matrix']=='TCF3.pfm']
        dfg = dfg[dfg['matrix']=='TCF3.pfm']

    dfc['Motif']=['CAGCTG' if i[4:10]=='CAGCTG' else 'non-CAGCTG' for i in dfc['seq'].tolist()]
    dfg['Motif']=['CAGCTG' if i[5:11]=='CAGCTG' else 'non-CAGCTG' for i in dfg['seq'].tolist()]

    slope, intercept, r, p, se = stats.linregress(dfg['Observed mutation frequency'],dfg['MOODS score'])
    sns.scatterplot(ax=axs[0],data=dfg,x='Observed mutation frequency',y='MOODS score',hue='Motif')
    mf = dfg['Observed mutation frequency'].tolist()
    axs[0].plot(mf,[slope*i+intercept for i in mf],c='r',lw=1)

    if supp==True:
        axs[0].set_title('A'+r'$\bf{G}$'+'CT, r = %.3f, p = %.3f' % (r,p))
    else:
        axs[0].set_title('A'+r'$\bf{G}$'+'CT, r = %.2f, $\it{P}$<10$^{%s}$' % (r,str(p).split('e')[1]))
    
    slope, intercept, r, p, se = stats.linregress(dfc['Observed mutation frequency'],dfc['MOODS score'])
    sns.scatterplot(ax=axs[1],data=dfc,x='Observed mutation frequency',y='MOODS score',hue='Motif')
    mf = dfc['Observed mutation frequency'].tolist()
    axs[1].plot(mf,[slope*i+intercept for i in mf],c='r',lw=1)

    if supp==True:
        axs[1].set_title('AG'+r'$\bf{C}$'+'T, r = %.3f, p = %.3f' % (r,p))
    else:
        axs[1].set_title('AG'+r'$\bf{C}$'+'T, r = %.2f, $\it{P}$<10$^{%s}$' % (r,str(p).split('e')[1]))

    axs[0].get_legend().remove()
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.05 , hspace=None)

    if supp==True:
        plt.savefig('figures/Figure_S1.png',dpi=300)    
    else:
        plt.savefig('figures/Figure_5.png',dpi=300)
    
    plt.show()

def fig_6():
    fig, axs = plt.subplots(1,2,sharey=True)

    df1 = pd.read_csv('data/cagctg_subreg.csv',index_col='IGHV')
    df2 = pd.read_csv('data/canntg_subreg.csv',index_col='IGHV')

    cmap = sns.color_palette('muted', 4)
    sns.heatmap(df1,cmap=sns.color_palette('muted', int(df1.values.max())+1),ax=axs[0],cbar=False,yticklabels=False)
    sns.heatmap(df2,cmap=sns.color_palette('muted', int(df2.values.max())+1),ax=axs[1],cbar=False,yticklabels=False)
    axs[0].axhline(39,lw=1,color='black',linestyle='dashed')
    axs[0].axhline(72,lw=1,color='black',linestyle='dashed')
    axs[0].axhline(155,lw=1,color='black',linestyle='dashed')
    axs[0].axhline(210,lw=1,color='black',linestyle='dashed')
    axs[0].axhline(215,lw=1,color='black',linestyle='dashed')
    axs[0].axhline(218,lw=1,color='black',linestyle='dashed')
    axs[1].axhline(39,lw=1,color='black',linestyle='dashed')
    axs[1].axhline(72,lw=1,color='black',linestyle='dashed')
    axs[1].axhline(155,lw=1,color='black',linestyle='dashed')
    axs[1].axhline(210,lw=1,color='black',linestyle='dashed')
    axs[1].axhline(215,lw=1,color='black',linestyle='dashed')
    axs[1].axhline(218,lw=1,color='black',linestyle='dashed')
    axs[0].set_title('CAGCTG')
    axs[1].set_title('CANNTG')
    axs[1].set_ylabel('')

    patches = [ mpatches.Patch(color=cmap[i],label="{l}".format(l=i) , edgecolor='b' ) for i in range(4)]
    legend=plt.legend(handles=patches, bbox_to_anchor=(1.05, 1),loc=2, borderaxespad=0.5, frameon=True,title='Number of E-box Motifs')

    plt.tight_layout()
    plt.savefig('figures/Figure_6.png',dpi=300)
    plt.show()

def fig_7():
    fig, axs = plt.subplots(1,2,sharey=True)

    dfr = pd.read_csv('data/ramos_e2achip_rpkm_500.csv')
    dfro = dfr[dfr['region']=='other']
    dfr = dfr[dfr['region']!='other']
    dfg = pd.read_csv('data/gm12878_e2achip_rpkm_500.csv')
    dfgo = dfg[dfg['region']=='other']
    dfg = dfg[dfg['region']!='other']

    sns.scatterplot(ax=axs[0],data=dfro,y='case',x='control',alpha=0.8,color='darkgray')
    sns.scatterplot(ax=axs[0],data=dfr,y='case',x='control',hue='region',hue_order = ['IGHV4-34','Eμ','TRBV20-1'])

    axs[0].plot([i for i in range(1000)],[i for i in range(1000)],color='black',lw=1)
    axs[0].set_ylim(0.01,1000)
    axs[0].set_xlim(0.01,1000)
    axs[0].set_title('Ramos E2A ChIP-seq')
    axs[0].set_ylabel('log(RPKM) E2A ChIP-seq (500bp bins)')
    axs[0].set_xlabel('log(RPKM) Control ChIP-seq (500bp bins)')
    axs[0].set_yscale('log')
    axs[0].set_xscale('log')
    axs[0].set_aspect('equal')

    sns.scatterplot(ax=axs[1],data=dfgo,y='case',x='control',alpha=0.8,color='darkgray')
    sns.scatterplot(ax=axs[1],data=dfg,y='case',x='control',hue='region',hue_order = ['IGHV3-21','Eμ','TRBV20-1'])

    axs[1].plot([i for i in range(1000)],[i for i in range(1000)],color='black',lw=1)
    axs[1].set_ylim(0.01,1000)
    axs[1].set_xlim(0.01,1000)
    axs[1].set_title('GM12878 E2A ChIP-seq')
    axs[1].set_ylabel('log(RPKM) E2A ChIP-seq (500bp bins)')
    axs[1].set_xlabel('log(RPKM) Control ChIP-seq (500bp bins)')
    axs[1].set_yscale('log')
    axs[1].set_xscale('log')
    axs[1].set_aspect('equal')

    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.1 , hspace=None)
    plt.savefig('figures/Figure_7.png',dpi=300)
    plt.show()

def preprocess_data(df):
    dfg = df[df['sequence id']%2!=0]    
    dfc = df[df['sequence id']%2==0] 

    dfg['Motif']=['CAGCTG' if i[5:11]=='CAGCTG' else 'non-CAGCTG' for i in dfg['sequence'].tolist()]
    dfg["5' Integrated gradients score"] = [ast.literal_eval(i)[5] for i in dfg['explanation']]
    dfg["3' Integrated gradients score"] = [ast.literal_eval(i)[10] for i in dfg['explanation']]
    dfg["AGCT IG Score"] = [ast.literal_eval(i)[6:10] for i in dfg['explanation']]
    dfg['Flanking Bases'] = [(i[5],i[10]) for i in dfg['sequence']] 

    dfc['Motif']=['CAGCTG' if i[4:10]=='CAGCTG' else 'non-CAGCTG' for i in dfc['sequence'].tolist()]
    dfc["5' Integrated gradients score"] = [ast.literal_eval(i)[4] for i in dfc['explanation']]
    dfc["3' Integrated gradients score"] = [ast.literal_eval(i)[9] for i in dfc['explanation']]
    dfc["AGCT IG Score"] = [ast.literal_eval(i)[5:9] for i in dfc['explanation']]
    dfc['Flanking Bases'] = [(i[4],i[9]) for i in dfc['sequence']] 

    return dfg, dfc

if __name__ == '__main__':
    sns.set(style='darkgrid',font='Arial',font_scale=0.8)
    df = pd.read_csv('data/AGCT_pos_truepred_igexps.csv')
    dfg, dfc = preprocess_data(df)

    fig_1()
    fig_2()
    fig_3()
    fig_4()
    fig_5()
    fig_5(supp=True)
    fig_6()
    fig_7()
