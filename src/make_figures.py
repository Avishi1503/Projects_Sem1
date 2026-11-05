"""Create figures and statistical summaries from completed pipeline outputs."""
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import stats
from common import ROOT,json_write

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'white','savefig.facecolor':'white'})
TEAL='#007f86'; ORANGE='#df7634'; INK='#183247'


def rescale_figure():
    path=ROOT/'results/rescale/query.rescaled.hist.tsv'
    frame=pd.read_csv(path,sep='\t')
    fig,ax=plt.subplots(figsize=(10,5.4),layout='constrained')
    q=frame.query_count/frame.query_count.sum()
    ax.plot(frame.length,q,label='Original query',color='#9b83b7',lw=1.5)
    ax.plot(frame.length,frame.reference_probability,label='Reference',color=TEAL,lw=2)
    ax.plot(frame.length,frame.sampled_probability,label='Sampled output',color=ORANGE,lw=1.5,alpha=.85)
    strict=ROOT/'results/rescale/query.strict.hist.tsv'
    if strict.exists():
        second=pd.read_csv(strict,sep='\t')
        ax.plot(second.length,second.sampled_probability,label='Strict shared-support output',color=INK,lw=1,alpha=.75,ls='--')
    ax.set(xlim=(0,700),xlabel='Fragment length (bp)',ylabel='Fraction of records',title='Rescaling the query toward the reference')
    ax.legend(frameon=False); ax.grid(axis='y',alpha=.15)
    fig.savefig(ROOT/'results/rescale/distribution_comparison.png',dpi=180); plt.close(fig)


def vplot_figure():
    base=ROOT/'results/vplot/vplot.matrix.tsv'
    matrix=np.load(base.with_suffix('.npy'))
    info=json.loads(base.with_suffix('.summary.json').read_text())
    fig,ax=plt.subplots(figsize=(10,6),layout='constrained')
    shown=matrix[:min(701,len(matrix))]
    image=ax.imshow(np.ma.masked_equal(shown,0),origin='lower',aspect='auto',extent=[info['offset_low']-.5,info['offset_high']+.5,-.5,len(shown)-.5],cmap='magma',norm=LogNorm(vmin=1,vmax=max(2,int(shown.max()))))
    ax.set(xlabel='Fragment midpoint offset from region midpoint (bp)',ylabel='Fragment length (bp)',title='V-plot | region strand orientation')
    ax.axvline(0,color='white',lw=.6,alpha=.5)
    fig.colorbar(image,ax=ax,label='Region / fragment pairs (log colour scale)')
    fig.savefig(ROOT/'results/vplot/vplot.png',dpi=180); plt.close(fig)


def tcga_figures():
    out=ROOT/'results/tcga'; frame=pd.read_csv(out/'NKX2-1.tsv',sep='\t')
    groups=['Solid Tissue Normal','Primary Tumor']
    # The assignment asks for the expression records; preserve and label this
    # descriptive view. Inference below uses patient-level matched observations.
    fig,ax=plt.subplots(figsize=(8,6),layout='constrained')
    values=[frame.loc[frame.sample_type==g,'log2_tpm_plus_1'].to_numpy() for g in groups]
    bp=ax.boxplot(values,labels=[f'Solid Tissue Normal\nn = {len(values[0])} files',f'Primary Tumor\nn = {len(values[1])} files'],patch_artist=True,showfliers=False,medianprops={'color':INK,'linewidth':2},widths=.5)
    rng=np.random.default_rng(2026)
    for i,(v,c,patch) in enumerate(zip(values,[TEAL,ORANGE],bp['boxes']),1):
        patch.set_facecolor(c); patch.set_alpha(.25)
        ax.scatter(rng.normal(i,.055,len(v)),v,s=9,alpha=.35,color=c,edgecolors='none')
    ax.set(ylabel='log2(TPM + 1)',title='NKX2-1 expression in TCGA-LUAD')
    ax.grid(axis='y',alpha=.15)
    fig.savefig(out/'NKX2-1_boxplot.png',dpi=180); plt.close(fig)
    samples=frame.groupby(['case_id','sample_id','sample_type'],as_index=False).agg(tpm=('tpm','median'),expression_files=('file_id','count'))
    samples['log2_tpm_plus_1']=np.log2(samples.tpm+1)
    samples.to_csv(out/'NKX2-1_biological_samples.tsv',sep='\t',index=False)
    patient=samples.groupby(['case_id','sample_type'],as_index=False).agg(tpm=('tpm','median'))
    patient['log2_tpm_plus_1']=np.log2(patient.tpm+1)
    wide=patient.pivot(index='case_id',columns='sample_type',values='log2_tpm_plus_1')
    paired=wide.dropna(subset=groups)
    paired.to_csv(out/'NKX2-1_matched_patients.tsv',sep='\t')
    delta=paired['Primary Tumor'].to_numpy()-paired['Solid Tissue Normal'].to_numpy()
    summary=dict(gene='NKX2-1',expression_files=len(frame),unique_biological_samples=len(samples),unique_patients=int(patient.case_id.nunique()),groups={g:dict(files=int((frame.sample_type==g).sum()),biological_samples=int((samples.sample_type==g).sum()),patients=int((patient.sample_type==g).sum()),median_tpm=float(frame.loc[frame.sample_type==g,'tpm'].median()),median_log2_tpm_plus_1=float(frame.loc[frame.sample_type==g,'log2_tpm_plus_1'].median())) for g in groups},matched_patients=len(paired),aggregation='Median TPM across files within sample; median sample TPM within patient/tissue; log2(TPM+1) after aggregation')
    if len(delta) and np.any(delta!=0):
        test=stats.wilcoxon(delta,alternative='two-sided',zero_method='wilcox',method='approx')
        summary['paired_wilcoxon']=dict(statistic=float(test.statistic),pvalue=float(test.pvalue),median_tumor_minus_normal_log2=float(np.median(delta)),nonzero_pairs=int(np.count_nonzero(delta)),method='two-sided normal approximation, Wilcox zero handling, no continuity correction')
    else:
        summary['paired_wilcoxon']=dict(pvalue=1.0 if len(delta) else None,reason='all differences zero' if len(delta) else 'no matched pairs')
    unpaired=wide.loc[~wide.index.isin(paired.index)]
    a=unpaired['Primary Tumor'].dropna().to_numpy(); b=unpaired['Solid Tissue Normal'].dropna().to_numpy()
    summary['unmatched_patients']=dict(tumor=len(a),normal=len(b))
    if len(a)>=2 and len(b)>=2:
        test=stats.mannwhitneyu(a,b,alternative='two-sided',method='asymptotic')
        summary['unmatched_mannwhitney']=dict(statistic=float(test.statistic),pvalue=float(test.pvalue))
    else: summary['unmatched_mannwhitney']=dict(pvalue=None,reason='This pipeline requires at least two unmatched patients in each group; one remaining normal patient is insufficient for a useful secondary group comparison')
    json_write(out/'statistics.json',summary)
    if len(paired):
        fig,ax=plt.subplots(figsize=(8,6),layout='constrained')
        for _,row in paired.iterrows(): ax.plot([0,1],[row[g] for g in groups],color='#9caeb7',alpha=.5,lw=.8)
        ax.scatter(np.zeros(len(paired)),paired[groups[0]],color=TEAL,s=18,zorder=3)
        ax.scatter(np.ones(len(paired)),paired[groups[1]],color=ORANGE,s=18,zorder=3)
        p=summary['paired_wilcoxon']['pvalue']
        ax.set(xticks=[0,1],xticklabels=['Solid Tissue Normal','Primary Tumor'],ylabel='log2(TPM + 1)',title=f'NKX2-1 | {len(paired)} matched patients\nPaired Wilcoxon p = {p:.3g}')
        ax.grid(axis='y',alpha=.15)
        fig.savefig(out/'NKX2-1_paired_patients.png',dpi=180); plt.close(fig)
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--project',type=int,choices=[2,3,4]); args=p.parse_args()
    for project,fn in [(2,rescale_figure),(3,vplot_figure),(4,tcga_figures)]:
        if args.project in (None,project): fn()
