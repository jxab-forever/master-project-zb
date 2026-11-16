"""Data-driven Matplotlib figures; no hard-coded solution values."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter

TECHS = ['Coal','GasCT','GasCC','Oil','Biomass','Nuclear','Hydro','Wind','SolarPV','SolarTH','Geo']
COLORS = ['#354052','#ef9f43','#cc6138','#8a6855','#b69870','#a28cc6','#4b88b5','#479e85','#ddc259','#71bfc0','#bd79a2']


def style():
    plt.rcParams.update({'font.family':['Microsoft YaHei','DejaVu Sans'],
        'font.size':10,'axes.titlesize':12,'axes.titleweight':'bold',
        'axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':'#b8c1cc',
        'axes.labelcolor':'#354052','text.color':'#253247','xtick.color':'#59677b',
        'ytick.color':'#59677b','grid.alpha':0.18,'grid.color':'#8795a7',
        'figure.facecolor':'#fafbfe','axes.facecolor':'white',
        'svg.fonttype':'none','axes.unicode_minus':False})


def save(fig, path):
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path.with_suffix('.png'),dpi=180,bbox_inches='tight',facecolor=fig.get_facecolor())
    fig.savefig(path.with_suffix('.svg'),bbox_inches='tight',facecolor=fig.get_facecolor())
    plt.close(fig)


def plot_prototype(folder):
    style()
    folder=Path(folder)
    cap=pd.read_csv(folder/'capacity_by_tech.csv')
    gen=pd.read_csv(folder/'generation_by_tech.csv')
    co2=pd.read_csv(folder/'co2_emissions.csv')
    dem=pd.read_csv(folder/'demand.csv')
    summary=json.loads((folder/'summary.json').read_text(encoding='utf8'))
    fig,axes=plt.subplots(2,2,figsize=(13,8.8),layout='constrained')
    fig.suptitle(f"Prototype 电力系统 · Python / HiGHS\n最优目标值 {summary['objective']:,.6f} 十亿美元  |  31 年 · 11 类技术",fontsize=16,fontweight='bold')
    ax=axes[0,0]
    ax.stackplot(cap.year,cap[TECHS].to_numpy().T,colors=COLORS,alpha=.92)
    ax.set(title='有效运行容量',xlabel='模型年份',ylabel='GW',xlim=(0,30))
    ax.grid(axis='y')
    ax=axes[0,1]
    ax.stackplot(gen.year,gen[TECHS].to_numpy().T,colors=COLORS,alpha=.92)
    demand_column=next(c for c in dem if c!='year')
    ax.plot(dem.year[1:],dem[demand_column][1:],color='#1b2940',ls='--',lw=1.7,label='Demand')
    ax.set(title='发电量与需求',xlabel='模型年份',ylabel='TWh / 年',xlim=(0,30))
    ax.grid(axis='y')
    ax=axes[1,0]
    ax.plot(co2.year,co2.co2_Mt,color='#b85551',lw=2.3)
    ax.fill_between(co2.year,co2.co2_Mt,color='#b85551',alpha=.12)
    ax.set(title=f"年度 CO₂ 排放 · 累计 {summary['total_co2']:,.1f} Mt",xlabel='模型年份',ylabel='Mt CO₂ / 年',xlim=(0,30))
    ax.grid(axis='y')
    ax=axes[1,1]
    positions=np.arange(len(TECHS))
    ax.bar(positions-.19,cap[TECHS].iloc[0],width=.38,color='#9cabbf',label='Year 0')
    ax.bar(positions+.19,cap[TECHS].iloc[-1],width=.38,color=COLORS,label='Year 30')
    ax.set(title='起点与终点容量',ylabel='GW',xticks=positions,xticklabels=TECHS)
    ax.tick_params(axis='x',labelrotation=50)
    ax.legend(frameon=False)
    ax.grid(axis='y')
    handles=[plt.Line2D([0],[0],color=c,lw=6,label=t) for t,c in zip(TECHS,COLORS)]
    axes[0,0].legend(handles=handles,ncol=4,fontsize=8,loc='upper left',frameon=False)
    save(fig,folder/'figures/prototype_dashboard')
    fig,ax=plt.subplots(figsize=(11,5.5),layout='constrained')
    ax.stackplot(gen.year,gen[TECHS].to_numpy().T,colors=COLORS,labels=TECHS)
    ax.plot(dem.year[1:],dem[demand_column][1:],color='#253247',ls='--',lw=2,label='Demand')
    ax.set(title='Prototype · 发电技术组合与需求',xlabel='模型年份',ylabel='TWh / 年',xlim=(0,30))
    ax.legend(ncol=4,loc='upper left',fontsize=9,frameon=False)
    ax.grid(axis='y')
    save(fig,folder/'figures/generation_mix')


def plot_softwarex(folder):
    style()
    folder=Path(folder)
    d=pd.read_csv(folder/'annual_results.csv')
    choices=pd.read_csv(folder/'technology_choices.csv')
    location=pd.read_csv(folder/'location_results.csv')
    summary=json.loads((folder/'summary.json').read_text(encoding='utf8'))
    years=d.year.to_numpy()
    gap=summary['mip_gap']*100
    status='已证明最优' if summary['status']=='Optimal' and gap<1e-6 else f"可行解 · 最优性差距 {gap:.2f}%"
    fig,axes=plt.subplots(2,2,figsize=(13,8.8),layout='constrained')
    fig.suptitle(f"SoftwareX 工业转型 · Python / HiGHS\n目标值 {summary['objective_MUSD']:,.2f} 百万美元  |  {status}",fontsize=16,fontweight='bold')
    ax=axes[0,0]
    block_years=np.r_[years,years[-1]+d.years.iloc[-1]]
    ax.stackplot(block_years,np.r_[d.existing_tonnes_per_year,d.existing_tonnes_per_year.iloc[-1]],
                 np.r_[d.new_tonnes_per_year,d.new_tonnes_per_year.iloc[-1]],
                 labels=['现有设施（含改造）','新建设施'],colors=['#4c87af','#60aa96'],step='post')
    ax.step(block_years,np.r_[d.demand_tonnes_per_year,d.demand_tonnes_per_year.iloc[-1]],
            where='post',color='#253247',ls='--',label='需求')
    ax.set(title='产出与需求',xlabel='年份（每个子期持续五年）',ylabel='吨 / 年',xlim=(years[0],block_years[-1]))
    ax.grid(axis='y');ax.legend(frameon=False,fontsize=9)
    ax=axes[0,1]
    cumulative=np.cumsum(d.co2_tonnes_per_year*d.years)
    ax.plot(years+5,cumulative,'o-',color='#b85551',lw=2)
    ax.axhline(summary['cumulative_co2_budget_tonnes'],color='#63758a',ls='--',label='全期碳预算')
    ax.set(title='累计 CO₂（直接排放与投入品上游排放）',xlabel='五年子期结束年',ylabel='吨 CO₂')
    ax.grid(axis='y');ax.legend(frameon=False,fontsize=9)
    ax=axes[1,0]
    ax.plot(years,d.electricity_MMBtu_per_year,'o-',color='#4c87af',label='实际用电')
    ax.plot(years,d.electricity_budget_MMBtu_per_year,'--',color='#aab5c4',label='供应上限')
    ax.set(title='电力消耗与供应约束',xlabel='五年子期起始年',ylabel='MMBtu / 年')
    ax.grid(axis='y');ax.legend(frameon=False,fontsize=9)
    ax=axes[1,1]
    captured=d.o_ep1gce_tonnes_per_year+d.o_ep1gcs_tonnes_per_year+d.n_ep1gce_tonnes_per_year+d.n_ep1gcs_tonnes_per_year
    ax.plot(years,captured,'o-',color='#60aa96',label='实际捕集')
    ax.axhline(d.capture_budget_tonnes_per_year.iloc[0],ls='--',color='#aab5c4',label='捕集与储存供应上限')
    ax.set(title='CO₂ 捕集及供应约束',xlabel='五年子期起始年',ylabel='吨 CO₂ / 年')
    ax.grid(axis='y');ax.legend(frameon=False,fontsize=9)
    save(fig,folder/'figures/softwarex_dashboard')
    label_translations={'No retrofit':'原工艺','5% energy efficiency':'能效改造',
        'Coal->NG':'燃料切换','Full electrification':'电气化','90% carbon capture':'碳捕集',
        'Alt. Pathway':'替代路径','Incumbent':'原工艺','New plant BAT':'先进工艺'}
    fig,axes=plt.subplots(2,1,figsize=(12,7.4),layout='constrained')
    for ax,form,title in zip(axes,['y_r','y_n'],['现有设施的改造技术选择','新建设施的技术选择']):
        z=choices[choices.form==form].copy()
        labels=dict(z[['technology','label']].drop_duplicates().itertuples(index=False,name=None))
        if form=='y_r':
            z=z.merge(location[['year','location','existing_online']],on=['year','location'])
            z.loc[z.existing_online<=.5,'activity']=0
        selected=z[z.activity>0.5]
        table=selected.pivot(index='location',columns='year',values='technology').reindex(index=range(1,6),columns=years).fillna(0)
        if form=='y_n': table=table.replace(1,0)
        from matplotlib.colors import ListedColormap,BoundaryNorm
        cmap=ListedColormap(['#edf0f5','#cad0da','#75a7b9','#d9b56c','#74b3a1','#a28cc6','#cb8c8c'])
        ax.imshow(table.to_numpy(),cmap=cmap,norm=BoundaryNorm(np.arange(-.5,7),7),aspect='auto')
        for (i,j),v in np.ndenumerate(table.to_numpy()):
            text='停用' if v==0 else label_translations.get(labels[int(v)],labels[int(v)])
            ax.text(j,i,text,ha='center',va='center',fontsize=9)
        ax.set(title=title,xticks=range(len(years)),xticklabels=years,
               yticks=range(5),yticklabels=[f'设施 {i}' for i in range(1,6)],xlabel='五年子期起始年')
    save(fig,folder/'figures/technology_pathways')
    return summary


def plot_all(root):
    root=Path(root)
    if (root/'prototype/summary.json').exists(): plot_prototype(root/'prototype')
    if (root/'softwarex/summary.json').exists(): plot_softwarex(root/'softwarex')
