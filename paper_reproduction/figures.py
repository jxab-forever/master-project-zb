"""Paper-style quantitative panels rebuilt from Python STREAM solution CSVs.

Visual grammar follows Thierry et al., SoftwareX32(2025)102355, Figs3-6.
All marks represent re-solved model quantities; no article pixels are plotted.
"""
from pathlib import Path
import json
import math
import string
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch,Polygon,Rectangle
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter,MaxNLocator

ROOT=Path(__file__).resolve().parent
FIG=ROOT/'figures'
DATA=ROOT/'source_data'
RESULTS=None
MM=1/25.4
TECH=['Incumbent','Efficiency','Coal to NG','Electrification','CCS','Alternative','BAT']
COLORS=dict(zip(TECH,['#f8b749','#e8cf4f','#b0ca40','#dce6af','#45b0cf','#58ba70','#f2dedb']))
FUEL=['Electricity','Coal','Fuel oil','Natural gas']
FUEL_COLORS=dict(zip(FUEL,['#2b80bd','#7b7b7b','#86917e','#db870f']))


def setup():
    mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
        'font.size':7,'axes.labelsize':7,'axes.titlesize':7,'xtick.labelsize':6,'ytick.labelsize':6,
        'legend.fontsize':6.2,'legend.frameon':False,'axes.linewidth':.5,
        'xtick.major.width':.4,'ytick.major.width':.4,'xtick.major.size':2,'ytick.major.size':2,
        'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','savefig.facecolor':'white',
        'figure.facecolor':'white','axes.facecolor':'white','hatch.linewidth':.32,'mathtext.default':'regular'})
    FIG.mkdir(exist_ok=True);DATA.mkdir(exist_ok=True)


def load(name):
    folder=(RESULTS or ROOT)/'runs'/name
    s=json.loads((folder/'summary.json').read_text(encoding='utf8'))
    d=pd.read_csv(folder/'trajectory.csv')
    d['inactive_candidate_capacity']=0.
    inactive=(d.form=='new')&(d.technology==1)
    d.loc[inactive,'inactive_candidate_capacity']=d.loc[inactive,'installed_capacity']
    d.loc[inactive,'installed_capacity']=0.
    d['tech']=[technology(f,k) for f,k in zip(d.form,d.technology)]
    return s,d,pd.read_csv(folder/'energy.csv'),pd.read_csv(folder/'annual.csv'),pd.read_csv(folder/'states.csv')


def technology(form,k):
    return ({1:'Incumbent',2:'Efficiency',3:'Coal to NG',4:'Electrification',5:'CCS',6:'Alternative'} if form=='existing'
            else {1:'None',2:'Incumbent',3:'BAT',4:'Electrification',5:'CCS',6:'Alternative'})[int(k)]


def panel(ax,label):
    ax.text(-.16,1.06,label,transform=ax.transAxes,fontsize=8.5,fontweight='bold',va='bottom')
    ax.set_xlim(.4,8.6);ax.set_xticks(range(1,9));ax.set_ylim(bottom=0)
    fmt=ScalarFormatter(useMathText=True);fmt.set_powerlimits((0,0));ax.yaxis.set_major_formatter(fmt)
    ax.yaxis.get_offset_text().set_size(5.7);ax.yaxis.set_major_locator(MaxNLocator(5))


def export(fig,stem):
    fig.savefig(FIG/f'{stem}.svg')
    fig.savefig(FIG/f'{stem}.pdf')
    fig.savefig(FIG/f'{stem}.png',dpi=600)
    fig.savefig(FIG/f'{stem}.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
    plt.close(fig)


def tech_legend(ax,loc='upper left',size=6.2):
    ax.axis('off')
    labels=['Incumbent','Efficiency (retrofit)','Coal to NG (retrofit)','Full electrification',
            '90% carbon capture\n& storage','Alternative pathway','Best available technology\n(new build)']
    handles=[Patch(facecolor=COLORS[t],edgecolor='black',linewidth=.35,label=l) for t,l in zip(TECH,labels)]
    ax.legend(handles=handles,title='Technologies',loc=loc,fontsize=size,title_fontsize=6.6,
              borderaxespad=0,labelspacing=.4,handlelength=1.2,handleheight=1.0)


def stack_tech(ax,d,field,by_location=False):
    bottom=np.zeros(8)
    for form in ['existing','new']:
        if by_location:
            groups=[(tech,d[(d.form==form)&(d.location==l)&(d.tech==tech)])
                    for l in range(1,6) for tech in TECH]
        else:
            groups=[(tech,d[(d.form==form)&(d.tech==tech)]) for tech in TECH]
        for tech,z in groups:
            values=z.groupby('subperiod')[field].sum().reindex(range(1,9),fill_value=0).to_numpy()
            ax.bar(range(1,9),values,bottom=bottom,width=.80,color=COLORS[tech],
                edgecolor='#252525',linewidth=.25,hatch='////' if form=='new' else None)
            bottom+=values
    ax.set_ylim(0,max(bottom.max()*1.06,1e-3))
    return bottom


def fig3():
    setup();fig=plt.figure(figsize=(7.204724,5.748031))
    gs=fig.add_gridspec(2,5,left=.092,right=.992,bottom=.285,top=.915,wspace=.69,hspace=.42)
    titles=['Aggregate capacity\nby technology','Production by facility\n& technology',
            'Aggregate energy\nuse by fuel','Direct emissions\nby technology','Emissions by facility\n& fate']
    for row,name in enumerate(['bau_nominal','nominal']):
        s,d,e,a,states=load(name)
        d.assign(case=name).to_csv(DATA/f'fig3_{name}_trajectory.csv',index=False)
        e.assign(case=name).to_csv(DATA/f'fig3_{name}_energy.csv',index=False)
        for col in range(5):
            ax=fig.add_subplot(gs[row,col]);panel(ax,string.ascii_lowercase[row*5+col])
            if row==0:ax.set_title(titles[col],pad=15,fontsize=6.5)
            else:ax.set_xlabel('Subperiod',labelpad=2)
            if col in [0,1,3]:
                field={0:'installed_capacity',1:'production',3:'direct_emissions'}[col]
                stack_tech(ax,d,field,by_location=col==1)
                ax.set_ylabel(r't yr$^{-1}$' if col!=3 else r't CO$_2$ yr$^{-1}$',labelpad=1)
                if col==1:ax.plot(a.subperiod,a.demand,ls='--',lw=.65,color='#9b2a31',zorder=6)
            elif col==2:
                bottom=np.zeros(8)
                order=[(form,'Electricity') for form in ['existing','new']]
                order += [(form,fuel) for form in ['existing','new'] for fuel in FUEL[1:]]
                for form,fuel in order:
                    z=e[(e.fuel==fuel)&(e.form==form)]
                    y=z.groupby('subperiod').paper_energy.sum().reindex(range(1,9),fill_value=0).to_numpy()
                    ax.bar(range(1,9),y,bottom=bottom,width=.8,color=FUEL_COLORS[fuel],
                        hatch='////' if form=='new' else None,edgecolor='#252525',linewidth=.25)
                    bottom+=y
                ax.set_ylim(0,bottom.max()*1.06);ax.set_ylabel(r'MMBtu yr$^{-1}$',labelpad=1)
            else:
                for field,color,sign in [('direct_emissions','#fa6858',1),('captured','#d48ab1',-1)]:
                    bottom=np.zeros(8)
                    for form in ['existing','new']:
                        for l in range(1,6):
                            z=d[(d.form==form)&(d.location==l)]
                            y=z.groupby('subperiod')[field].sum().reindex(range(1,9),fill_value=0).to_numpy()*sign
                            ax.bar(range(1,9),y,bottom=bottom,width=.8,color=color,
                                hatch='////' if form=='new' else None,edgecolor='#252525',linewidth=.25)
                            bottom+=y
                positive=d.groupby('subperiod').direct_emissions.sum().max()
                negative=d.groupby('subperiod').captured.sum().max()
                ax.set_ylim(-negative*1.06 if negative>1e-8 else 0,positive*1.06)
                ax.axhline(0,color='black',lw=.35)
                ax.set_ylabel(r't CO$_2$ yr$^{-1}$',labelpad=1)
    fig.text(.016,.765,'Business as usual',rotation=90,ha='center',va='center',fontsize=7,fontweight='bold')
    fig.text(.016,.407,'50% emission target',rotation=90,ha='center',va='center',fontsize=7,fontweight='bold')
    tech_legend(fig.add_axes([.10,.005,.38,.235]))
    ax=fig.add_axes([.48,.07,.17,.155]);ax.axis('off')
    ax.legend(handles=[Patch(facecolor=FUEL_COLORS[f],edgecolor='black',lw=.3,label=f) for f in FUEL],
        title='Fuels',title_fontsize=6.6,loc='upper left',borderaxespad=0,handlelength=1.1)
    ax=fig.add_axes([.66,.07,.17,.155]);ax.axis('off')
    ax.legend(handles=[Patch(facecolor='#fa6858',edgecolor='black',lw=.3,label=r'CO$_2$ emitted'),
         Patch(facecolor='#d48ab1',edgecolor='black',lw=.3,label=r'CO$_2$ captured')],
         title='Emission fate',title_fontsize=6.6,loc='upper left',borderaxespad=0,handlelength=1.1)
    ax=fig.add_axes([.83,.07,.16,.155]);ax.axis('off')
    ax.legend(handles=[Patch(facecolor='white',edgecolor='black',lw=.3,label='Existing facility'),
         Patch(facecolor='white',edgecolor='black',lw=.3,hatch='////',label='New-build facility'),
         Line2D([0],[0],color='#9b2a31',lw=.65,ls='--',label='Demand')],
         title='Facility type',title_fontsize=6.6,loc='upper left',borderaxespad=0,handlelength=1.1)
    export(fig,'Figure_3_nominal_trajectories')


def albers(lon,lat):
    lon=np.asarray(lon)*np.pi/180;lat=np.asarray(lat)*np.pi/180
    p1,p2,p0=np.deg2rad([29.5,45.5,37.5]);l0=np.deg2rad(-96)
    n=.5*(np.sin(p1)+np.sin(p2));c=np.cos(p1)**2+2*n*np.sin(p1)
    rho=np.sqrt(c-2*n*np.sin(lat))/n;rho0=np.sqrt(c-2*n*np.sin(p0))/n
    theta=n*(lon-l0)
    return rho*np.sin(theta),rho0-rho*np.cos(theta)


def state_polygons():
    geo=json.loads((ROOT/'reference/us_states.geojson').read_text())
    out=[]
    for feature in geo['features']:
        if feature['properties']['name'] in ['Alaska','Hawaii','Puerto Rico']:continue
        geometry=feature['geometry'];polys=geometry['coordinates']
        if geometry['type']=='Polygon':polys=[polys]
        for poly in polys:
            for ring in poly:
                xy=np.asarray(ring);x,y=albers(xy[:,0],xy[:,1]);out.append(np.column_stack([x,y]))
    return out


def fig4():
    setup();s,d,energy,a,states=load('nominal')
    facilities=pd.read_csv(ROOT/'reference/facilities.csv').iloc[:5]
    facilities['location']=np.arange(1,6)
    facilities.to_csv(DATA/'fig4_facilities.csv',index=False)
    states.to_csv(DATA/'fig4_states.csv',index=False)
    matrices=[np.zeros((5,8),int) for _ in range(4)]
    event_rows=[]
    for l in range(1,6):
        z=states[states.location==l].sort_values('subperiod')
        prior_online=1;prior_tech=1;prior_new=1;prior_exp=0.
        for row in z.itertuples():
            t=int(row.subperiod)-1
            events=[row.expansion>prior_exp+1e-5,prior_online>.5 and row.online<.5,
                    row.online>.5 and row.existing_technology!=prior_tech,
                    row.new_technology!=prior_new and row.new_technology!=1]
            for event,flag in enumerate(events):
                if flag:
                    kind=int(row.new_technology if event==3 else row.existing_technology)
                    # Technology of a retirement is the last operating technology.
                    if event==1:kind=prior_tech
                    matrices[event][l-1,t]=TECH.index(technology('new' if event==3 else 'existing',kind))+1
                    event_rows.append(dict(event=event,location=l,subperiod=t+1,technology=matrices[event][l-1,t]))
            prior_online=row.online;prior_tech=row.existing_technology;prior_new=row.new_technology;prior_exp=row.expansion
    pd.DataFrame(event_rows).to_csv(DATA/'fig4_events.csv',index=False)
    fig=plt.figure(figsize=(7.204724,7.047244))
    gs=fig.add_gridspec(4,3,width_ratios=[.78,1.18,1.18],left=.075,right=.995,bottom=.265,top=.965,wspace=.12,hspace=.32)
    from matplotlib.colors import ListedColormap,BoundaryNorm
    cmap=ListedColormap(['white']+[COLORS[t] for t in TECH]);norm=BoundaryNorm(np.arange(-.5,8),8)
    for i,(title,matrix) in enumerate(zip(['Capacity expansion switch','Retirement switch','Retrofit switch','New-build switch'],matrices)):
        ax=fig.add_subplot(gs[i,0]);ax.imshow(matrix,cmap=cmap,norm=norm,aspect='auto',interpolation='none',extent=[.5,8.5,5.5,.5])
        ax.set_xticks(range(1,9));ax.set_yticks(range(1,6));ax.set_ylabel('Plant',labelpad=1)
        ax.set_title(title,pad=4);ax.text(-.17,1.08,string.ascii_lowercase[i],transform=ax.transAxes,fontweight='bold',fontsize=8.5)
        if i==3:ax.set_xlabel('Subperiod')
    polygons=state_polygons();allpoints=np.vstack(polygons)
    xlim=(allpoints[:,0].min()-.02,allpoints[:,0].max()+.02);ylim=(allpoints[:,1].min()-.02,allpoints[:,1].max()+.02)
    for t in range(1,9):
        row,col=(t-1)//2,(t-1)%2+1;ax=fig.add_subplot(gs[row,col])
        for xy in polygons:ax.add_patch(Polygon(xy,facecolor='white',edgecolor='#565656',lw=.25,zorder=0))
        z=d[d.subperiod==t]
        for f in facilities.itertuples():
            x,y=albers(f.LON,f.LAT)
            for form in ['existing','new']:
                zz=z[(z.location==f.location)&(z.form==form)]
                cap=zz.installed_capacity.sum();prod=zz.production.sum()
                if cap<1e-6:continue
                tech=zz.loc[zz.installed_capacity.idxmax(),'tech']
                hatch='////' if form=='new' else None
                ax.scatter([x],[y],s=cap*.72,facecolors=COLORS[tech],edgecolors='#333333',lw=.35,
                           hatch=hatch,zorder=3,alpha=.55 if prod<cap*.985 else 1.)
                if prod<cap*.985:
                    ax.scatter([x],[y],s=prod*.72,facecolors=COLORS[tech],edgecolors='white',lw=.4,hatch=hatch,zorder=4)
                if t==1:ax.annotate(str(f.location),(x,y),xytext=(4,-6),textcoords='offset points',fontsize=5,zorder=8)
            if matrices[1][f.location-1,t-1]:ax.scatter([x],[y],s=8,marker='s',c='#ce1b29',edgecolors='white',lw=.2,zorder=7)
        ax.set(xlim=xlim,ylim=ylim);ax.set_aspect('equal');ax.axis('off')
        ax.set_title(f'Subperiod {t}',fontsize=7,pad=1)
        ax.text(.0,1.02,string.ascii_lowercase[t+3],transform=ax.transAxes,fontweight='bold',fontsize=8.5)
    tech_legend(fig.add_axes([.08,.027,.40,.20]),size=6.2)
    ax=fig.add_axes([.48,.045,.22,.18]);ax.axis('off')
    ax.legend(handles=[Line2D([0],[0],marker='o',color='black',ls='',markersize=math.sqrt(v*.72),label=f'{v} t yr$^{{-1}}$') for v in [50,100,150,200]],
        title='Installed capacity',title_fontsize=6.6,loc='upper left',borderaxespad=0,labelspacing=1.1)
    ax=fig.add_axes([.71,.045,.28,.18]);ax.axis('off')
    ax.legend(handles=[Patch(facecolor='white',edgecolor='black',lw=.3,label='Existing facility'),
        Patch(facecolor='white',edgecolor='black',lw=.3,hatch='////',label='New-build facility'),
        Line2D([0],[0],marker='s',color='#ce1b29',markersize=3,ls='',label='Retired facility'),
        Line2D([0],[0],marker='o',color='#888',mfc='white',ms=5,ls='',label='Inner circle = production')],
        title='Facility type',title_fontsize=6.6,loc='upper left',borderaxespad=0,labelspacing=.9)
    export(fig,'Figure_4_facility_decisions_maps')


def tradeoff_data():
    bau=load('bau_nominal')[0];points=[]
    for i in range(10):
        s=load(f'target_{i:02d}')[0]
        points.append(dict(target=s['target'],npv=s.get('npv',s['objective']),emissions=s['total_co2'],
            gap=s['mip_gap'],status=s['status']))
    p=pd.DataFrame(points).sort_values('target')
    p['incremental_cost']=p.npv-bau.get('npv',bau['objective'])
    p['emissions_avoided']=bau['total_co2']-p.emissions
    prev_c=np.r_[bau.get('npv',bau['objective']),p.npv.to_numpy()[:-1]]
    prev_e=np.r_[bau['total_co2'],p.emissions.to_numpy()[:-1]]
    denom=prev_e-p.emissions.to_numpy()
    p['marginal_cost_MUSD_per_t']=np.divide(p.npv-prev_c,denom,out=np.full(len(p),np.nan),where=denom>1e-7)
    p.to_csv(DATA/'fig5_tradeoff.csv',index=False)
    return bau,p


def fig5():
    setup();bau,p=tradeoff_data()
    fig,axes=plt.subplots(1,3,figsize=(7.204724,2.992126))
    fig.subplots_adjust(left=.09,right=.99,bottom=.23,top=.91,wspace=.50)
    for ax,label in zip(axes,'abc'):
        ax.text(-.2,1.08,label,transform=ax.transAxes,fontweight='bold',fontsize=8.5)
        ax.yaxis.set_major_locator(MaxNLocator(5));ax.tick_params(labelsize=6.5)
    ax=axes[0];ax.plot(p.emissions,p.npv,'o-',color='#b32880',lw=.8,ms=3)
    ax.scatter([bau['total_co2']],[bau.get('npv',bau['objective'])],color='#a73086',marker='*',s=24)
    ax.annotate('BAU',(bau['total_co2'],bau.get('npv',bau['objective'])),xytext=(-17,6),textcoords='offset points',fontsize=6)
    ax.set(xlabel=r'Cumulative emissions (t CO$_2$)',ylabel='Net present cost (MUSD)')
    for i,offset in [(0,(-14,-14)),(9,(3,-14))]:
        ax.annotate(f'{p.target.iloc[i]:.0%}',(p.emissions.iloc[i],p.npv.iloc[i]),xytext=offset,
            textcoords='offset points',fontsize=6.5,arrowprops={'arrowstyle':'-','lw':.4})
    ax=axes[1];ax.plot(p.emissions_avoided,p.incremental_cost,'o-',color='#004c70',lw=.8,ms=3)
    ax.set(xlabel=r'Emissions avoided (t CO$_2$)',ylabel='Incremental cost (MUSD)')
    for i,offset in [(0,(3,6)),(9,(-17,-14))]:
        ax.annotate(f'{p.target.iloc[i]:.0%}',(p.emissions_avoided.iloc[i],p.incremental_cost.iloc[i]),xytext=offset,textcoords='offset points',fontsize=6.5)
    ax=axes[2]
    ax.plot(p.target*100,p.marginal_cost_MUSD_per_t,'o-',color='#008647',lw=.9,ms=3)
    ax.set(xlabel='Emissions target (%)',ylabel=r'Marginal cost (MUSD t$^{-1}$ CO$_2$)',xlim=(-2,52))
    ax.set_xticks([0,10,20,30,40,50]);ax.set_ylim(bottom=0)
    undefined=p.marginal_cost_MUSD_per_t.isna()
    if undefined.any():
        ax.scatter(p.loc[undefined,'target']*100,np.full(undefined.sum(),.025),marker='x',s=10,
                   color='#777',transform=ax.get_xaxis_transform(),clip_on=False,zorder=7)
    # The bands follow the article's visual structure; no phase boundaries are inferred.
    for lo,hi in [(0,10),(20,30),(40,50)]:ax.axvspan(lo,hi,color='#dedede',alpha=.55,zorder=0)
    export(fig,'Figure_5_cost_emission_tradeoff')


def fig6():
    setup();fig=plt.figure(figsize=(7.204724,5.472441))
    gs=fig.add_gridspec(2,5,width_ratios=[1,1,1,1,1.1],left=.065,right=.90,bottom=.11,top=.955,hspace=.40,wspace=.40)
    summaries=[]
    for i in range(8):
        n=i+1;s,d,_,_,_=load(f'scenario_{n}');b=load(f'scenario_{n}_bau')[0]
        ax=fig.add_subplot(gs[i//4,i%4]);panel(ax,string.ascii_lowercase[i]);y=stack_tech(ax,d,'installed_capacity')
        ax.set_title(f'Scenario {n}',pad=3);ax.set_ylim(0,y.max()*1.3)
        if i//4==1:ax.set_xlabel('Subperiod',labelpad=1)
        if i%4==0:ax.set_ylabel(r't yr$^{-1}$',labelpad=1)
        for row,(label,v,high) in enumerate([('Demand',s['demand_growth'],.02),('Elec. supply',s['electricity_growth'],.05),('CCS cost',s['ccs_cost_factor'],2)]):
            ay=.955-row*.067
            ax.text(.03,ay,label,transform=ax.transAxes,fontsize=5.2,va='top')
            color='#dc4848' if v==high else '#4b9c74'
            ax.add_patch(Rectangle((.78,ay-.046),.09,.045,transform=ax.transAxes,facecolor=color,edgecolor='none'))
            ax.text(.91,ay-.022,'+' if v==high else '-',transform=ax.transAxes,ha='center',va='center',fontsize=5)
        cost=s.get('npv',s['objective'])-b.get('npv',b['objective']);avoid=b['total_co2']-s['total_co2']
        summaries.append(dict(scenario=n,incremental_cost=cost,emissions_avoided=avoid,
            average_cost=cost/avoid if avoid>1e-7 else np.nan,gap=s['mip_gap'],bau_gap=b['mip_gap']))
        d.assign(scenario=n).to_csv(DATA/f'fig6_scenario_{n}_capacity.csv',index=False)
    leg=fig.add_subplot(gs[0,4]);tech_legend(leg,size=5.4)
    tech_handles=leg.get_legend();leg.add_artist(tech_handles)
    leg.legend(handles=[Patch(facecolor='white',edgecolor='black',lw=.3,label='Existing facility'),
        Patch(facecolor='white',edgecolor='black',lw=.3,hatch='////',label='New-build facility')],
        title='Facility type',title_fontsize=6,fontsize=5.4,loc='lower left',borderaxespad=0,handlelength=1.1)
    metrics=pd.DataFrame(summaries);metrics.to_csv(DATA/'fig6_sensitivity_metrics.csv',index=False)
    ax=fig.add_subplot(gs[1,4]);ax2=ax.twinx();ax3=ax.twinx();ax3.spines['right'].set_position(('axes',1.53))
    ax.text(-.20,1.07,'i',transform=ax.transAxes,fontweight='bold',fontsize=8.5)
    ax.bar(metrics.scenario-.16,metrics.incremental_cost,width=.28,color='#168e91')
    ax2.bar(metrics.scenario+.16,metrics.emissions_avoided,width=.28,color='#a7449f')
    ax3.plot(metrics.scenario,metrics.average_cost,ls='',marker='v',ms=2.4,mec='#555',mfc='white',mew=.5)
    ax.set_xticks(range(1,9));ax.tick_params(labelsize=5);ax.set_xlabel('Scenario',labelpad=1)
    for a,color,label in [(ax,'#168e91','Incremental cost (MUSD)'),(ax2,'#a7449f',r'Emissions avoided (t CO$_2$)'),(ax3,'#555',r'Average cost (MUSD t$^{-1}$ CO$_2$)')]:
        a.set_ylabel(label,color=color,fontsize=5.4,labelpad=1);a.tick_params(axis='y',labelsize=5,colors=color,pad=1)
        a.yaxis.set_major_locator(MaxNLocator(4));a.set_ylim(bottom=0)
        fmt=ScalarFormatter(useMathText=True);fmt.set_powerlimits((0,0));a.yaxis.set_major_formatter(fmt);a.yaxis.get_offset_text().set_size(5)
    export(fig,'Figure_6_factorial_sensitivity')


def main():
    global RESULTS, FIG, DATA
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('--figures',nargs='+',type=int,default=[3,4,5,6])
    parser.add_argument('--results',type=Path,default=ROOT)
    args=parser.parse_args()
    RESULTS=args.results.resolve()
    FIG=RESULTS/'figures'
    DATA=RESULTS/'source_data'
    for n in args.figures:
        globals()[f'fig{n}']();print(f'Exported Figure {n}: SVG/PDF/PNG/TIFF',flush=True)


if __name__=='__main__':main()
