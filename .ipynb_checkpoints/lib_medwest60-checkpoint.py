#!/usr/bin/env  python
#=======================================================================
"""
StepHANIE Leroux
Collection of my "customed" tools related to  MEDWEST60 analysis...
"""
#=======================================================================


## standart libraries
import os,sys
import numpy as np

from scipy.signal import argrelmax
from scipy.stats import linregress

# xarray
import xarray as xr

# plot
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.colors import Colormap
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
import matplotlib.cm as cm
import matplotlib.dates as mdates
import matplotlib.ticker as mticker
from matplotlib.colors import from_levels_and_colors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmocean


# custom tools
import lib_medwest60

def readmedwestens_1mb(diriprefix='/Users/leroux/DATA/MEDWEST60_DATA/',
                       CONFIGCASEmed='MEDWEST60-GSL04',
                       ens='ens01',
                       mb='001',
                       CONFIGCASEref='MEDWEST60-BLBT02',
                       typ="gridT-2D",
                       varna='sosstsst',diriprefixref="/Users/leroux/DATA/MEDWEST60_DATA/extracted_eNATL60/1h/",
                       bathyfilepath='/Users/leroux/DATA/MEDWEST60_DATA/MEDWEST60-I/MEDWEST60_Bathymetry_v3.3.nc4'):

    
    '''
    Goal: read data from one ensemble member.
    Take parameters:
    diriprefix='/Users/leroux/DATA/MEDWEST60_DATA/',
    CONFIGCASEmed='MEDWEST60-GSL04',
    ens='ens01',
    mb='001',
    CONFIGCASEref='MEDWEST60-BLBT02',
    typ="gridT-2D",
    varna='sosstsst',diriprefixref="/Users/leroux/DATA/MEDWEST60_DATA/extracted_eNATL60/1h/",
    bathyfilepath='/Users/leroux/DATA/MEDWEST60_DATA/MEDWEST60-I/MEDWEST60_Bathymetry_v3.3.nc4'
    
    Returns:  nav_lat,nav_lon,bathy,data,varname,latexvarname
    '''

    dirimed = diriprefix+CONFIGCASEmed+"-S/"+ens+"/1h/"+typ+"/"
    print(dirimed)
    #diriref = diriprefixref+typ+"/"+"tmp/"  

    filiprefix = mb+CONFIGCASEmed+"-"+ens+"_1h_*"+typ+"_"
    print(filiprefix+"*.nc")
   
    bathy =  xr.open_dataset(bathyfilepath)["Bathymetry"]
    # longitude
    nav_lon = xr.open_dataset(bathyfilepath)['nav_lon']
    # latitude
    nav_lat = xr.open_dataset(bathyfilepath)['nav_lat']
    
    data   = xr.open_mfdataset(dirimed+filiprefix+"*.nc",concat_dim='time_counter',decode_times=True)[varna]
        
    varname,latexvarname = latexvaname(varna)
    
    return nav_lat,nav_lon,bathy,data,varname,latexvarname


def readallmbs(NMBtot=1,typ="gridT-2D", varna='sossheig',CONFIGCASEmed='MEDWEST60-GSL04',ens='ens01',CONFIGCASEref='MEDWEST60-GSL04'):

    ie=1

    for ie in range(1,NMBtot+1):

        if ie<10:
            mbn="00"+str(ie)
        if ie>9:
            mbn="0"+str(ie)

        nav_lat,nav_lon,bathy,tmpMB,varname,latexvarname = readmedwestens_1mb(typ=typ, mb=mbn,varna=varna,CONFIGCASEmed=CONFIGCASEmed,ens=ens,CONFIGCASEref=CONFIGCASEref)

        if (ie!=1):
            concdata = xr.concat([concdata,tmpMB], dim='e')
        else:
            concdata = tmpMB
    return nav_lat,nav_lon,bathy,concdata,varname,latexvarname 

def readmedwesttwinexp(diriprefix='/Users/leroux/DATA/MEDWEST60_DATA/',
                       CONFIGCASEmed='MEDWEST60-GSL03',
                       CONFIGCASEref='MEDWEST60-BLBT02',
                       typ="gridT-2D",
                       varna='sosstsst',diriprefixref="/Users/leroux/DATA/MEDWEST60_DATA/extracted_eNATL60/1h/",
                       bathyfilepath='/Users/leroux/DATA/MEDWEST60_DATA/MEDWEST60-I/MEDWEST60_Bathymetry_v3.3.nc4'):

    dirimed_el = diriprefix+CONFIGCASEmed+"-S/1h/"+typ+"/"
    dirimed_lf = diriprefix+CONFIGCASEmed+"bis-S/1h/"+typ+"/"
    #diriref = diriprefixref+typ+"/"+"tmp/"  

    if typ=='curloverf':
        typ2='curloverF'
        filiprefix_el = CONFIGCASEmed+"_1h_20100???_20100???_"+typ2
        filiprefix_lf = CONFIGCASEmed+"bis_1h_20100???_20100???_"+typ2
    else:
        filiprefix_el = CONFIGCASEmed+"_1h_20100???_20100???_"+typ
        filiprefix_lf = CONFIGCASEmed+"bis_1h_20100???_20100???_"+typ

    
    bathy =  xr.open_dataset(bathyfilepath)["Bathymetry"]
    # longitude
    nav_lon = xr.open_dataset(bathyfilepath)['nav_lon']
    # latitude
    nav_lat = xr.open_dataset(bathyfilepath)['nav_lat']
    
    data_el   = xr.open_mfdataset(dirimed_el+filiprefix_el+"*_merg.nc",concat_dim='time_counter',decode_times=True)[varna]
    data_lf   = xr.open_mfdataset(dirimed_lf+filiprefix_lf+"*_merg.nc",concat_dim='time_counter',decode_times=True)[varna]

    if (data_el.time_counter[0]==data_lf.time_counter[0]):
        print("Dates ok. READ : "+varna)
        diff = data_el - data_lf
    else:
        print("!! dates not ok !! CANT COMPUTE DIFF")
        diff = 0
        
    if varna=='sosstsst':
        varname='SST'
        latexvarname=varname
    if varna=='sossheig':
        varname='SSH'
        latexvarname=varname
    if varna=='socurloverf':
        varname='curloverf'
        latexvarname="$\zeta/f$"
    
    return nav_lat,nav_lon,bathy,data_el,data_lf,diff,varname,latexvarname



def readonlybathy(bathyfilepath='/Users/leroux/DATA/MEDWEST60_DATA/MEDWEST60-I/MEDWEST60_Bathymetry_v3.3.nc4'):

    bathy =  xr.open_dataset(bathyfilepath)["Bathymetry"]
    # longitude
    nav_lon = xr.open_dataset(bathyfilepath)['nav_lon']
    # latitude
    nav_lat = xr.open_dataset(bathyfilepath)['nav_lat']
    
    return nav_lat,nav_lon,bathy
    
    
    
    

def plotmapMEDWEST(fig1,ehonan,nav_lon,nav_lat,cmap,norm,plto='tmp_plot',typlo='pcolormesh',coastL=False,coastC=False,coastLand=False,xlim=(0,10), ylim=(0,10),su='b',so='k',loncentr=0.,latcentr=0.,labelplt="",incrgridlon=5,incrgridlat=5,edgcol1='#585858',edgcol2='w',mk="o",mks=0.1,scattcmap=True,scattco='k'):

        # Projection
        trdata  = ccrs.PlateCarree() 
        
        ax = plt.axes(projection= ccrs.PlateCarree())
        ax.outline_patch.set_edgecolor(edgcol2)

        gridl=True
        if gridl:
            # gridlines
            gl = ax.gridlines(draw_labels=True,linewidth=1, color='#585858', alpha=0.2, linestyle='--')
            # grid labels
            #label_style = {'size': 12, 'color': 'black', 'weight': 'bold'}
            label_style = {'size': 12, 'color': '#BDBDBD', 'weight': 'normal'}
            
            gl.xlabel_style = label_style
            gl.xlabels_bottom = False
            gl.xlocator = mticker.FixedLocator(np.arange(-180,180,incrgridlon,dtype=float))
            gl.ylabel_style = label_style
            gl.ylabels_right = False
            gl.ylocator = mticker.FixedLocator(np.arange(-90,90,incrgridlat,dtype=float))

        # Add Coastlines and or plain continents
        if coastC:
            ax.add_feature(ccf.COASTLINE, facecolor='w', edgecolor='none')
        if coastLand:
            ax.add_feature(ccf.LAND, facecolor='w', edgecolor='none')
        if coastL:
            #ax.coastlines(color='#585858',linewidth=1)
             ax.coastlines(color='w',linewidth=1)
        
        ### PLOTS: 
        if typlo=='pcolormesh':
            cs  = plt.pcolormesh(nav_lon, nav_lat, ehonan,cmap=cmap,transform=trdata,norm=norm)
        
        if typlo=='contourf':
            cs  = plt.contourf(nav_lon, nav_lat, ehonan,transform=trdata,levels=levels,norm=norm,cmap=cmap,extend='both')

        # geographical limits
        plt.xlim(xlim)
        plt.ylim(ylim) 
        
        
        
def plotmapMEDWEST_gp(fig3,ax,data2plot,cmap,norm,plto='tmp_plot',gridpts=True,gridptsgrid=False,gridinc=200,gstyle='lightstyle'): 
    
    cs  = ax.pcolormesh(data2plot,cmap=cmap,norm=norm)

    #ax = plt.gca()
    # Remove the plot frame lines. 
    ax.spines["top"].set_visible(False)  
    ax.spines["bottom"].set_visible(False)  
    ax.spines["right"].set_visible(False)  
    ax.spines["left"].set_visible(False)  

    ax.tick_params(axis="both", which="both", bottom="False", top="False",  
                labelbottom="False", labeltop='False',left="False", right="False", labelright="False",labelleft="False")  

    
    if gridpts:
    # show gridpoint on axes
        ax.tick_params(axis="both", which="both", bottom="False", top="False",  
                labelbottom="True", labeltop='False',left="False", right="False", labelright="False",labelleft="True")  
        plto = plto+"_wthgdpts"

    if gridptsgrid:
        lstylegrid=(0, (5, 5)) 
        if (gstyle=='darkstyle'):
            cmap.set_bad('#424242')
            lcolorgrid='w'#"#585858" # "#D8D8D8"
            tcolorgrid='#848484'#"#848484"
            
        if (gstyle=='ddarkstyle'):
            cmap.set_bad('#424242')
            lcolorgrid='w'#"#585858" # "#D8D8D8"
            tcolorgrid='w'#'#848484'#"#848484"
        if (gstyle=='lightstyle'):
            cmap.set_bad('w')
            lcolorgrid="#585858" # "#D8D8D8"
            tcolorgrid='#848484'#"#848484"            

        lalpha=0.2
        lwidthgrid=1.
        #ax = plt.gca()
        ax.xaxis.set_major_locator(mticker.MultipleLocator(gridinc))
        ax.yaxis.set_major_locator(mticker.MultipleLocator(gridinc))   
        ax.tick_params(axis='x', colors=tcolorgrid)
        ax.tick_params(axis='y', colors=tcolorgrid)
        ax.grid(which='major',linestyle=lstylegrid,color=lcolorgrid,alpha=lalpha,linewidth=lwidthgrid)
        ax.axhline(y=1.,xmin=0, xmax=883,zorder=10,color=lcolorgrid,linewidth=lwidthgrid,linestyle=lstylegrid,alpha=lalpha )
    
    return cs,ax



def plottsts(var1,var2,xl1=0,yl1=0,xl2=-1,yl2=-1,diro='./',namo='plotts.png',dpifig=300):
    if xl2==-1:
        xl2=xl1
    if yl2==-1:
        yl2=yl1
        
    ts1 = var1.isel(x=xl1,y=yl1)
    ts2 = var2.isel(x=xl2,y=yl2)

    fig2 = plt.figure(figsize=([15,8]),facecolor='white')
    
    plt.plot(ts1,color='blue',marker='o')
    plt.plot(ts2,color='k',marker='+')

    plt.show()
    saveplt(plt,fig2,diro,namo,dpifig=dpifig)
    return plt,fig2


def saveplt(fig,diro,namo,dpifig=300):
    fig.savefig(diro+namo, facecolor=fig.get_facecolor(), edgecolor='none',dpi=dpifig,bbox_inches='tight', pad_inches=0)
    plt.close(fig) 

def mycolormap(levbounds,cm_base='Spectral_r',cu='w',co='k'):
    lmin = levbounds[0]
    lmax = levbounds[1]
    incr = levbounds[2]
    levels = np.arange(lmin,lmax,incr)
    nice_cmap = plt.get_cmap(cm_base)
    colors = nice_cmap(np.linspace(0,1,len(levels)))[:]
    cmap, norm = from_levels_and_colors(levels, colors, extend='max')
    cmap.set_under(cu)
    cmap.set_over(co)
    return cmap,norm


def addcolorbar(fig,cs,ax,levbounds,levincr=1,tformat="%.2f",tlabel='',shrink=0.45,facmul=1.,orientation='vertical',pad=0.03,tc='k',loc='lower right'):
    lmin = levbounds[0]
    lmax = levbounds[1]
    incr = levbounds[2]
    levels = np.arange(lmin,lmax,incr)
    cblev = levels[::levincr]
    
    if orientation =='horizontal':
        axins1 = inset_axes(ax,
                        height="15%",  # height : 5%
                            width="50%",
                        loc=loc,
                        bbox_to_anchor=(0.08, 0.1,0.9,0.2),
                        bbox_transform=ax.transAxes,
                        borderpad=0)
        
    if orientation =='vertical':
        axins1 = inset_axes(ax,
                        height="50%",  # height : 5%
                            width="2%",
                        loc='center left',
                       borderpad=2)

    cb = fig.colorbar(cs,cax=axins1,pad=pad,
                                    extend='both',                   
                                    ticks=cblev,
                                    spacing='uniform',
                                    orientation=orientation,
                                    )
    
    new_tickslabels = [tformat % i for i in cblev*facmul]
    cb.set_ticklabels(new_tickslabels)
    cb.ax.set_xticklabels(new_tickslabels, rotation=70,size=10,color=tc)
    cb.ax.tick_params(labelsize=10,color=tc) 
    cb.set_label(tlabel,size=14,color=tc)
    
    
    return cb,axins1

def textunit(varname):
    if varname=='SST':
        suffix=" (ºC)"
    if varname=='SSH':
        suffix=" (m)"
    if varname=='curloverf':
        suffix=""
    return suffix

def textunitfac(varname,faclab):
    if ((faclab=="1")|(faclab=="")):
        if varname=='SST':
            suffix=" (ºC)"
        elif varname=='SSH':
            suffix=" (m)"
        elif varname=='curloverf':
            suffix=""
        elif ((varname=='e1t')|(varname=='e2t')):
            suffix=" (m)"
        else :
            suffix=""
    else :
        if varname=='SST':
            suffix=" ("+faclab+" ºC)"
        elif varname=='SSH':
            if faclab=='10$^{-3}$':
                suffix=" (mm)"       
            else:
                suffix=" ("+faclab+" m)"
        elif varname=='curloverf':
            suffix=" x("+faclab+")"
        else :
            suffix=" x("+faclab+")"    
    return suffix

def getslope(data,it1,it2,fac=1,tconv=24):
    # tconv conversion from data time freq to days
    from scipy.stats import linregress

    truc1=data*fac
    # regression linear on the log plot from it1 to it2
    resreg1=linregress(np.arange(it1,it2), np.log(truc1.isel(time_counter=slice(it1,it2))))

    # doubling time: td=ln(2)/k where k is the slope of ln(y) = ln(yo) + kt
    slope=(np.log(2)/resreg1.slope)/tconv  # in days
    #print(td1)
    return(slope)


def fxRMSEreg(diff,region=[0,0,0,0]):
    #reg=[x0,x1,y0,y1]; if unrequested then compute over all domain
    if np.array(region).sum()==0:
        region=[0,diff.shape[2]-1,0,diff.shape[1]-1]
    diffreg = diff.isel(x=slice(region[0],region[1]),y=slice(region[2],region[3]))
    diffsqreg = diffreg*diffreg
    diffsqstreg = diffsqreg.stack(z=('x', 'y'))
    RMSEreg = np.sqrt(diffsqstreg.mean(dim='z')).load()
    return RMSEreg


def fillnacorrection(dat):
    dat = dat.where(dat!=0.,-9999.)
    dat = dat.fillna(-9999.)
    dat = dat.where(dat!=-9999.)
    return dat


def stdens(varIN):
    STDens    = varIN.std(dim='e').load()
    STDens    = fillnacorrection(STDens)
    STDensdom = STDens.stack(z=('x', 'y'))
    STDensdom = STDensdom.mean(dim='z').load()
    return STDens,STDensdom

def latexvarname(varna):
    if varna=='sosstsst':
        varname='SST'
        latexvarname=varname
    elif varna=='sossheig':
        varname='SSH'
        latexvarname=varname
    elif varna=='socurloverf':
        varname='curloverf'
        latexvarname="$\zeta/f$"
    else :
        varname=varna
        latexvarname=varna
    return varname,latexvarname