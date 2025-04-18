#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:29:38 2022

@author: gao
"""

"""
f: fraction of germ cells
p1: p_{g->s}
p2: p_{s->g}
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D

#-------------------------------------------------------------------------------
# set display width
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)

#-------------------------------------------------------------------------------
# fraction of germ-role with respect to n
def fra(n,p1,p2):
	return ((1-p1-p2)**n*p1+p2)/(p1+p2)

# growth time 
def time(alpha,beta,b,c,p1,p2,n):	
	t=0
	for i in range(1,n+1):
		f=fra(i-1,p1,p2)
		t_i=(1+c*(f*p1+beta*(1-f)*p2))/(1+b*(1-f)**alpha)
		t=t+t_i
	return t

# growth rate
def grate(alpha,beta,b,c,p1,p2,n):
	f=fra(n,p1,p2)
	dno=np.log((2**n)*f)
	numo=time(alpha,beta,b,c,p1,p2,n)
	return (dno)/(numo)

#----figure for changing p2----------------------------
from matplotlib.colors import LinearSegmentedColormap

import matplotlib as mpl
import numpy as np

def colorFader(c1,c2,mix=0): #fade (linear interpolate) from color c1 (at mix=0) to c2 (mix=1)
    c1=np.array(mpl.colors.to_rgb(c1))
    c2=np.array(mpl.colors.to_rgb(c2))
    return mpl.colors.to_hex((1-mix)*c1 + mix*c2)

def gra_c_list(c1,c2,c3,num_c):
	
	color_list=[]
	for x in range(num_c+1):
		color_list.append(colorFader(c1,c2,x/num_c))
		
	for x in range(num_c+1):
		color_list.append(colorFader(c2,c3,x/num_c))
	return color_list

col1=gra_c_list('#e0e0e0','#bababa','#4d4d4d',101) # grey
col2=gra_c_list('#d1e5f0','#4393c3','#2166ac',101) # blue
#col1=gra_c_list('#fddbc7','#d6604d','#67001f',101) #redd

color1=col1[::-1]+col2
#color1=col1=gra_c_list('#bababa','#4393c3','#2166ac',101)
					 
'''Parameters'''
n=200
alpha=beta=b=c=1

fig, ax = plt.subplots(1, 2, gridspec_kw={'width_ratios': [5,0.1]},figsize=(5.5, 4.7))
#plt.title("alpha=%s; beta=%s; b=%s; c=%s; n=%s"%(alpha,beta,b,c,n) )
start=8
x=np.linspace(1,n,n)
x1=np.linspace(-15,n+16,n+16+start+7)
x2=np.linspace(1,n,3)

cmap1 = LinearSegmentedColormap.from_list('mycmap', color1)
cmap2 = LinearSegmentedColormap.from_list('mycmap', color1)
c1 = np.arange(1,11)
norm1 = plt.Normalize(vmin=c1.min(), vmax=c1.max())
cmap1 = mpl.cm.ScalarMappable(norm=norm1, cmap=cmap1)
cmap1.set_array([])

"boundries"
#----------------- bound of growth rate; now it is based on isd
"isd and rd data"
sets=10
scale=1/sets

#--- grate for isd
rate_max_isd=[]
for j in range(1,sets):
	p1=scale*j 
	y_isd=[]
	for i in range(1,n+1):
		y_isd.append(grate(alpha,beta,b,c,p1,0,i))
	rate_max_isd.append(y_isd[-1])
	ax[0].plot(x, y_isd, ':',label="",c=cmap1.to_rgba( j+ 1),linewidth=2,alpha=0.8)
	rate_max=[]
	for k in range(1,2):
		p2=scale*k 
		y_rd=[]
		for i in range(1,n+1):
			y_rd.append(grate(alpha,beta,b,c,p1,p2,i))
		rate_max.append(y_rd[-1])
		# red dashes, blue squares and green triangles
		ax[0].plot(x, y_rd,'-',c=cmap1.to_rgba(j+1),linewidth=1,alpha=1)

# interval should be 0.1
M_grate=(1+b)*np.log(2*(1-scale))             # p_g->g=0.9 
m_grate=(1+b)*np.log(2*scale)             # but for p_{g->g}<0.5

# boundaries of isd
n_posi=170
alpha_bond=0.4
linewidth0=1.2
ax[0].plot(x1, [m_grate for i in x1], '-',c='k',linewidth=linewidth0,alpha=alpha_bond)			 
# lower boundary of rd and nd
h=0.1
grate_rd_mini=np.log(4)/(2+c*(1+beta))
ax[0].plot(x1, [M_grate for i in x1], '-',c='k',linewidth=linewidth0,alpha=alpha_bond)
ax[0].plot(x1, [grate_rd_mini for i in x1], '-',c='k',linewidth=linewidth0,alpha=alpha_bond)
ax[0].plot(x, [np.log(2) for i in x], '-.',c='k',linewidth=linewidth0*0.5,alpha=1)

ax[0].text(100,np.log(2)-0.1,r"$\lambda_{ND}$",color='k',fontsize=12,alpha=1)
ax[0].text(-12,M_grate+h,r"max$\lambda_{RD_{R}^{(S)}}$",fontsize=12,alpha=alpha_bond)			 
ax[0].text(-12,m_grate+h,r"min$\lambda_{RD_{R}^{(S)}}$",fontsize=12,alpha=alpha_bond)
ax[0].text(156,M_grate+h,r"max$\lambda_{RD_{RS}^{(S)}}$",fontsize=12,alpha=alpha_bond)			 
ax[0].text(159,grate_rd_mini+h,r"min$\lambda_{RD_{RS}^{(S)}}$",fontsize=12,alpha=alpha_bond)

# frame color
color_f="darkgray"
ax[0].spines['bottom'].set_color(color_f)
ax[0].spines['top'].set_color(color_f)
ax[0].spines['right'].set_color(color_f)
ax[0].spines['left'].set_color(color_f)
ax[0].tick_params(direction='in', length=2, width=1, colors='k')     # inout

ax[0].set_ylabel(r'Reproduction rate, $\lambda(n)$',fontsize=16)
ax[0].set_xlabel(r'Organismal size, $2^n$',fontsize=16)

#ax[0] = plt.gca()
ax[0].set_xlim([x[0]-2*start, x[-1]+2*start])
ax[0].set_ylim([m_grate-0.3, M_grate+0.35])

import matplotlib.lines as mlines
nd_line = mlines.Line2D([], [], color='k', linestyle='-.',label=r'$ND$')
rd_line = mlines.Line2D([], [], color='k', linestyle='-',label=r'$RD_{RS}^{(S)}$')
isd_line = mlines.Line2D([], [], color='k', linestyle=':',label=r'$RD_{R}^{(S)}$')

pos = ax[0].get_position()
ax[0].set_position([pos.x0, pos.y0, pos.width, pos.height])
#ax[0].legend(handles=[nd_line,rd_line,isd_line],loc='upper center', 
#		bbox_to_anchor=(0.5, 1.09),ncol=3,frameon=False,fontsize=11.5)
ax[0].legend(handles=[nd_line,rd_line,isd_line],loc='upper center', 
		bbox_to_anchor=(0.86, 0.35),ncol=1,frameon=False,fontsize=10)

#  color bar
img = ax[1].imshow(np.array([[0,1]]), cmap=cmap2,alpha=0.8)
img.set_visible(False)
ax[1].axis('off')
cbar_ax0 = fig.add_axes([0.85, 0.15, 0.016, 0.7])
cbar0 =fig.colorbar(img, cax=cbar_ax0,orientation="vertical")
cbar0.ax.tick_params(labelsize=8,length=2,direction='in')
cbar0.outline.set_visible(False)
cbar0.ax.set_ylabel(r'Differentiation probabilities, $p_{S \to R}=0.1, p_{R \to S}$', rotation=90,fontsize=12)

plt.show()
fig.savefig('./fig4_bounds_static.pdf', bbox_inches = 'tight')   # save figures


