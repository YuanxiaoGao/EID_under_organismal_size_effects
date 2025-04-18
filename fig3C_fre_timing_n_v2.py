#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 14:40:40 2018
@author: gao
"""
import os
import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
#-------------------------------------------------------------------------------
# set display width
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)
#-------------------------------------------------------------------------------
'''Parameters'''
alpha=1
beta=1
n_cluster=30
num_dup=100
bc_points=[1,10]
#-------------------------------------------------------------------------------
'''read data '''
grid_num=2
result=[[[[[]for n in range(2,n_cluster)]for c in range(2)] for b in range(2)] for i in range(4)] 
#result[stra_cluster][b][c][n-2][dup]
# save the four categories ; stra_cluster    0=RD; 1=IGD, 2=ISD, 3=IGSD

for dup in range(num_dup):  
	for n in range(2,n_cluster):  
		for b in range(0,grid_num):  
			for c in range(0,grid_num): 
				for stra_cluster in range(4): 
					if os.path.exists('./data_bc110/%s_%s_%s_%s_%s.txt'%(dup,n,b,c,stra_cluster)):
						data=np.loadtxt('./data_bc110/%s_%s_%s_%s_%s.txt'%(dup,n,b,c,stra_cluster))						
						result[stra_cluster][b][c][n-2].append(data)
					else:
						print('%s_%s_%s_%s_%s'%(dup,n,b,c,stra_cluster))

##-------------DRAW FIGURE A B C D----------------------------------------------------------------
import matplotlib.pyplot as plt
mark_list=['8','^','8','x']
label_list=['RD','IGD','ISD','IGSD']
color_list=['#238b45','#2171b5','#984ea3','#ff7f00','#cb181d']
size=[i for i in range(2,n_cluster)]

#================color list============================
import matplotlib as mpl
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import FormatStrFormatter

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

#color2=gra_c_list('w','#4292c6','#08306b',31)[5:]
color0=gra_c_list('w','#bf812d','#543005',101)[10:]     # brown
cmap0 = LinearSegmentedColormap.from_list('mycmap', color0)
color00=gra_c_list('w','#4292c6','#08306b',101)[10:]
cmap00 = LinearSegmentedColormap.from_list('mycmap', color00)


#gra_c_list('w','#bf812d','#543005',101)[60]
#gra_c_list('w','#4292c6','#08306b',101)[60]			  
# ====== set fig frame ========
fig, ax = plt.subplots(2, 2, gridspec_kw={'width_ratios': [1,1],'height_ratios':[1,1],
														'wspace':0.01, 'hspace':0.01},figsize=(7,3.5))
	
for i in range(2):
	for j in range(2):
		ax[i][j].spines['bottom'].set_color('gray')
		ax[i][j].spines['top'].set_color('gray') 
		ax[i][j].spines['right'].set_color('gray')
		ax[i][j].spines['left'].set_color('gray')
		ax[i][j].tick_params(axis='both', colors='k', labelsize=8)
		ax[i][j].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))		
ax[0][0].xaxis.set_ticklabels([])
ax[0][1].xaxis.set_ticklabels([])
ax[0][1].set_yticklabels([])
ax[1][1].set_yticklabels([])

fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.ylabel(r'Mean fraction of cell division times' '\n' 'when cell irreversibility occurs',fontsize=11,labelpad=-17)
#-----colormap--------------
size=[i for i in range(2,n_cluster)]
alpha0=0.85
markers=5

for b in range(0,2):  
	for c in range(0,2):    #grid_num
		g_mean=[]
		g_std=[]
		igsd_g=[[]for i in range(2,n_cluster)]
		s_mean=[]
		s_std=[]
		igsd_s=[[]for i in range(2,n_cluster)]
		for stra_cluster in [3]: 
			for n in range(2,n_cluster):
				for dup in range(num_dup):	
					item1=result[stra_cluster][b][c][n-2][dup][:-1,:]										
					deep=[]
					deep_s=[]
					for i in range(1,n+1):
						if all(x==1 for x in item1[-i:,5]): 
							deep.append((n+1-i)/n)
						if all(x==1 for x in item1[-i:,0]):
							deep_s.append((n+1-i)/n)
					igsd_g[n-2].append(np.min(np.array(deep)))
					igsd_s[n-2].append(np.min(np.array(deep_s)))
				g_mean.append(np.mean(np.array(igsd_g[n-2])))
				g_std.append(np.std(np.array(igsd_g[n-2])))
				s_mean.append(np.mean(np.array(igsd_s[n-2])))
				s_std.append(np.std(np.array(igsd_s[n-2])))
		c1 = np.arange(1, 101)
		norm1 = mpl.colors.Normalize(vmin=c1.min(), vmax=c1.max())
		cmap1 = mpl.cm.ScalarMappable(norm=norm1, cmap=cmap0)
		cmap1.set_array([])
		
		cmap2 = mpl.cm.ScalarMappable(norm=norm1, cmap=cmap00)
		cmap2.set_array([])
		
		color_list=[]
		color_list1=[]
		for k in range(len(g_mean)):
			co=g_mean[k]*100
			color_list.append(cmap1.to_rgba(co))
			co1=s_mean[k]*100
			color_list1.append(cmap2.to_rgba(co1))
	
		ax[b,c].plot(size,g_mean,linestyle='-',linewidth=1,color="k",alpha=0.6)
		ax[b,c].plot(size,s_mean,linestyle='-',linewidth=1,color="k",alpha=0.6)
		
		ax[b,c].tick_params(direction='in', length=1, width=1, colors='k')
		ax[b,c].tick_params(axis='both', which='major', labelsize=8)
		for l in range(len(g_mean)):
			ax[b,c].errorbar(size[l],g_mean[l],g_std[l],linestyle='-',linewidth=0.3,markersize=markers,marker="h",mfc=color_list[l], color=color_list[l],alpha=1,zorder=5)
			eb2=ax[b,c].errorbar(size[l],s_mean[l],s_std[l],linestyle='-',errorevery=2,linewidth=1,markersize=markers,marker="h",mfc=color_list1[l], color=color_list1[l],alpha=1,zorder=5)
			eb2[-1][0].set_linestyle('dotted')
#		ax[b,c].set_title(r"$b=%s$, $c=%s$"%(bc_points[int(b)],bc_points[int(c)]),fontsize=10)
		ax[b,c].text(22,0.1,r"$b=%s$, $c=%s$"%(bc_points[int(b)],bc_points[int(c)]),fontsize=8)
		# lim 
		ax[b,c].set_ylim(0,1.17)
		ax[b,c].set_xlim(0.8,30.5)
		if b==1:
			ax[b,c].set_xlabel(r'Organismal size, $2^n$',fontsize=11,labelpad=4)
			
col1=gra_c_list('w','#bf812d','#543005',101)[150]
col2=gra_c_list('w','#4292c6','#08306b',101)[150]

#from matplotlib.lines import Line2D
sc1 = plt.scatter([],[],s=20,facecolors=col1, edgecolors=col1)
sc2 = plt.scatter([],[],s=20,facecolors=col2, edgecolors=col2)
ax[1][1].legend([sc1, sc2], ["R cell", "S cell"],fontsize=8,loc="lower left",frameon=False, ncol=2)	


		##------x and y labels and limits--------
plt.tight_layout()
plt.show()
#fig.savefig('./fig3C_irre-stage-b%s_c%s.pdf'%(bc_points[b],bc_points[c]), bbox_inches = 'tight')   # save figures


#fig, ax = plt.subplots(1,1,figsize=(4,4))
#
#norm1 = plt.Normalize(0, 1)
#img = ax.imshow(np.array([[0,1]]), cmap=cmap0,alpha=1)
#img.set_visible(False)
#ax.axis('off')
#cbar_ax0 = fig.add_axes([0.2, 0.3, 0.28, 0.02])
#cbar0 =fig.colorbar(img, cax=cbar_ax0,orientation="horizontal",norm=norm1, boundaries=None)
#cbar0.ax.tick_params(labelsize=8,length=2,direction='in')
#cbar0.outline.set_visible(False)
#cbar0.ax.set_xlabel(r'Fraction of R cell', rotation=0,fontsize=9)
#cbar0.outline.set_visible(False)
#cbar0.ax.xaxis.set_label_position('top')
#cbar0.ax.xaxis.set_major_locator(plt.MaxNLocator(1))
##cbar0.ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
#plt.show()
#fig.savefig('./legand_r.pdf')
#
#fig, ax = plt.subplots(1,1,figsize=(4,4))
#
#norm1 = plt.Normalize(0, 1)
#img = ax.imshow(np.array([[0,1]]), cmap=cmap00,alpha=1)
#img.set_visible(False)
#ax.axis('off')
#cbar_ax0 = fig.add_axes([0.2, 0.3, 0.28, 0.02])
#cbar0 =fig.colorbar(img, cax=cbar_ax0,orientation="horizontal",norm=norm1, boundaries=None)
#cbar0.ax.tick_params(labelsize=8,length=2,direction='in')
#cbar0.outline.set_visible(False)
#cbar0.ax.set_xlabel(r'Fraction of S cell', rotation=0,fontsize=9)
#cbar0.outline.set_visible(False)
#cbar0.ax.xaxis.set_label_position('top')
#cbar0.ax.xaxis.set_major_locator(plt.MaxNLocator(1))
#plt.show()
#fig.savefig('./legand_s.pdf')


