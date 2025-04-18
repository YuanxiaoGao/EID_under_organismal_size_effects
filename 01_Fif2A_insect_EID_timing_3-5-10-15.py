#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  2 14:40:40 2018
@author: gao
"""
"""DP: six developement probabilities for cells
		[S->SS, S->GS, S->GG, G->SS, G->GS, G->GG] (np array: float elements)
"""
import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
#-------------------------------------------------------------------------------
# set display width
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)
#-------------------------------------------------------------------------------
'''Parameters'''
alpha=[-1,2]
beta=[-1,2]
#n=3
#M,SAM=[1000,100]
num_dup=20
#-------------------------------------------------------------------------------
'''read data '''
grid_num=41
alpha_expo_range=np.array([alpha[0],alpha[1]])                   # exponential distribution of alpha
grid_alpha=np.linspace(alpha_expo_range[0],alpha_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
alpha_list=10**grid_alpha

beta_expo_range=np.array([beta[0],beta[1]])
grid_beta=np.linspace(beta_expo_range[0],beta_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
beta_list=10**grid_beta

data_list=[[],[]]
for n in [3,5,10,15]:
	result=[[[] for i in range(grid_num)] for i in range(grid_num)]
	all_result=[[[[np.nan] for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
	
	for alpha_cluster in range(0,grid_num):                               # how many figures or cell divisions
		for beta_cluster in range(0,grid_num):                           # how many figures or cell divisions
	
			max_dup=[[] for i in range(num_dup)]
			for dup in range(0,num_dup):                               # how many figures or cell divisions
				max_dup[dup]=np.loadtxt('./data_n%s_m1000_sam100/%s_%s_%s_%s.txt'%(n,dup,n,alpha_cluster,beta_cluster))
				all_result[dup][alpha_cluster][beta_cluster]=max_dup[dup]
	
			"find the maximum over duplicates"
			grate_list=[]
			for dup in range(0,num_dup):
				grate_list.append(max_dup[dup][n,0])
	
			"save maximum"
			result[alpha_cluster][beta_cluster]=max_dup[np.argmax(grate_list)]
	
	"""classify the duplicates: store the 10 duplicates' data"""
	'''categaries: ND -1, RD=0; igd=1; isd=2, id=igsd=3  g->gg=1 no differentiation;'''
	all_result_classify=[[[0 for i in range(grid_num)] for i in range(grid_num)]for dup in range(num_dup)]   # set all as RD
	all_result_grate=[[[0 for i in range(grid_num)] for i in range(grid_num)]  for dup in range(num_dup)]
	
	igsd_s_classify=[[[0 for i in range(grid_num)] for i in range(grid_num)]for dup in range(num_dup)]   # set all as RD
	igsd_g_classify=[[[0 for i in range(grid_num)] for i in range(grid_num)]for dup in range(num_dup)]   # set all as RD
	
	for dup in range(0,num_dup):                               # how many figures or cell divisions
		for alpha_cluster in range(0,grid_num):                               # how many figures or cell divisions
			for beta_cluster in range(0,grid_num):                           # how many figures or cell divisions
	
				item =all_result[dup][alpha_cluster][beta_cluster]
				all_result_grate[dup][beta_cluster][alpha_cluster]=item[n][0]
	
				item1=item[:n,:]
	#			for i in range(1,n+1):
				for i in range(1,2):
					# ND
					if all(x==1 for x in item1[:,5]):                          # check whether last column f_g->gg are 1
						all_result_classify[dup][beta_cluster][alpha_cluster]=-1  
					# ISD	
					elif all(x==1 for x in item1[-i:,0]) and item1[-1:,5]!=1: 
						all_result_classify[dup][beta_cluster][alpha_cluster]=2 
					# IGD
					elif all(x==1 for x in item1[-i:,5]) and item1[-1:,0]!=1: 
						all_result_classify[dup][beta_cluster][alpha_cluster]=1 	
					# ID=IGSD
					elif all(x==1 for x in item1[-i:,0])and all(x==1 for x in item1[-i:,5]) : 
						all_result_classify[dup][beta_cluster][alpha_cluster]=3 
						for i in range(1,n+1):
							if all(x==1 for x in item1[-i:,0]):
								igsd_s_classify[dup][beta_cluster][alpha_cluster]=i
								
							if all(x==1 for x in item1[-i:,5]) : 
								igsd_g_classify[dup][beta_cluster][alpha_cluster]=i
	
	all_classify_narry=np.array([[np.array(i) for i in item] for item in all_result_classify])
	all_grate_narry=np.array([[np.array(i) for i in item] for item in all_result_grate])
	igsd_s_classify_arr=np.array([[np.array(i) for i in item] for item in igsd_s_classify])
	igsd_g_classify_arr=np.array([[np.array(i) for i in item] for item in igsd_g_classify])
	
	"percentage of 10"
	frac_igsd=[[[] for i in range(grid_num)] for i in range(grid_num)]
	fra_s=[[[] for i in range(grid_num)] for i in range(grid_num)]
	fra_g=[[[] for i in range(grid_num)] for i in range(grid_num)]
	
	for alpha_cluster in range(0,grid_num):                               # how many figures or cell divisions
		for beta_cluster in range(0,grid_num):                           # how many figures or cell divisions
	#		alpha_cluster,beta_cluster=30,30
			dup_num=np.count_nonzero(~np.isnan(all_classify_narry[:,alpha_cluster,beta_cluster]))
			# number of igsd
			igsd_num=np.count_nonzero(all_classify_narry[:,alpha_cluster,beta_cluster]==3) # num of igsd
			if igsd_num==0:
				frac_igsd[alpha_cluster][beta_cluster]=0
				fra_s[alpha_cluster][beta_cluster]=0
				fra_g[alpha_cluster][beta_cluster]=0
			else:
				frac_igsd[alpha_cluster][beta_cluster]=igsd_num/dup_num #fraction of igsd
				fra_s[alpha_cluster][beta_cluster]=n+1-np.sum(igsd_s_classify_arr[:,alpha_cluster,beta_cluster])/igsd_num
				fra_g[alpha_cluster][beta_cluster]=n+1-np.sum(igsd_g_classify_arr[:,alpha_cluster,beta_cluster])/igsd_num
	
	frac_igsd_arr=np.array([np.array(l) for l in frac_igsd])
	fra_s_arr=np.array([np.array(l1) for l1 in fra_s])
	fra_g_arr=np.array([np.array(l) for l in fra_g])
	#
	frac_igsd_arr[frac_igsd_arr == 0] = 'nan'
	fra_s_arr[fra_s_arr == 0] = 'nan'
	fra_g_arr[fra_g_arr == 0] = 'nan'
	data_list[0].append(fra_g_arr)
	data_list[1].append(fra_s_arr)
	
data_list=np.array([np.array(i) for i in data_list])
#np.save("data_list",data_list)

with open('./data_list.npy', 'rb') as f:
    data_list = np.load(f)

#================color list============================
import matplotlib as mpl
import matplotlib.pyplot as plt
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

#def fig(c1,c2,c3,num_c):
#	fig, ax = plt.subplots(figsize=(7.3, 1))
#	for x in range(2*num_c+2):
#	    ax.axvline(x, color=gra_c_list(c1,c2,c3,num_c)[x], linewidth=7) 
#	plt.show()
#	return 

#fig('#f5f5f5','#bf812d','#543005',31)
#-------------------------------------------------------------------------------
"draw figures"
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator

color1=gra_c_list('w','#bf812d','#543005',101)[10:]     # brown
color2=gra_c_list('w','#4292c6','#08306b',101)[10:]
cmap1 = LinearSegmentedColormap.from_list('mycmap', color1)
cmap2 = LinearSegmentedColormap.from_list('mycmap', color2)

max_fra=np.nanmax(1)
norm1 = plt.Normalize(0, max_fra)
cmap1.set_under('w')
cmap2.set_under('w')
cmap=[cmap1,cmap2]

im=[i for i in range(4)]
#pos_list1=[0.9, 0.47]
pos_list2=[0.055, 0.29,0.52,0.755]
#nbins_list=[2,2,2,3]
n_list=[3,5,10,15]

fig, ax = plt.subplots(1, 4, gridspec_kw={'width_ratios': [1,1,1,1],'wspace':0, 'hspace':0},figsize=(5, 3))   # 10.  5.5

for i in range(4):
#	for j in range(4):
		n=n_list[i]
		im[i]=ax[i].imshow(data_list[0][i], interpolation=None, origin='lower',cmap=cmap[0], norm=plt.Normalize(1, n))
		ax[i].tick_params(direction='in', length=1.3, width=1, colors='k')
		ax[i].tick_params(axis='both', which='major', labelsize=8)
		ax[i].spines['bottom'].set_color('gray')
		ax[i].spines['top'].set_color('gray') 
		ax[i].spines['right'].set_color('gray')
		ax[i].spines['left'].set_color('gray')
		ax[i].tick_params(axis='both', colors='k')		
		
		ax[i].xaxis.set_ticklabels([])
		ax[i].set_yticklabels([])

		cbar_ax1 = fig.add_axes([pos_list2[i],0.32, 0.2, 0.02])
		cbar1=fig.colorbar(im[i], cax=cbar_ax1, orientation="horizontal",norm=norm1, boundaries=None)
		cbar1.ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		cbar1.ax.tick_params(labelsize=8,length=0,direction='in')
		cbar1.ax.xaxis.set_ticks_position('top')
		cbar1.ax.xaxis.set_label_position('top')
		cbar1.ax.xaxis.set_tick_params(pad=-1)
		cbar1.outline.set_visible(False)

plt.show()
fig.savefig('./EID_timing_3-5-10-15_R.pdf', bbox_inches='tight')   # save figures

#==================
fig, ax = plt.subplots(1, 4, gridspec_kw={'width_ratios': [1,1,1,1],'wspace':0, 'hspace':0},figsize=(5, 3))   # 10.  5.5

for i in range(4):
#	for j in range(4):
		n=n_list[i]
		im[i]=ax[i].imshow(data_list[1][i], interpolation=None, origin='lower',cmap=cmap[1], norm=plt.Normalize(1, n))
		ax[i].tick_params(direction='in', length=1.3, width=1, colors='k')
		ax[i].tick_params(axis='both', which='major', labelsize=8)
		ax[i].spines['bottom'].set_color('gray')
		ax[i].spines['top'].set_color('gray') 
		ax[i].spines['right'].set_color('gray')
		ax[i].spines['left'].set_color('gray')
		ax[i].tick_params(axis='both', colors='k')		
		
		ax[i].xaxis.set_ticklabels([])
		ax[i].set_yticklabels([])

		cbar_ax1 = fig.add_axes([pos_list2[i],0.32, 0.2, 0.02])
		cbar1=fig.colorbar(im[i], cax=cbar_ax1, orientation="horizontal",norm=norm1, boundaries=None)
		cbar1.ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		cbar1.ax.tick_params(labelsize=8,length=0,direction='in')
		cbar1.ax.xaxis.set_ticks_position('top')
		cbar1.ax.xaxis.set_label_position('top')
		cbar1.ax.xaxis.set_tick_params(pad=-1)
		cbar1.outline.set_visible(False)

plt.show()
#fig.savefig('./EID_timing_3-5-10-15_S.pdf', bbox_inches='tight')   # save figures
