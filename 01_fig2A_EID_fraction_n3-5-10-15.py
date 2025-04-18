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
n=4
M,SAM=[1000,100]
num_dup=5
#-------------------------------------------------------------------------------
'''read data '''
grid_num=41
alpha_expo_range=np.array([alpha[0],alpha[1]])                   # exponential distribution of alpha
grid_alpha=np.linspace(alpha_expo_range[0],alpha_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
alpha_list=10**grid_alpha

beta_expo_range=np.array([beta[0],beta[1]])
grid_beta=np.linspace(beta_expo_range[0],beta_expo_range[1],num=grid_num,endpoint=True) # split log(a) into grid_num points
beta_list=10**grid_beta

frac_list=[0.4, 0.45, 0.45, 0.55, 0.6]
Matrix_igsd=[]
Matrix_g=[]
Matrix_s=[]
n_list=[3,5,10,15]
for n in n_list:

	result=[[[] for i in range(grid_num)] for i in range(grid_num)]
	all_result=[[[[np.nan] for i in range(grid_num)] for i in range(grid_num)] for dup in range(num_dup)]
	
	for alpha_cluster in range(0,grid_num):                               # how many figures or cell divisions
		for beta_cluster in range(0,grid_num):                           # how many figures or cell divisions
	
			max_dup=[[] for i in range(num_dup)]
			for dup in range(0,num_dup):                               # how many figures or cell divisions
				max_dup[dup]=np.loadtxt('/Users/gao/Desktop/dds_size/Code/06_analytic/00_dynamic_bc/data_n%s_m1000_sam100/%s_%s_%s_%s.txt'%(n,dup,n,alpha_cluster,beta_cluster))
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
	
	# convert to fraction of cell development
	fra_s_arr=fra_s_arr/n
	fra_g_arr=fra_g_arr/n
	#
	frac_igsd_arr[frac_igsd_arr == 0] = 'nan'
	fra_s_arr[fra_s_arr == 0] = 'nan'
	fra_g_arr[fra_g_arr == 0] = 'nan'
	frac_list.append(np.nanmax(frac_igsd_arr))
	
	Matrix_igsd.append(frac_igsd_arr)
	Matrix_g.append(fra_g_arr)
	Matrix_s.append(fra_s_arr)
	
np.save("Matrix_igsd",Matrix_igsd)
np.save("Matrix_g",Matrix_g)
np.save("Matrix_s",Matrix_s)



#======================draw figures==========================================
with open('Matrix_igsd.npy', 'rb') as f:
    Matrix_igsd = np.load(f)
with open('Matrix_g.npy', 'rb') as f:
    Matrix_g = np.load(f)
with open('Matrix_s.npy', 'rb') as f:
    Matrix_s = np.load(f)
	 
#================color list============================
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker

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
M_igsd=[]
for i in [0,1,2,3]:
	M_igsd.append(Matrix_igsd[i])


#-------------------------------------------------------------------------------
"draw figures"
from matplotlib.colors import LinearSegmentedColormap
# ----set fig frame--------
#	fig, ax = plt.subplots(1, 1, figsize=(3.5, 3.5))
fig, ax = plt.subplots(1, 4, gridspec_kw={'width_ratios': [1,1,1,1]},figsize=(13, 3.3))
plt.subplots_adjust(wspace=0, hspace=0)
#-----colormap--------------
color1=gra_c_list('w','#4292c6','#08306b',31)[3:]
#color0=color1=color2=gra_c_list('w','#9ecae1','#08306b',31)[5:-3]
cmap1 = LinearSegmentedColormap.from_list('mycmap', color1)
#	max_fra=np.nanmax(frac_igsd_arr)
norm1 = plt.Normalize(0, 0.6)	
cmap1.set_under('w')

n0=[3,5,10,15]
im=[i for i in range(4)]
for i in range(4):
	#----figure------------------
	im[i] = ax[i].imshow(M_igsd[i], interpolation=None, origin='lower',cmap=cmap1, norm=norm1)
	# title
#	ax[i].set_title(r'Organismal size, $n=$%s'%n0[i],fontsize=13)
	ax[i].set_title(r'Organismal size $2^n, n=%s$'%n0[i],fontsize=13)

	# x and y label
	ax[i].set_xlabel(r'Differentiation benefit, $b$',fontsize=14)
	if i==0:
		ax[i].set_ylabel(r'Differentiation cost, $c$',fontsize=14,labelpad=0)
	
	# xy ticks
	ax[i].tick_params(direction='in', length=3, width=1, colors='k',pad=1.5 )
	# artifical x and y ticks
	
	if i==0:
		test0=np.linspace(0,grid_num-1,4,endpoint=True)
		ax[i].set_xticks( test0, minor=False)
		test=np.linspace(alpha_expo_range[0],alpha_expo_range[1],4,endpoint=True)
		x_label=[r'$10^{'+str(int(i))+'}$' for i in test]
		ax[i].xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

		ax[i].set_yticks( test0, minor=False)
		y_test=np.linspace(beta_expo_range[0],beta_expo_range[1],4,endpoint=True)
		y_label=[r'$10^{'+str(int(i))+'}$' for i in y_test]
		ax[i].yaxis.set_major_formatter(mpl.ticker.FixedFormatter(y_label))
		ax[i].tick_params(axis='both', which='in',pad=1.5 ,labelsize=10)
	else:
		test0=np.linspace(0,grid_num-1,4,endpoint=True)
		ax[i].set_xticks( test0[1:], minor=False)
		test=np.linspace(alpha_expo_range[0],alpha_expo_range[1],4,endpoint=True)
		x_label=[r'$10^{'+str(int(i))+'}$' for i in test[1:]]
		ax[i].xaxis.set_major_formatter(mpl.ticker.FixedFormatter(x_label))

		ax[i].set_yticks([] )

	pos_list=0.905

	cbar_ax1 = fig.add_axes([pos_list, 0.15, 0.007, 0.71])
	cbar1=fig.colorbar(im[i], cax=cbar_ax1, orientation="vertical",norm=norm1, boundaries=None)
	cbar1.outline.set_visible(False)
	tick_locator = ticker.MaxNLocator(nbins=7)
	cbar1.locator = tick_locator
	cbar1.update_ticks()
	cbar1.ax.set_ylabel(r'Fraction of $EID$ being optimal strategy', rotation=90,fontsize=11)
	cbar1.ax.tick_params(labelsize=8,length=2,direction='in')

plt.show()
	
#fig.savefig('./fig2A_EID_n3-15_dup%s_n%s_V2.pdf'%(num_dup,n), bbox_inches='tight')   # save figures
	
