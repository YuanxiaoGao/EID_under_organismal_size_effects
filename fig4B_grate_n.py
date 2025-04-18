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
size_list=np.array([2,n_cluster],dtype=int)      # alpha beta =-1.5 or 1.5
bc_points=[1,10]
#-------------------------------------------------------------------------------
'''read data '''
grid_num=2
result=[[[[[]for n in range(2,n_cluster)]for c in range(2)] for b in range(2)] for i in range(4)]   # save the four categories

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

#-------------DRAW FIGURE A B C D----------------------------------------------------------------
import matplotlib.pyplot as plt
result[3][0][1][0]  #IGSD
result[2][1][1][0]  #ISD
result[1][1][1][20]   #IGD
result[0][1][1][3]   #RD

mark_list=['8','^','8','x']
label_list=[r'$RD_{RS}$',r'$RD_S$',r'$RD_R$',r'$EID$']  # RD_s--soma reversible
color_list=['grey','#8c510a','#bf812d','#dfc27d','#2166ac']
size=[i for i in range(2,n_cluster)]

#------------------one figure----------------
fig, ax = plt.subplots(1, 1,figsize=(5.5, 5))
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.grid(False)
plt.xlabel(r'Organismal size, $2^n$',fontsize=16)
plt.ylabel(r'Reproduction rate, $\lambda(n)$',fontsize=16)
#-----colormap--------------
alpha0=0.75
markers=8

for b in range(0,2):  
	for c in range(0,2):    #grid_num
		# for others		
		for stra_cluster in [0,1,2,3]: 
			grate_mean=[]
			grate_std=[]
			for n in range(2,n_cluster):
				grate1=[]
				for dup in range(num_dup):	
					grate1.append(result[stra_cluster][b][c][n-2][dup][-1][0])
				grate_mean.append(np.mean(np.array(grate1)))
				grate_std.append(np.std(np.array(grate1)))
			# errorbar
			ax.errorbar(size,grate_mean,grate_std,linestyle='-',linewidth=0.6,marker="h",
				markersize=markers,label=label_list[stra_cluster],color=color_list[stra_cluster+1],alpha=alpha0)
		# frame color
		color_f="darkgray"
		ax.spines['bottom'].set_color(color_f)
		ax.spines['top'].set_color(color_f) 
		ax.spines['right'].set_color(color_f)
		ax.spines['left'].set_color(color_f)
		ax.tick_params(direction='in', length=2, width=1, colors='k')
		# legend for b c
		ax.text(size[18-3*b],grate_mean[18+c]*1.2-b*0.45,r"$b=%s$, $c=%s$"%(bc_points[int(b)],bc_points[int(c)]),rotation=b*25,rotation_mode='anchor',fontsize=11,alpha=1)
# ND
ax.errorbar(size,[np.log(2) for i in range(len(size))],[0 for i in range(len(size))],
	linestyle='-',linewidth=1,marker="h",markersize=markers,color=color_list[0],alpha=alpha0*0.8 )							
#------x and y labels and limits--------
# artificial legend
import matplotlib.lines as mlines
legend_list=[]
legend_list1=mlines.Line2D([], [], color=color_list[0], linestyle='None',marker="h",markersize=9,label=r'$ND$',alpha=alpha0)
legend_list.append(legend_list1)
for i in range(4):
	target=mlines.Line2D([], [], color=color_list[i+1], linestyle='None',marker="h",markersize=9,label=label_list[i],alpha=alpha0) 
	legend_list.append(target)
	
pos = ax.get_position()
ax.set_position([pos.x0, pos.y0, pos.width, pos.height])
ax.legend(handles=legend_list,loc='upper left', bbox_to_anchor=(-0.01, 0.99),ncol=1,fontsize=11,frameon=False)

plt.tight_layout()
plt.show()
fig.savefig('./fig4_lambda_n_infinite_prediction.pdf', bbox_inches = 'tight')   # save figures
