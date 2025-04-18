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
num_dup=20
alpha=1
beta=1
n_cluster=30 #maimum n
n_size=29   # nth
size_list=np.array([2,n_size],dtype=int)      # alpha beta =-1.5 or 1.5
grid_num=2
bc_points=[1,10]

#-------------------------------------------------------------------------------
'''read data '''
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

"get the frequency and time of the random five strategies"
  
def Fration(n,dp):                              # dp is a sp series with number n
	p12=[np.array([0,0])]                       #p_g->s & p_s->g     total n+1 items
	fgs=[np.array([1,0])]                       #total n+1 items
	for i in range(n):
		p1=dp[i][3]+0.5*dp[i][4]
		p2=dp[i][2]+0.5*dp[i][1]
		p12.append(np.array([p1,p2]))
	for i in range(n+1):
		matrix=np.array([[1-p12[i][0],p12[i][1]],[p12[i][0],1-p12[i][1]]])
		fgs.append(matrix.dot(fgs[i]))
	del fgs[0]	
	
	return [p12,fgs]  

def Growth_num(n,dp,para):
	p12,fgs=Fration(n,dp)
	t_list=[0]
	t=0
	for i in range(1,n+1):
		cost=1+para[1][0]*(fgs[i-1][0]*p12[i][0]+para[1][1]*fgs[i-1][1]*p12[i][1])
		bene=1+para[0][0]*(fgs[i-1][1])**para[0][1]
		t_i=cost/bene
		t_list.append(t_i)
		t=t+t_i
#		print(t_i)
	growth_rate=np.log(2**n*fgs[n][0])/t
	return [growth_rate,t_list]	

#-------------DRAW FIGURE A B C D----------------------------------------------------------------
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.ticker import FormatStrFormatter

#label_list=['ND','RD','IGD','ISD','IGSD']
label_list=['ND','RD','IGD','ISD','IGSD']
color_list=['#969696','#8c510a','#bf812d','#dfc27d','#2166ac']
size=[i for i in range(0,n_size+1)]
		
for b_in in range(1,2): 
	for c_in in range(1): 
		b=bc_points[b_in]
		c=bc_points[c_in]		
		fig, (ax1, ax2) = plt.subplots(1,2)
		plt.rcParams["figure.figsize"] = [11, 5]
		plt.rcParams["figure.autolayout"] = True
		plt.rcParams['axes.linewidth'] =1		
		ax1.patch.set_edgecolor('r')  
		ax1.patch.set_linewidth('0.01')  		
		# frame color
		color_f="lightgrey"
		ax1.spines['bottom'].set_color(color_f)
		ax1.spines['top'].set_color(color_f) 
		ax1.spines['right'].set_color(color_f)
		ax1.spines['left'].set_color(color_f)
		
		ax2.spines['bottom'].set_color(color_f)
		ax2.spines['top'].set_color(color_f) 
		ax2.spines['right'].set_color(color_f)
		ax2.spines['left'].set_color(color_f)
		ax1.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
		ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
		
		ax1.tick_params(direction='in', length=2, width=1, colors='k',
	               grid_color='r', grid_alpha=1)
		ax2.tick_params(direction='in', length=2, width=1, colors='k',
	               grid_color='r', grid_alpha=1)	
		ax1.locator_params(axis='y', nbins=5)
		ax2.locator_params(axis='y', nbins=4)
		
		ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
		ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
		labelpad=0.2
		ax1.set_xlabel(r'$i$th Cell division',fontsize=16, labelpad=labelpad+2)
		ax2.set_xlabel(r'$i$th Cell division',fontsize=16, labelpad=labelpad+2)
		ax1.set_ylabel(r'Frequency of cells, $f_{R}(i)$, $f_{S}(i)$',fontsize=16, labelpad=labelpad)
		ax2.set_ylabel(r'Cell division rate, $r(i)$',fontsize=16, labelpad=labelpad)
		
		marker_s=3.5
	    # ND
		for n in range(0,n_size):
			# germ
			ax1.scatter(n, 1,s=50,marker="o",color=color_list[0],alpha=1)
			ax1.plot(size,[1 for i in range(n_size+1)] , linestyle = '-',linewidth=2,color=color_list[0],alpha=1,zorder=-2)		
			# soma
			ax1.scatter(n, 0,s=50,marker="o",facecolors="w",edgecolors=color_list[0],alpha=1,zorder=-1)
			ax1.plot(size,[0 for i in range(n_size+1)] , linestyle = '-',linewidth=2,color=color_list[0],alpha=1,zorder=-2)		
			
			ax2.scatter(n, 1,s=50,marker="o",color=color_list[0],alpha=1,zorder=-1)
			ax2.plot(size,[1 for i in range(n_size+1)] , linestyle = '-',linewidth=2,color=color_list[0],alpha=1,zorder=-2)	
		
		alp_small=0.2
		# for others		
		for stra_cluster in range(4): 	
			all_stra=result[stra_cluster][b_in][c_in][n_size-2]
	
			fg_list=[]
			t_list=[]
			for run in range(num_dup):
				one=all_stra[run]
				dp=one[0:n_size,:]
				fg=np.array(Fration(n_size,dp)[1])[:,0]
				t=np.array(Growth_num(n_size,dp,np.array([[b,1],[c,1]]))[1])
				fg_list.append(fg)
				t_list.append(t)
				# lines
				# germ
				ax1.plot(size,fg, linestyle = '-',linewidth=0.1,color=color_list[stra_cluster+1],alpha=alp_small)	
				ax1.scatter(size,fg,s=1,marker="o",label=label_list[stra_cluster+1],color=color_list[stra_cluster+1],alpha=alp_small)
				#soma
				ax1.plot(size,1-fg, linestyle = '-',linewidth=0.1,color=color_list[stra_cluster+1],alpha=alp_small)	
				ax1.scatter(size,1-fg,s=1,marker="o",label=label_list[stra_cluster+1],facecolors="w",edgecolors=color_list[stra_cluster+1],alpha=alp_small)

				ax2.plot(size,1/t, linestyle = '-',linewidth=0.1,color=color_list[stra_cluster+1],alpha=alp_small,zorder=1)	
				ax2.scatter(size,1/t,s=1,marker="o",label=label_list[stra_cluster+1],color=color_list[stra_cluster+1],alpha=alp_small,zorder=3)

			fg_np=np.array(fg_list)
			fg_mean=np.mean(fg_np,axis=0)
			fg_std=np.std(fg_np,axis=0)
			fs_std=np.std(1-fg_np,axis=0)
			
			t_np=np.array(t_list)
			t_mean=np.mean(1/t_np,axis=0)
			t_std=np.std(1/t_np,axis=0)
			
			# shades
			alpha_l=0.9
#			#germ
			ax1.scatter(size,fg_mean ,s=50,marker="o",color=color_list[stra_cluster+1],alpha=1,zorder=2)
			ax1.plot(size,fg_mean, linestyle = '-',linewidth=2, label=label_list[stra_cluster+1], color=color_list[stra_cluster+1],alpha=alpha_l,zorder=1)		
			ax1.fill_between(size,fg_mean-fg_std,fg_mean+fg_std,color=color_list[stra_cluster+1],alpha=0.1,zorder=3)

			# soma
			ax1.scatter(size,1-fg_mean ,s=50,marker="o",facecolors="w",edgecolors=color_list[stra_cluster+1],alpha=1,zorder=4)
			ax1.plot(size,1-fg_mean, linestyle = '-',linewidth=2,label=label_list[stra_cluster+1],color=color_list[stra_cluster+1],alpha=alpha_l,zorder=1)		
			ax1.fill_between(size,1-fg_mean-fs_std,1-fg_mean+fs_std,color=color_list[stra_cluster+1],alpha=0.1,zorder=5)

			ax2.scatter(size,t_mean ,s=50,marker="o",color=color_list[stra_cluster+1],alpha=1)
			ax2.plot(size,t_mean, linestyle = '-',linewidth=2,color=color_list[stra_cluster+1],alpha=alpha_l)		
			ax2.fill_between(size,t_mean-t_std,t_mean+t_std,color=color_list[stra_cluster+1],alpha=0.1)

			if stra_cluster==3:
				ax1.text(25,0.45,r"$f_R(i)$",fontsize=14)
				ax1.text(25,0.55,r"$f_S(i)$",fontsize=14)
				
		plt.tight_layout()
		plt.show()
		fig.savefig('./fg4c_r_n%s_b%s_c%s_%s.pdf'%(n_size,b,c,num_dup), bbox_inches = 'tight')   # save figures


