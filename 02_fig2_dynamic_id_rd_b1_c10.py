#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29 13:29:38 2022
@author: gao
"""

"""
f: fraction of germ cells
p1: p_{g->s}
p2: q_{s->g}
DP=np.array([[s_ss,s_gs,s_gg];[g_ss,g_gs,g_gg]]);
RD,IGD,ISD,IGSD : stra_cluster 0,1,2,3
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
#-------------------------------------------------------------------------------
# set display width
np.get_printoptions()['linewidth']
np.set_printoptions(linewidth=2000)
#-------------------------------------------------------------------------------
'''Parameters'''
alpha=1
beta=1
n_list=[3,5,10,15]
num_dup=100
bc_points=[1,10]
#-------------------------------------------------------------------------------
'''read data'''
grid_num=2
result=[[[]for stra in range(4)] for n in range(len(n_list))]   # save the four categories
for dup in range(num_dup):  
	for nth in range(len(n_list)): 
		n= n_list[nth]
#		for bth in range(0,grid_num):   
		for bth in range(1):   

			b=bc_points[bth]
#			for cth in range(0,grid_num): 
			for cth in range(1,2): 
#			for cth in range(1): 
				c=bc_points[cth]
				for stra_cluster in range(4): 
					if os.path.exists('./data_bc110/%s_%s_%s_%s_%s.txt'%(dup,n,bth,cth,stra_cluster)):
						data=np.loadtxt('./data_bc110/%s_%s_%s_%s_%s.txt'%(dup,n,bth,cth,stra_cluster))						
						result[nth][stra_cluster].append(data)
# result=[ [n=3; [rd],[igd],[isd],[igsd] ];  [n=5; []   ....]       ]
result1=[[[]for stra in range(2)] for n in range(len(n_list))]    # only rd and id
for nth in range(4):
	result1[nth][0]=result[nth][0]+result[nth][1]+result[nth][2]
	result1[nth][1]=result[nth][3]
	
result2=[[[]for stra in range(2)] for n in range(len(n_list))]    # only rd and id
for nth in range(4): 
	result2[nth][1]=result1[nth][1]
	grate_list=[]
	for i in range(3*num_dup):
		data=result1[nth][0][i]
		grate_list.append(data[-1][0])
	ind=sorted(range(len(grate_list)), key=lambda i: grate_list[i])[-num_dup:]

	for j in ind:
		result2[nth][0].append(result1[nth][0][j])

#---------------------------------------------------------------------------
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
	t_list=[]
	b_list=[]
	c_list=[]
	t=0
	for i in range(1,n+1):
		cost=1+para[1][0]*(fgs[i-1][0]*p12[i][0]+para[1][1]*fgs[i-1][1]*p12[i][1])
		bene=1+para[0][0]*(fgs[i-1][1])**para[0][1]
		b_list.append(bene)
		c_list.append(cost)
		t_i=cost/bene
		t_list.append(t_i)
		t=t+t_i
	growth_rate=np.log(2**n*fgs[n][0])/t
	r_list=1/np.array(t_list)
	return [b_list,c_list,r_list,np.array(fgs)[:,0],growth_rate]

para=np.array([[b,alpha],[c,beta]])                    # para=n	np.array([[b,alpha,...],[c,beta,...]])
#---------------------------------------------------------------------------
"draw figures"
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator
# ----set fig frame--------
fig, ax = plt.subplots(1, 5, gridspec_kw={'width_ratios': [1,1,1,1,0.2]},figsize=(13, 3)) #,constrained_layout = True)
#fig.suptitle(r'b=%s, c=%s'%(b,c), fontsize=14,y=1.12)
fig.tight_layout(pad=1.2)
alp_small=0.9
#-----colormap--------------
color_list=[['#dfc27d','#bf812d','#8c510a','#673b07'],['#c6dbef','#9ecae1','#4292c6','#08306b']]   # blue brown
linestyle_list=["-","-"]
label_list=["RD","EID"]
ylabel_list=[r"Benefit composition, $B(i)$",r"Cost composition, $C(i)$",r"Cell division rate, $r(i)$",r"Fraction of R cells, $f_{R}(i)$",r"Reproduction rate, $\lambda(i)$"]
for i in range(4):   # n 
	n=n_list[i]
	x = np.linspace(1,n_list[i],n_list[i])
	for sta in range(2):   # rd    id
		bcrg_list=[[],[],[],[],[]]    # b  c  r  grate      all data
		for dup in range(num_dup):
			dp=result2[i][sta][dup][0:-1,:]
			bcrg=Growth_num(n,dp,para)			
			for k in range(5):
				bcrg_list[k].append(np.array(bcrg[k]))
		mean_list=[]
		std_list=[]
		for l in range(5):
			data=np.array(bcrg_list[l])
			mean_list.append(np.mean(data, axis=0))
			std_list.append(np.std(data, axis=0))
			
		x1=np.linspace(0,n_list[i],n_list[i]+1)
		if n==15:
			for j in range(4):   # r, b ,c ,grate
		    # ND			
				if j==3:
					ax[j].scatter(x1, [1 for i in range(16)],s=45,marker="o",color="#969696",alpha=1)
					ax[j].plot(x1,[1 for i in range(16)], linestyle = '-',linewidth=2,color="#969696",alpha=1,zorder=-3)		
				else:
					ax[j].scatter(x1[1:], [1 for i in range(15)],s=45,marker="o",color="#969696",alpha=1)
					ax[j].plot(x1[1:],[1 for i in range(15)], linestyle = '-',linewidth=2,color="#969696",alpha=1,zorder=-3)		
					
		for j in range(5):   # r, b ,c ,grate
			if j==4:            # growth rate figure
				#ND
				ax[j].scatter(-1,np.log(2),s=80,marker="H",color="#969696",alpha=1)
				# EID  RD
				ax[j].scatter([sta],mean_list[-1],s=80,marker="H",facecolors=color_list[sta][i], edgecolors=color_list[sta][i],alpha=alp_small)
				ax[j].errorbar([sta],mean_list[-1],std_list[-1],color=color_list[sta][i])
				ax[j].set_xlim(-1.5,1.5)
				ax[j].set_xlabel('ND  RD  EID',fontsize=12, labelpad=9)
			elif j==3:         # fraction of germ-like cells
				ax[j].plot(x1,mean_list[j], linestyle = linestyle_list[sta],linewidth=2,color=color_list[sta][i],alpha=alp_small,zorder=3)	
				ax[j].fill_between(x1,mean_list[j]-std_list[j],mean_list[j]+std_list[j],color=color_list[sta][i],alpha=alp_small*0.1,zorder=2)
				ax[j].scatter(x1,mean_list[j],s=50,marker="o",facecolors=color_list[sta][i], edgecolors=color_list[sta][i],linestyle='-',linewidths=1,alpha=alp_small,zorder=4)
				ax[j].set_xlabel(r'$i$th cell division',fontsize=14)
			else:	
				ax[j].plot(x,mean_list[j], linestyle = linestyle_list[sta],linewidth=2,color=color_list[sta][i],alpha=alp_small,zorder=3)	
				ax[j].fill_between(x,mean_list[j]-std_list[j],mean_list[j]+std_list[j],color=color_list[sta][i],alpha=alp_small*0.1,zorder=2)
				ax[j].scatter(x,mean_list[j],s=50,marker="o",facecolors=color_list[sta][i], edgecolors=color_list[sta][i],linestyle='-',linewidths=1,alpha=alp_small,zorder=4)
				ax[j].set_xlabel(r'$i$th cell division',fontsize=14)
			# x and y label
			ax[j].set_ylabel(ylabel_list[j],fontsize=14,labelpad=0)
			ax[j].tick_params(axis='y', which='major', pad=1.5)				
			# xy ticks
			ax[j].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
			ax[j].tick_params(direction='in', length=2, width=1, colors='k',grid_color='r', grid_alpha=1)     # inout
			ax[j].locator_params(axis='x', nbins=n_list[i])
			ax[j].locator_params(axis='y', nbins=5)
			ax[j].xaxis.set_major_locator(MaxNLocator(integer=True))
		ax[4].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
		ax[0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
		ax[3].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
		ax[4].set_xticks([])
   # legend
	colo=['#dfc27d','#bf812d','#8c510a','#673b07','#c6dbef','#9ecae1','#4292c6','#08306b']
	sc_list=[]
	label_list=[r"$RD, n=3$", r"$RD, n=5$",r"$RD, n=10$",r"$RD, n=15$",r"$EID, n=3$", r"$EID, n=5$",r"$EID, n=10$",r"$EID, n=15$"]
	for i in range(8):
		sc1=plt.scatter([],[],s=100,facecolors=colo[i],edgecolors=colo[i],linestyle='-')
		sc_list.append(sc1)
	ax[2].legend(sc_list, label_list,numpoints=1,ncol=8,frameon=False,loc='upper', bbox_to_anchor=(2.7, 1.18), shadow=False,fontsize = 12)
	if b==1 and c==10:
		ax[1].text(8,1.6,r"$b=%s,c=%s$"%(b,c),fontsize=12, style='italic')
		ax[0].text(13,1.005,r"$ND$",fontsize=12, color='#969696',style='italic')
	 
plt.show()
fig.savefig('./fig2B_dynamic_id_rd-b%s_c%s.pdf'%(b,c), bbox_inches = 'tight')   # save figures
