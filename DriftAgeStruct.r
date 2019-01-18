# **********************************************************************
# R Script for DriftAgeStruct
# By ADJ Overall
# University of Brighton
# Last Update January 2019
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# **********************************************************************
#
#
# To run this script, copy into the working directory containing your input files and type 
# > source(“DriftAgeStruct.r”).
# 
# leslie.drift = function(infile,B,N,T,iteration,initial,s)
  #  Arguments:
  #    infile is the input file
  #    B is the starting frequency of the mutation
  #    N is the population size
  #    T is the number of generations the simulation runs for
  #    iteration is the number of simulation repeats
  #    initial is the number of iterations of the leslie matrix run to stabilise
  #    s is the selection coefficient
  #  Example of use:
  #  > leslie.drift(infile=“datafile.txt”    
  #  ,B=0.5,N=1000,T=1000,iteration=1000,initial=100,s=0.005)
# 

leslie.drift = function(infile,B,N,T,iteration,initial,recessive=T,s){  set.seed(12345)	
# read in lx and mx data
  ds = read.table(infile, header = TRUE)# define variables    HET = matrix(0,iteration,T) # matrix of heterozygosities  mean.HET=matrix(0,T)	# matrix of mean heterozygosities  dummyX=c(1:T)  # dummy for plot
  het.Dist=matrix(0,iteration,T) # heterozygosity at end time point
 
# B is specified starting allele frequency 

  A = 1 - B  NewB = B # NewB tracks the change in allele frequency   Nx = nrow(ds) # number of rows in dataset# Input data is lx and mx. lx is converted to Px
# From lx to Px  Px = matrix(0,Nx) # matrix of survival probabilities

# constructs Px matrix   Px[1] = 1  for(i in 2:Nx){    Px[i]=ds$lx[i]/prod(Px[1:(i-1)])  }   for(i in 1:(Nx-1)){    Px[i] = Px[i+1]  }# fertilities derived from mx data
  fx=matrix(0,Nx)  #Post-birth pulse  for (i in 2:(Nx)){    fx[i] = Px[i-1] * ds$mx[i]  }# construct leslie matrix
  leslie = matrix(0, Nx, Nx)  for(i in 1: Nx){ 
    leslie[1,i]= fx[i]
  }  for(i in 2: Nx){ 
    leslie[i,(i-1)]= Px[i-1]
  }# leslie matrix run for ds$initial iterations
# Let population age-structure stabilise to set nt.0  nt.0 = matrix(0,Nx)  nt.0[1] = N  for(i in 1:initial){      nt.1 = leslie %*% nt.0    nt.0 = nt.1  }  total = sum(nt.0)  for(i in 1:Nx){
    nt.0[i] = nt.0[i]/total
  }  nt.0 = nt.0 * N  NewBorn.Size = nt.0[1]# outputs the numbers for age-structured population.
  
  print(nt.0)

# A barplot can be produced:
  barplot(nt.0,beside=T,col="light blue",names=as.character(1:Nx),main="Age-Structure of Population",xlab="Age-class(x)",ylab="Frequency")

# If so, then put break in to pause simulation

  readline(prompt="Press [enter] to continue")

# Prepare figure for simulations

  plot(dummyX,type="n",xlab="Time",ylab="Mean Heterozygosity",cex.lab=1.5, cex.axis=1.5,xlim=c(0,T),ylim=c(0,0.5))
  
# Now introduce genotypes and drift to leslie matrix  
  for(j in 1:iteration){  # ds$iteration specifies number of generations    B = NewB  
    A = 1 - B
    Ny = 3 * Nx	#treble matrix to allow for each age class to have three 		#genotypes	    nt.y = matrix(0,Ny)	# number of genotypes per age-class    Px.y = matrix(0,Ny) # survival probability for each genotype    fx.y = matrix(0,Ny) # fertility distribution for each genotype    p = matrix(0,Ny)    allele_freq = matrix(0,T)    #rewrite leslie matrix     start=1    for(i in 1:Nx){      nt.y[start]=nt.0[i]; nt.y[start+1]=nt.0[i]; nt.y[start+2]=nt.0[i]; start = start +3    }    #mutliply each age-class by genotype frequencies    for(i in seq(1,by=3,Ny)){
      nt.y[i]=A^2*nt.y[i]
    }    for(i in seq(2,by=3,Ny)){
      nt.y[i]=2*A*B*nt.y[i]
    }    for(i in seq(3,by=3,Ny)){
      nt.y[i]=B^2*nt.y[i]
    }    #n111,n121,n221,n112,n122,n222 etc      #Assign each survival (Px) to a genotype    start=1    for(i in 1:(Nx -1)){      Px.y[start]=Px[i]; Px.y[start+1]=Px[i]; Px.y[start+2]=Px[i]; start = start +3    }    #S = matrix(0,Nx)
    k = 1
    for(i in 1:Nx){

      if(ds$mx[i]>0){
        Px.y[k+2] = (1 - s) * Px.y[k+2]
      }
      k = k+3
    }

    #Assign each fecundity (fx) to a genotype    start=1    for(i in 1:Nx){      fx.y[start]=fx[i]; fx.y[start+1]=fx[i]; fx.y[start+2]=fx[i]; start = start +3    }    mating = matrix(0,Ny,Ny)    for (i in 4:Ny){
      mating[i,i]=1
    }    leslie.y = matrix(0,Ny,Ny)    for(i in 1:Ny){
      leslie.y[1,i]=fx.y[i]
    }    for(i in 4:Ny){
      leslie.y[i,(i-3)]=Px.y[(i-3)]
    }    b = matrix(0,Ny)    hom=matrix(0,Ny)    het=matrix(0,Ny)    x=c(0,0,0)    for(ii in 1:T){
      for(iii in seq(4,by=3,(Nx-1))){        TOTAL = nt.y[iii]+ nt.y[iii+1]+ nt.y[iii+2]        total =Px.y[iii]*nt.y[iii]+Px.y[iii+1]*nt.y[iii+1]+ Px.y[iii+2]*nt.y[iii+2]        total = round(TOTAL- total)        if(total > 0){
          x=rmultinom(1,total,c(round(nt.y[iii])/total,round(nt.y[iii+1])/total,round(nt.y[iii+2])/total))
        }           nt.y[iii] = Px.y[iii]*nt.y[iii] + x[1]          nt.y[iii+1] = Px.y[iii+1]*nt.y[iii+1] + x[2]          nt.y[iii+2] = Px.y[iii+2]*nt.y[iii+2] + x[3]      }
      b = sum(fx.y * nt.y)      for(i in 1:Ny){
        p[i]=fx.y[i]*nt.y[i]
      }      for(i in seq(1,by=3,Ny)){
        hom[i]=p[i]
      }      for(i in seq(2,by=3,Ny)){
        het[i]=p[i]
      }      A = (sum(hom)+0.5*sum(het))/b      B = 1- A      if (A > 1){
        A = 1
      }      mating[1,1] = A^2      mating[2,1] = 2*A*(1-A)      mating[3,1] = (1-A)^2      nt.y1 = mating %*% leslie.y %*% nt.y      nt.y = nt.y1      allele_freq[ii]=(0.5*nt.y[2]+nt.y[3])/(nt.y[1]+nt.y[2]+nt.y[3])
      HET[j,ii] = 2*allele_freq[ii]*(1 - allele_freq[ii])      total = nt.y[1]+nt.y[2]+ nt.y[3]      total = round(total)      x=rmultinom(1,NewBorn.Size,c(round(nt.y[1])/total,round(nt.y[2])/total,round(nt.y[3])/total))      x = x/NewBorn.Size      # total maintains stable population size      nt.y[1] = x[1] * NewBorn.Size      nt.y[2] = x[2] * NewBorn.Size      nt.y[3] = x[3] * NewBorn.Size
    }
    
  }
  
  for (i in 1:T){
    mean.HET[i] = mean(HET[1:iteration,i])
  }

  lines(mean.HET,col="black",lwd=1,lty=1)
  print(mean.HET)
	
# output the final heterozygosity values for each simulation

het.Dist = HET[1:iteration,T]
print(cat(het.Dist,sep="\n"))
}  