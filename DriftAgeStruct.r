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

leslie.drift = function(infile,B,N,T,iteration,initial,recessive=T,s){
# read in lx and mx data

  het.Dist=matrix(0,iteration,T) # heterozygosity at end time point
 
# B is specified starting allele frequency 

  A = 1 - B


# constructs Px matrix 


    leslie[1,i]= fx[i]
  }
    leslie[i,(i-1)]= Px[i-1]
  }
# Let population age-structure stabilise to set nt.0
    nt.0[i] = nt.0[i]/total
  }
  
  print(nt.0)

# A barplot can be produced:
  barplot(nt.0,beside=T,col="light blue",names=as.character(1:Nx),main="Age-Structure of Population",xlab="Age-class(x)",ylab="Frequency")

# If so, then put break in to pause simulation

  readline(prompt="Press [enter] to continue")

# Prepare figure for simulations

  plot(dummyX,type="n",xlab="Time",ylab="Mean Heterozygosity",cex.lab=1.5, cex.axis=1.5,xlim=c(0,T),ylim=c(0,0.5))
  
# Now introduce genotypes and drift to leslie matrix
  for(j in 1:iteration){  # ds$iteration specifies number of generations
    A = 1 - B
    Ny = 3 * Nx	#treble matrix to allow for each age class to have three 
      nt.y[i]=A^2*nt.y[i]
    }
      nt.y[i]=2*A*B*nt.y[i]
    }
      nt.y[i]=B^2*nt.y[i]
    }
    k = 1
    for(i in 1:Nx){

      if(ds$mx[i]>0){
        Px.y[k+2] = (1 - s) * Px.y[k+2]
      }
      k = k+3
    }

    #Assign each fecundity (fx) to a genotype
      mating[i,i]=1
    }
      leslie.y[1,i]=fx.y[i]
    }
      leslie.y[i,(i-3)]=Px.y[(i-3)]
    }
      for(iii in seq(4,by=3,(Nx-1))){
          x=rmultinom(1,total,c(round(nt.y[iii])/total,round(nt.y[iii+1])/total,round(nt.y[iii+2])/total))
        } 

        p[i]=fx.y[i]*nt.y[i]
      }
        hom[i]=p[i]
      }
        het[i]=p[i]
      }
        A = 1
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