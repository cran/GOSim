normalize.kernel = function(Ker, method="none"){
	if(method != "none"){
		if(method == "sqrt"){ # result between -1 and 1
			Kd<-sqrt(diag(Ker) + 1e-10)
			Ker<-Ker/(Kd%*%t(Kd))			
		}
		else if(method == "Lin"){ # result: diagonal = 1				
			for(i in 1:ncol(Ker)){					
				if(i > 1){
					for(j in 1:(i-1)){
						Ker[i,j] = 2*Ker[i,j] / (Ker[i,i] + Ker[j,j])
						Ker[j,i] = Ker[i,j]
					}
				}
			}			
			diag(Ker) = 1	
		}
		else if(method == "Tanimoto"){ 
			for(i in 1:ncol(Ker)){					
				if(i > 1){
					for(j in 1:(i-1)){
						Ker[i,j] = Ker[i,j] / (Ker[i,i] + Ker[j,j] - Ker[i,j])
						Ker[j,i] = Ker[i,j]
					}
				}
			}			
			diag(Ker) = 1	
		}								
		else
			stop(paste("Unknown normalization method", method))
	}
	Ker
}
