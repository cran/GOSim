normalize.kernel = function(Ker, method="none"){
	if(method != "none"){
		if(method == "sqrt"){ # result between -1 and 1
			Kd<-sqrt(diag(Ker) + 1e-10)
			Ker<-Ker/(Kd%*%t(Kd))			
		}
		else if(method == "Lin"){ # result: diagonal = 1		
			Kd = diag(Ker)
			Ker = 2*Ker / outer(Kd, Kd, "+")
		}
		else if(method == "Tanimoto"){ 
			Kd = diag(Ker)
			Ker = Ker / (outer(Kd, Kd, "+") - Ker)
		}								
		else
			stop(paste("Unknown normalization method", method))
		diag(Ker) = 1
	}
	Ker
}
