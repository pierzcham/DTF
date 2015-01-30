dtf <- function(eeg,fsamp,plot=FALSE,f=0,...) {
    eeg <- as.matrix(eeg)
    B <- length(eeg[1,])
    N <- length(eeg[,1])
    eeg.ar <- ar(eeg)
    p <- eeg.ar$order
    test <- matrix(eeg.ar$ar,(p*B),B,byrow=FALSE)
    A <- list()
    for(i in 1:p){
        A[i] <- list(matrix(test[seq(i,(p*B-p+i),p),],B,B))
    }
    s3=matrix(0,(p+1),N)
    for(n in 0:p){
        for(k in 0:(N-1)){
            s3[(n+1),(k+1)]=exp(-1i*2*pi*n*k/N)
        }
    }

    H <- list()
    H0 <- list()
    temp=0
    for(i in 1:N){
        for(j in 1:p){
            temp <- temp + A[[j]]*s3[(j+1),i]
        }
        H0[i] <- list(diag(B)*s3[1,i]-temp)
        H[i] <- list(solve(H0[[i]]))
        temp=0
    }
    
    H00 <- matrix(list(),B,B)
    temp <- rep(NULL,N)
    for(i in 1:B){
        for(j in 1:B){
            for(k in 1:N){
                temp[k] <- H[[k]][i,j]
            }
            H00[i,j] <- list(temp)
        }
    }

    fNyq <- fsamp/2 # Nyquist freq
    Nf=N/2
    nyfreqs = seq(0,fNyq,fsamp/N)
    Se00 <- matrix(list(),B,B)
    for(i in 1:B){
        for(j in 1:B){
            X <- H00[[i,j]][1:(Nf+1)]# transformata jest symetryczna 
            Sa <- Mod(X) # Amplitude spectrum
            Se <- 1/(fsamp*N)*Sa^2
            Se[2:Nf] <- Se[2:Nf]*2 
            Se00[i,j] <- list(Se)
        }
    }
    
    S1 <- rep(0,(Nf+1))
    for(i in 1:B){
        for(j in 1:B){
            S1 <- S1+Se00[[i,j]]
        }
        for(j in 1:B){
            Se00[i,j] <- list(Se00[[i,j]]/S1)
        }
        S1 <- rep(0,(Nf+1))
    }

   
    f.length <- length(f)
    S <- Se00
    S.f1 <- matrix(0,B,B)
    S.f2 <- matrix(0,B,B)
    S.f3 <- matrix(0,B,B)
    if(f.length>1){
        for(i in 1:B){
            for(j in 1:B){
                S.zoo <- zoo(S[[i,j]],nyfreqs)
                S.zoo.f <- S.zoo[index(S.zoo)>=f[1]&index(S.zoo)<=f[2]]
                S.f1[i,j] <- sum(coredata(S.zoo.f))/length(coredata(S.zoo.f))
                S.zoo.f <- S.zoo[index(S.zoo)>=f[3]&index(S.zoo)<=f[4]]
                S.f2[i,j] <- sum(coredata(S.zoo.f))/length(coredata(S.zoo.f))
                S.zoo.f <- S.zoo[index(S.zoo)>=f[5]&index(S.zoo)<=f[6]]
                S.f3[i,j] <- sum(coredata(S.zoo.f))/length(coredata(S.zoo.f))
            }
        }
        S.f <- list("S.f1"=S.f1,"S.f2"=S.f2,"S.f3"=S.f3)    
    }
    else{
        S.f <- list("S"=S)
    }

    if(plot){
        par(mfrow=c(B,B),mar=c(0.1,0.1,0.1,0.1),...)
        for(i in 1:B){
            for(j in 1:B){
                plot(nyfreqs,S[[i,j]],type='l',ylim=c(0,1),axes=FALSE)
                polygon(c(nyfreqs,rev(nyfreqs)),c(rep(0,length(S[[i,j]])),rev(S[[i,j]])),col="gray")
            }
        }
    }
    S.f
}

