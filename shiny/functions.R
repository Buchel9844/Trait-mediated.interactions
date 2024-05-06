#To test funcitons
#Amin <- matrix(c(1,1,1,1), ncol=2)
#Aslope <- matrix(c(1,2,3,4), ncol=2)
#C<- matrix(c(1,1,1,1), ncol=2)
#N <-  matrix(c(1,1,1,1), ncol=2)
#No <- matrix(c(1,1,1,1), ncol=2)
#lambda <- c(1,2)
#N0=c(1,2)
# b.alpha=c(0.1,0.2)
# b.lambda=c(0.1,0.1)

sigmoidal.function <- function(Amin, Aslope,C,N,No,b.alpha){
  e = exp(Aslope*(N-No)) # c is stretching the graph horizontally 
  a = C*(1-e) #stretching the graph vertically
  d = Amin + b.alpha
  alpha = (a/(1 + e)) + d
  
  return(alpha)
}


functions.plot <- function(Amin, Aslope, C ,N,No,b.alpha){
  MatA <- NULL
  for(i in N){
  MatA.i <- data.frame(alpha =c(sigmoidal.function(Amin, Aslope,C ,i,No,b.alpha)), 
                       focal = c("i","i","j","j"), 
                       neigh= c("i","j","i","j"),
                       type = c("ii","ij","ji","jj"), 
                       density= i)
  MatA <- bind_rows(MatA,MatA.i)
 }

  plot.interaction <- ggplot(MatA, 
               aes(y=alpha, x= density, color=focal, group=type))+
    labs(title="Interactions effect") + 
    geom_point(alpha=0.5,aes(shape=type,colour = after_stat(y < 0)))+     
    guides(color="none") +
    geom_hline(yintercept=0,linetype="dashed") +
    xlab("Neighbour density") + ylab("per capita effect of neigh on focal")+
    theme_bw()+
    if(all( MatA$alpha< 0)){
      scale_colour_manual(
        values = c("red")) 
    }else{scale_colour_manual(
      values = c("limegreen","red")) }
  plot.interaction
  
}


fecundity.plot <- function(lambda,Amin, Aslope, C ,N, No,b.alpha,b.lambda){
  MatA.df  <- NULL
  for(i in N){
    MatA <- sigmoidal.function(Amin, Aslope,C ,i,No,b.alpha)
    MatA.df.i <- data.frame(fecundity = c((lambda[1]*exp(sum(MatA[1,]*i))*b.lambda[1]),
                                      (lambda[2]*exp(sum(MatA[2,]*i))*b.lambda[2])),
                         focal = c("i","j"), 
                        density= i)
    MatA.df <- bind_rows(MatA.df,MatA.df.i)
  }

  
  fecundity.plot <- MatA.df %>% 
    ggplot(aes(y=fecundity, x= density, 
               color=focal))+
    labs(title="Fecundity distribution") + 
    geom_point() +
    geom_line(aes(color=focal))+
    xlab("Density of neighbour species") + 
    ylab("Fecundity of focal species")+
    theme_bw()
  
  return(fecundity.plot)
}



Abundance.plot<- function(lambda,Amin, Aslope, C , N0,No, time,b.alpha,b.lambda){
    
  Nt<- data.frame(Ni = N0[1],Nj=N0[2], time=1)

  for (t in 1:time){

     matA <-  sigmoidal.function(Amin, Aslope, C ,
                                        matrix(rep(as.numeric(Nt[t,c(1:2)]),each=2),
                                               ncol=2),
                                    No,b.alpha)


    
    Nt[t+1,1] <- Nt[t,1]*(lambda[1]*exp(sum(matA[1,]*Nt[t,c(1,2)]))*b.lambda[1])
    Nt[t+1,2] <- Nt[t,2]*(lambda[2]*exp(sum(matA[2,]*Nt[t,c(2,1)]))*b.lambda[2])
    Nt[t+1,3]<- t+1

    
  }
  
  Nt.plot <-  Nt %>%
    gather(Ni,Nj, key="focal", value="abundance") %>%
    ggplot(aes(y=abundance, x= time, color=focal))+
    labs(title="Abundance distribution") + 
    geom_point() +
    geom_line(aes()) + 
    xlab("Generation") + ylab("Abundance of focal species)")+
    theme_bw()

  return(Nt.plot)
}

