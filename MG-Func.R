##### build-in functions for experiments of M-G time series

genmfs<-function(ID,fis,Type){
  fis = addmf(fis, Type, ID, "S3", "trimf", 
              c(0.0,0.25,0.5))
  fis = addmf(fis, Type, ID, "S2", "trimf", 
              c(0.25,0.5,0.75))
  fis = addmf(fis, Type, ID, "S1", "trimf", 
              c(0.5,0.75,1))
  fis = addmf(fis, Type, ID, "CE", "trimf", 
              c(0.75,1,1.25))
  fis = addmf(fis, Type, ID, "B1", "trimf", 
              c(1,1.25,1.5))
  fis = addmf(fis, Type, ID, "B2", "trimf", 
              c(1.25,1.5,1.75))
  fis = addmf(fis, Type, ID, "B3", "trimf", 
              c(1.5, 1.75, 2))
  fis
}
Pairs<-function(Data,n){
  pair<-matrix(nrow=length(Data)-n,ncol = n+1)
  for (i in 1:(length(Data)-n)) {
    pair[i,1]<-Data[i]
    pair[i,2]<-Data[i+1]
    pair[i,3]<-Data[i+2]
    pair[i,4]<-Data[i+3]
    pair[i,5]<-Data[i+4]
    pair[i,6]<-Data[i+5]
  }
  return(pair)
}

Pairs2<-function(Data,n){
  pair<-matrix(nrow=length(Data)-n,ncol = n+1)
  for (i in 1:(length(Data)-n)) {
    pair[i,1]<-Data[i]
    pair[i,2]<-Data[i+1]
    pair[i,3]<-Data[i+2]
    pair[i,4]<-Data[i+3]
    pair[i,5]<-Data[i+4]
    pair[i,6]<-Data[i+5]
    pair[i,7]<-Data[i+6]
    pair[i,8]<-Data[i+7]
    pair[i,9]<-Data[i+8]
    pair[i,10]<-Data[i+9]
  }
  return(pair)
}

###generate whole rule list
rulelist<-function(Dataset){
  rules<-Dataset
  degree<-Dataset
  for (i in 1:nrow(Dataset)) {
    for (j in 1:ncol(Dataset)) {
      if(Dataset[i,j]<=0.375){
        rules[i,j]<-1
      }
      else if(Dataset[i,j]<=0.625){
        rules[i,j]<-2
      }
      else if(Dataset[i,j]<=0.875){
        rules[i,j]<-3
      }
      else if(Dataset[i,j]<=1.125){
        rules[i,j]<-4
      }
      else if(Dataset[i,j]<=1.375){
        rules[i,j]<-5
      }
      else if(Dataset[i,j]<=1.625){
        rules[i,j]<-6
      }
      else{
        rules[i,j]<-7
      }
      degree[i,j]<-evalmf(Dataset[i,j],fis$input[[1]]$mf[[rules[i,j]]]$type,
                          fis$input[[1]]$mf[[rules[i,j]]]$params)
    }
  }
  Deg<-apply(degree, 1, prod)
  rl<-cbind(rules,Deg)
  return(rl)
}


###keep the rules with highest degree
relist<-function(rl){
  dup<-duplicated(rl[,1:5])
  RL<-matrix(nrow = 0,ncol = 7)
  for (i in 1:nrow(rl)) {
    if(dup[i]==FALSE){
      RL<-rbind(RL,rl[i,])
    }
  }
  RuleList<-RL
  for (i in 1:nrow(RL)) {
    IP<-RL[i,1:5]
    for (j in 1:nrow(rl)) {
      IP2<-rl[j,1:5]
      if(identical(IP,IP2)){
        if(rl[j,7]>RL[i,7]){
          RL[i,6]<-rl[j,6]
          RL[i,7]<-rl[j,7]
        }
      }
    }
  }
  return(RL)
}


relist2<-function(rl){
  dup<-duplicated(rl[,1:9])
  RL<-matrix(nrow = 0,ncol = 11)
  for (i in 1:nrow(rl)) {
    if(dup[i]==FALSE){
      RL<-rbind(RL,rl[i,])
    }
  }
  RuleList<-RL
  for (i in 1:nrow(RL)) {
    IP<-RL[i,1:9]
    for (j in 1:nrow(rl)) {
      IP2<-rl[j,1:9]
      if(identical(IP,IP2)){
        if(rl[j,11]>RL[i,11]){
          RL[i,10]<-rl[j,10]
          RL[i,11]<-rl[j,11]
        }
      }
    }
  }
  return(RL)
}

tstFIS<-function(){
  fis = newfis("MG")
  
  fis = addvar(fis, "input", "T-8", c(0.0, 2))
  fis = addvar(fis, "input", "T-7", c(0.0, 2))
  fis = addvar(fis, "input", "T-6", c(0.0, 2))
  fis = addvar(fis, "input", "T-5", c(0.0, 2))
  fis = addvar(fis, "input", "T-4", c(0.0, 2))
  fis = addvar(fis, "input", "T-3", c(0.0, 2))
  fis = addvar(fis, "input", "T-2", c(0.0, 2))
  fis = addvar(fis, "input", "T-1", c(0.0, 2))
  fis = addvar(fis, "input", "T", c(0.0, 2))
  
  fis = addvar(fis, "output", "T+1", c(0.0, 2))
  
  fis=genmfs(1,fis,'input')
  plotmf(fis,"input",1,main = "T-8")
  fis=genmfs(2,fis,'input')
  plotmf(fis,"input",2,main = "T-7")
  fis=genmfs(3,fis,'input')
  plotmf(fis,"input",3,main = "T-6")
  fis=genmfs(4,fis,'input')
  plotmf(fis,"input",4,main = "T-5")
  fis=genmfs(5,fis,'input')
  plotmf(fis,"input",5,main = "T-4")
  fis=genmfs(6,fis,'input')
  plotmf(fis,"input",6,main = "T-3")
  fis=genmfs(7,fis,'input')
  plotmf(fis,"input",7,main = "T-2")
  fis=genmfs(8,fis,'input')
  plotmf(fis,"input",8,main = "T-1")
  fis=genmfs(9,fis,'input')
  plotmf(fis,"input",9,main = "T")
  
  
  fis=genmfs(1,fis,'output')
  plotmf(fis,"output",1,main = "T+1")
  fis
}

tstFIS2<-function(){
  fis = newfis("MG")
  
  fis = addvar(fis, "input", "T-8", c(0.0, 1.8))
  fis = addvar(fis, "input", "T-7", c(0.0, 1.8))
  fis = addvar(fis, "input", "T-6", c(0.0, 1.8))
  fis = addvar(fis, "input", "T-5", c(0.0, 1.8))
  fis = addvar(fis, "input", "T-4", c(0.0, 1.8))
  
  fis = addvar(fis, "output", "T+1", c(0.0, 1.8))
  
  fis=genmfs(1,fis,'input')
  plotmf(fis,"input",1,main = "T-4")
  fis=genmfs(2,fis,'input')
  plotmf(fis,"input",2,main = "T-3")
  fis=genmfs(3,fis,'input')
  plotmf(fis,"input",3,main = "T-2")
  fis=genmfs(4,fis,'input')
  plotmf(fis,"input",4,main = "T-1")
  fis=genmfs(5,fis,'input')
  plotmf(fis,"input",5,main = "T")
  
  
  fis=genmfs(1,fis,'output')
  plotmf(fis,"output",1,main = "T+1")
  fis
}


try<-function(n){
  mse<-matrix(nrow = 0,ncol = 6)
  for (i in 1:n) {
    tst<-CompSN()
    mse<-rbind(mse,tst)
    print(mse)
  }
  return(mse)
}

try3<-function(n){
  mse<-matrix(nrow = 0,ncol = 5)
  seed=31
  for (i in 1:n) {
    print(seed)
    tst<-Compnsfis(seed)
    mse<-rbind(mse,tst)
    print(mse)
    seed=seed+1
  }
  return(mse)
}

MSE<-function(a,b){
  sum((a-b)^2)/length(a)
}