library(FuzzyR)
source("MG-Func.R")


##### different datasets and rule base#####
fis<-tstFIS()

  X<-read.csv('MG.csv')
  ds<-X[1:999,]

  ##NF data set
  dsNF<-Pairs2(ds,9)
  plot(t<-(1:990),dsNF[,10],type = "l")

  ##add noise
  signal<-sd(ds)
  N5dB<-signal/10^(5/20)
  N10dB<-signal/10^(10/20)
  #####trn noise generation#####
  noise10<-runif(999,-sqrt(3)*N10dB,sqrt(3)*N10dB)
  noise5<-runif(999,-sqrt(3)*N5dB,sqrt(3)*N5dB)
  ##add noise
  ds10<-ds+noise10
  ds5<-ds+noise5

  par(mfrow=c(3,1))
  plot(1:999,ds,type = "l",main = "NF data")
  plot(1:999,ds10,type = "l",main = "10dB data")
  plot(1:999,ds5,type = "l",main = "5dB data")
  ##data pair
  ds10<-Pairs2(ds10,9)
  ds5<-Pairs2(ds5,9)

  ##trn and tst data set
  Num<-1:700
  trnNF<-dsNF[Num,]
  trn10<-ds10[Num,]
  trn5<-ds5[Num,]
  Num2<-701:990
  tstNF<-dsNF[Num2,]
  tst10<-ds10[Num2,]
  tst5<-ds5[Num2,]
  ##create rulelist
  rlnf<-rulelist(trnNF)
  rl10<-rulelist(trn10)
  rl5<-rulelist(trn5)

  ##filter rules by their degrees
  rlNF<-relist2(rlnf)[,1:10]
  rl10dB<-relist2(rl10)[,1:10]
  rl5dB<-relist2(rl5)[,1:10]

  ###testing data
  ##noise free testing set(ground truth)

  #par(mfrow=c(3,1))
  #plot(t<-(701:990),tstNF[,10],type = "l",main = "NF testing data")
  #plot(t<-(701:990),tst10[,10],type = "l",main = "10dB testing data")
  #plot(t<-(701:990),tst5[,10],type = "l",main = "5dB testing data")


  ####rulelist for both system####
  rl<-rl5dB
  rl<-cbind(rl,matrix(c(1,1),nrow=nrow(rl),ncol = 2))
  print(nrow(rl))

  ##setting input data and fuzzification deviation
  Input<-tst5[,1:9]
  noiseSD<-N5dB



  ######Different NSFLSs######

  
  nsfls<-function(fls,index,fuzzification.method,noiseSD,firing.method){
    for (i in 1:index) {
      fls$input[[i]]$fuzzification.method<-fuzzification.method
      fls$input[[i]]$fuzzification.params<-noiseSD
      fls$input[[i]]$firing.method<-firing.method
    }
    return(fls)
  }

  
  par(mfrow=c(1,1))
  STDfis<-tstFIS()
  STDfis<-nsfls(STDfis,9,"gauss",noiseSD,"tnorm.min.max")
  STDfis = addrule(STDfis, rl)

  CENfis<-tstFIS()
  CENfis<-nsfls(CENfis,9,"gauss",noiseSD,"tnorm.min.defuzz.centroid")
  CENfis = addrule(CENfis, rl)

  SIMfis<-tstFIS()
  SIMfis<-nsfls(SIMfis,9,"gauss",noiseSD,"similarity.set")
  SIMfis = addrule(SIMfis, rl)


  ##evaluation
  STDNon<-evalfis(Input,STDfis,point_n = 101)
  CENNon<-evalfis(Input,CENfis)
  SIMNon<-evalfis(Input,SIMfis)

  t<-795:951
  plot(t,ds[795:951],type="l",xlab="Time (t)",ylab="M-G Values")
  #points(t,tstNF[86:242,10],type = "l",col="tomato")
  points(t,ds5[795:951],type="l",lty=2,col="brown")
  points(t,STDNon[86:242],type = "l",col="tomato",lwd=2)
  points(t,CENNon[86:242],type = "l",col="forestgreen",lwd=2)
  points(t,SIMNon[86:242],type = "l",col="blue4",lwd=2)
  legend("topleft", pch=16, legend =c("Noise Free","5dB Noisy","Sim-NS","Cen-NS","Sta-NS"),
         col =c("black","brown","blue4","forestgreen","tomato"), cex = 0.7)
  title(main = "Noisy training, Noisy input (NSR=5dB)",cex.main=1.1)
 par(mfrow=c(1,1))

