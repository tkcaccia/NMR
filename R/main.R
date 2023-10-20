
read_multiple_address = function(read){
  addr=list()
  addr$n=0
  for(i in 1:nrow(read)){
    Address=lapply(read,function(x) x[i])
    names(Address)=colnames(read)
    temp=read_address(Address)
    addr$n=temp$n+addr$n
    addr$Names= c(addr$Names,temp$Names)
    addr$ID= c(addr$ID,rep(temp$ID,temp$n))
    addr$Source= c(addr$Source,temp$Source)
    addr$loc= rbind(addr$loc,temp$loc)
    addr$deleted= c(addr$deleted,temp$deleted)
    addr$FID= c(addr$FID,temp$FID)
  }  
  addr
}
reset_allignment_spectra = 
  function(Address,reset=14.75){
    name=0
    read=read_address(Address,name)
    SNames=read$Names
    Experiments=read$Source
    
    rbnmr <- function (i)
    { 
      setwd (Experiments[i])
      con_proc=file("procs")
      Procs <- readLines(con_proc, n = -1)
      pma=pmatch("##$OFFSET=",Procs)
      OFFSET <- as.numeric(gsub(".*[=]", "",Procs[pma]))
      Procs[pma]=paste("##$OFFSET=",reset)
      writeLines(Procs,con=con_proc)
      i=i+1
      close(con_proc)
      
    }
    Spec <- lapply (1:length(Experiments),rbnmr)
    
  }

write_spectra =
  function(spectra,Address,name=0){
    read=read_address(Address,name)
    SNames=read$Names
    Experiments=read$Source
    loc=read$loc
    deleted=read$deleted
    
    
    rbnmr <- function (i){ 
      
      setwd (Experiments[i])
      
      con_proc=file("procs")
      
      Procs <- readLines(con_proc, n = -1)
      BYTORD <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$BYTORDP=",Procs)]))
      if (BYTORD=="0") {ENDIAN="little"} else {ENDIAN="big"}
      NC_proc <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$NC_proc=",Procs)])) 
      OFFSET <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$OFFSET=",Procs)]))
      SW_p <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$SW_p=",Procs)]))
      SF <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$SF=",Procs)]))
      SI <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$SI=",Procs)]))
      
      DataR = NULL
      yy=as.numeric(spectra$Spectra[[i]]$y)
      yy=as.integer(round(yy*(2^(-NC_proc))))
      close(con_proc)
      writeBin (yy,"1r",endian="little")
    }
    Spec <- lapply (1:length(Experiments),rbnmr)
    
  }




matrix_creation = function(addr,right=0,left=10,delta=0.0001,split=1000,span=0.07){
  split_size=(left-right)/split
  half_split_size=split_size/2
  new_x=seq(right,left,delta)
  len_new_x=length(new_x)
  allign=matrix(nrow=addr$n,ncol=len_new_x)
  colnames(allign)=as.character(new_x)
  rownames(allign)=addr$Names
  for(k in 1:addr$n){
    spectrum=addr$Source[k]
    spectra=read_spectra(spectrum)
    y0=spectra$Spectra[[1]]$y
    x0=spectra$Spectra[[1]]$x
    z=new_x
    for(j in 1:split){
      deltaj=(j-1)*split_size
      right2=right+deltaj
      left2=right+deltaj+split_size
      right1=right2-half_split_size
      left1=left2+half_split_size
      sel1=x0>right1 & x0<left1
      x1=x0[sel1]
      y1=y0[sel1]
      
      sel2=new_x>=right2 & new_x<left2
      
      z[sel2]=predict(loess(y1~x1,span = span,control=loess.control(surface="direct")),new_x[sel2])
    }
    z[len_new_x]=predict(loess(y1~x1,span = span,control=loess.control(surface="direct")),new_x[len_new_x])
    allign[k,]=z
    print(k)
  } 
  allign
}


read_spectra =
  function(Address,name=0){
    
    read=read_address(Address,name)
    SNames=read$Names
    Experiments=read$Source
    FIDs=read$FID
    loc=read$loc
    deleted=read$deleted
    rbnmr <- function (i){
      DataR = NULL
      setwd (FIDs[i])
      con_acqu=file("acqus")
      Acqus <- readLines(con_acqu, n = -1)
      DATE=as.numeric(gsub(".*[=]", "",Acqus[pmatch("##$DATE= ",Acqus)]))
      NS <- as.numeric(gsub(".*[=]", "",Acqus[pmatch("##$NS=",Acqus)]))
      
      class(DATE) = c('POSIXt','POSIXct')
      
      setwd (Experiments[i])
      con_proc=file("procs")
      Procs <- readLines(con_proc, n = -1)
      BYTORD <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$BYTORDP=",Procs)]))
      
      if (BYTORD=="0") {ENDIAN="little"} else {ENDIAN="big"}
      NC_proc <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$NC_proc=",Procs)]))
      OFFSET <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$OFFSET=",Procs)]))
      SW_p <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$SW_p=",Procs)]))
      SF <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$SF=",Procs)]))
      SI <- as.numeric(gsub(".*[=]", "",Procs[pmatch("##$SI=",Procs)]))
      con = file("1r", open="rb")
      
      DataR$y <- readBin (con,what="int",n=SI,endian=ENDIAN)
      close (con)
      close(con_proc)
      close(con_acqu)
      DataR$x <- seq (OFFSET,OFFSET-SW_p/SF,length=SI)
      DataR$y <- DataR$y/(2^(-NC_proc))
      DataR$date=DATE
      DataR$ns=NS
      DataR$SF=SF
      DataR$delta=(SW_p/SF)/(SI-1)
      DataR
    }
    Spec <- lapply (1:length(Experiments),rbnmr)
    return(list(Spectra=Spec, 
                Names=SNames, 
                n=length(Experiments),
                Source=Experiments,
                loc=loc,
                deleted=deleted))
    
  }






read_address =
  function(Address,nick=0){
    name=nick
    if(is.list(Address)){
      if(Address$Partition=="")   stop("You have to declare the partition (e.g., C:, D:)")
      DIR <- paste(Address$Partition,"/",sep="")
      SNames=NULL
      
      if(Address$Data_directory!=""){
        if(Address$Data_directory=="*"){
          if(nick==0){
            name=1
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_1=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_1,"/",sep="")
            
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_1
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_1,sep="")
            }
            
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_1=Address$Data_directory
          DIR=paste(DIR,Address$Data_directory,"/",sep="")
          if(any(name==1)){
            if(is.null(SNames)){
              SNames=DIR_1
            }else{
              SNames=paste(SNames,"_",DIR_1,sep="")
            }
          }   
        }
      }
      
      if(Address$User_name!=""){
        if(Address$User_name=="*"){
          if(nick==0){
            name=2
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_2=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_2,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_2
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_2,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          
          DIR_2=Address$User_name
          DIR=paste(DIR,Address$User_name,"/",sep="")
          if(any(name==2)){
            if(is.null(SNames)){
              SNames=DIR_2
            }else{
              SNames=paste(SNames,"_",DIR_2,sep="")
            }
          }   
        }
      }
      
      
      if(Address$Spectroscopy!=""){
        if(Address$Spectroscopy=="*"){
          if(nick==0){
            name=3
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_3=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_3,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_3
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_3,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_3=Address$Spectroscopy
          DIR=paste(DIR,Address$Spectroscopy,"/",sep="")
          if(any(name==3)){
            if(is.null(SNames)){
              SNames=DIR_3
            }else{
              SNames=paste(SNames,"_",DIR_3,sep="")
            }
          }   
        }
      }
      
      if(Address$Dataset_name!=""){
        if(Address$Dataset_name=="*"){
          if(nick==0){
            name=4
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_4=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_4,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_4
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_4,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_4=Address$Dataset_name
          DIR=paste(DIR,Address$Dataset_name,"/",sep="")
          if(any(name==4)){
            if(is.null(SNames)){
              SNames=DIR_4
            }else{
              SNames=paste(SNames,"_",DIR_4,sep="")
            }
          }   
        }
      }
      
      if(Address$Esperiment_number!=""){
        if(Address$Esperiment_number=="*"){
          if(nick==0){
            name=5
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_5=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_5,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_5
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_5,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_5=Address$Esperiment_number
          DIR=paste(DIR,Address$Esperiment_number,"/",sep="")
          if(any(name==5)){
            if(is.null(SNames)){
              SNames=DIR_5
            }else{
              SNames=paste(SNames,"_",DIR_5,sep="")
            }
          }   
        }
      }
      FID <- DIR
      if(Address$Pdata!=""){
        if(Address$Pdata=="*"){
          if(nick==0){
            name=6
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_6=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_6,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_6
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_6,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_6=Address$Pdata
          DIR=paste(DIR,Address$Pdata,"/",sep="")
          if(any(name==6)){
            if(is.null(SNames)){
              SNames=DIR_6
            }else{
              SNames=paste(SNames,"_",DIR_6,sep="")
            }
          }   
        }
      }
      
      if(Address$Processing_number!=""){
        if(Address$Processing_number=="*"){
          if(nick==0){
            name=7
          }
          temp1=list()
          temp2=list()
          for(i in 1:length(DIR)){
            DIR_7=dir(DIR[i])
            temp1[[i]]=paste(DIR[i],DIR_7,"/",sep="")
            if(is.null(SNames[i])){
              temp2[[i]]=DIR_7
            }else{
              temp2[[i]]=paste(SNames[i],"_",DIR_7,sep="")
            }
          }
          DIR=unlist(temp1)
          SNames=unlist(temp2)
        }else{
          DIR_7=Address$Processing_number
          DIR=paste(DIR,Address$Processing_number,"/",sep="")
          if(any(name==7)){
            if(is.null(SNames)){
              SNames=DIR_7
            }else{
              SNames=paste(SNames,"_",DIR_7,sep="")
            }
          }   
        }
        
      }
    }else{
      DIR=Address
      naddr=length(Address)
      str=strsplit(Address,"/")
      
      if(nick==0){
        SNames=Address
      }else{
        
        SName=NULL
        for(jj in 1:naddr){
          
          SName[jj]=paste(str[[jj]][nick+1],collapse = "_")
          
        }
      }
      ii_fid=1
      while(ii_fid<length(str[[1]]) & !file.exists(paste(paste(str[[1]][1:ii_fid],collapse = "/"),"/fid",sep=""))){
        ii_fid=ii_fid+1
      }
      ii_ser=1
      while(ii_ser<length(str[[1]]) & !file.exists(paste(paste(str[[1]][1:ii_ser],collapse = "/"),"/ser",sep=""))){
        ii_ser=ii_ser+1
      }
      ii=min(c(ii_fid,ii_ser))
      FID=NULL
      for(jj in 1:naddr){
        
        FID[jj]=paste(paste(str[[jj]][1:ii],collapse = "/"),"/",sep="")
      }
      FID
    }
    
    
    
    Experiments <- DIR
    
    
    sel=file.exists(paste(DIR,"1r",sep=""))
    deleted=Experiments[!sel]
    if(any(sel)){
      Experiments <- Experiments[sel]
      FID <- FID[sel]
      SNames = SNames[sel]
      aa=strsplit(Experiments,"/")
      loc=t(matrix(unlist(aa),ncol=length(aa)))
      Spec=NULL
      return(list(Names=SNames, n=length(Experiments),Source=Experiments,loc=loc,deleted=deleted,FID=FID))
    }else{
      return(list(Names=NULL, n=0,Source=NULL,loc=NULL,deleted=deleted,FID=NULL))
    }
  }





















opt_bucket <- function (X,size_bucket, slackness=0) {
  xaxis=as.numeric(colnames(spetri_allineati_NOESY))
  p = nrow(X)
  q = ncol(X)
  b = abs(xaxis[1]-xaxis[2])
  a = size_bucket/b
  a = round(a)
  l = slackness*a
  l = round(l)
  R = colMeans(X)
  v = NULL
  
  for (t in seq(1+a,q-a, a)) {
    I = which.min( R[(t-l):(t+l)] )
    f = ((t-(1+a))/a)+1
    v[f] = I+a*(f-1)+(a-l)
  }
  
  z = unique(v)
  v = c(1, z, q)
  Z = matrix(0, nrow=p, ncol=(length(v)-1) )
  
  for (j in 1:p) {
    for (n in 1:(length(v)-1)) {
      x_prov = X[j,v[n]:v[n+1]]
      Z[j,n] = trapezq(x_prov)
    }
  }
  
  vv = matrix(0, nrow=p, ncol=length(v))
  
  for (j in 1:p) {
    vv[j,] = X[j,v]
  }
  
  for (tt in 2:(length(v)-1)) {
    Z[,tt] = Z[,tt] - vv[,tt]
  }
  
  ZNN = Z
  
  Z = 100*Z/rowSums(abs(Z))
  
  A = NULL
  for (k in 1:length(v)) {
    A[k]=xaxis[v[k]]
  }
  
  s = length(A)
  I_b = cbind(A[1:(s-1)], A[2:s])
  
  for (k in 1:(s-1)) {
    T[k] = A[k]-A[k+1]
    S_b = T
  }
  
  return (list(ZNN=ZNN,Z=Z,I_b=I_b,S_b=S_b))
}

trapezq <- function (y) {
  n = length(y)
  sums = y[1] + 2*sum(y[2:(n-1)]) + y[n]
  q = sums/2
  return (q)
}





write_1r =
  function(spectra,address,NC_proc){
    txt=paste(address,"/1r",sep="")
    yy=as.numeric(spectra$Spectra[[1]]$y)
    yy=as.integer(round(yy*(2^(-NC_proc))))
    writeBin (yy,txt,endian="little")
  }


NC_proc.calculation = function(spectra){
  yy=as.numeric(spectra$Spectra[[1]]$y)
  NC_proc=-ceiling(log2((2^30)/max(yy)))
  NC_proc
}

SI.calculation = function(spectra){
  SI=length(spectra$Spectra[[1]]$y)
  SI
}

par.calculation = function(spectra,SF){
  xx=as.numeric(spectra$Spectra[[1]]$x)
  rr=c(xx[1],xx[length(xx)])
  OFFSET=max(c(xx[1],xx[length(xx)]))
  ratio=abs(xx[1]-xx[length(xx)])
  SW_p=SF*ratio
  return(list(OFFSET=OFFSET,SF=SF,SW_p=SW_p))
}



write_proc = function(address,SI,NC_proc,OFFSET,SW_p,SF){
  txt=paste(address,"/proc",sep="")
  fileConn<-file(txt)
  writeLines(c("##TITLE= Parameter file, TopSpin 4.1.3",
               "##JCAMPDX= 5.0",
               "##DATATYPE= Parameter Values",
               "##NPOINTS= 1	$$ modification sequence number",
               "##ORIGIN= Bruker BioSpin GmbH",
               "##OWNER= nmruser",
               "##$ABSF1= 1000",
               "##$ABSF2= -1000",
               "##$ABSG= 1",
               "##$ABSL= 3",
               "##$ALPHA= 0",
               "##$AQORDER= 0",
               "##$ASSFAC= 0",
               "##$ASSFACI= 0",
               "##$ASSFACX= 0",
               "##$ASSWID= 0",
               "##$AUNMP= <proc_jres>",
               "##$AXLEFT= 0",
               "##$AXNAME= <F2>",
               "##$AXNUC= <off>",
               "##$AXRIGHT= 0",
               "##$AXTYPE= 0",
               "##$AXUNIT= <>",
               "##$AZFE= 0.1",
               "##$AZFW= 0.1",
               "##$BCFW= 0.2",
               "##$BC_mod= 4","##$BYTORDP= 0","##$COROFFS= 0","##$CY= 15","##$DATMOD= 1","##$DC= 2",
               "##$DFILT= <>","##$DTYPP= 0","##$ERETIC= no","##$F1P= 0","##$F2P= 0","##$FCOR= 0.5",
               "##$FTSIZE= 1024","##$FT_mod= 4","##$GAMMA= 1","##$GB= 0","##$INTBC= 1","##$INTSCL= 1",
               "##$ISEN= 128","##$LB= 0.3","##$LEV0= 0","##$LPBIN= 0","##$MAXI= 10000","##$MC2= 0",
               "##$MEAN= 0","##$ME_mod= 0","##$MI= 0",
               "##$NCOEF= 0",
               paste("##$NC_proc=",NC_proc),
               "##$NLEV= 6","##$NOISF1= 0","##$NOISF2= 0","##$NSP= 1","##$NTH_PI= 0",
               "##$NZP= 0",
               paste("##$OFFSET=",OFFSET),
               "##$PC= 1","##$PHC0= 0","##$PHC1= 0","##$PH_mod= 0","##$PKNL= yes","##$PPARMOD= 0",
               "##$PPDIAG= 0","##$PPIPTYP= 0","##$PPMPNUM= 100","##$PPRESOL= 1","##$PSCAL= 4","##$PSIGN= 0",
               "##$PYNMP= <proc.py>","##$RDF1= 0","##$RDF2= 0","##$RDFWHM= 0","##$RDINT= 0","##$RDPOS= 0",
               "##$REVERSE= no",
               paste("##$SF=",SF),
               paste("##$SI=",SI),
               "##$SIGF1= 0","##$SIGF2= 0","##$SINO= 400",
               "##$SIOLD= 16384","##$SPECTYP= <>","##$SREF_mod= 1","##$SREGLST= <1H.SERUM_OR_PLASMA>",
               "##$SSB= 0","##$STSI= 0","##$STSR= 0",
               paste("##$SW_p=",SW_p),
               "##$SYMM= 0","##$S_DEV= 0","##$TDeff= 0",
               "##$TDoff= 0","##$TI= <Title>","##$TILT= no","##$TM1= 0","##$TM2= 0","##$TOPLEV= 0","##$USERP1= <user>",
               "##$USERP2= <user>","##$USERP3= <user>","##$USERP4= <user>","##$USERP5= <user>","##$WDW= 3",
               "##$XDIM= 8192","##$YMAX_p= 0","##$YMIN_p= 0","##END="), fileConn)
  close(fileConn)
  
}



write_procs = function(address,SI,NC_proc,OFFSET,SW_p,SF){
  txt=paste(address,"/procs",sep="")
  fileConn<-file(txt)
  writeLines(c("##TITLE= Parameter file, TopSpin 4.1.3","##JCAMPDX= 5.0","##DATATYPE= Parameter Values","##NPOINTS= 3	$$ modification sequence number",
               "##ORIGIN= Bruker BioSpin GmbH","##OWNER= nmruser","##$ABSF1= 13.00973","##$ABSF2= -3.655826","##$ABSG= 1","##$ABSL= 3","##$ALPHA= 0",
               "##$AQORDER= 0","##$ASSFAC= 0","##$ASSFACI= 0","##$ASSFACX= 0","##$ASSWID= 0","##$AUNMP= <proc_jres>","##$AXLEFT= 0","##$AXNAME= <F2>",
               "##$AXNUC= <1H>","##$AXRIGHT= 0","##$AXTYPE= 0","##$AXUNIT= <>","##$AZFE= 0.1","##$AZFW= 0.1","##$BCFW= 0","##$BC_mod= 4","##$BYTORDP= 0",
               "##$COROFFS= 0","##$CY= 15","##$DATMOD= 1","##$DC= 0","##$DFILT= <>","##$DTYPP= 0","##$ERETIC= no","##$F1P= 0","##$F2P= 0","##$FCOR= 1",
               "##$FTSIZE= 16384","##$FT_mod= 6","##$GAMMA= 0","##$GB= 0","##$INTBC= 1","##$INTSCL= 1","##$ISEN= 128","##$LB= 0","##$LEV0= 0","##$LPBIN= 0",
               "##$MAXI= 10000","##$MC2= 0","##$MEAN= 0","##$ME_mod= 0","##$MI= 0",
               "##$NCOEF= 0",
               paste("##$NC_proc=",NC_proc),
               "##$NLEV= 6","##$NOISF1= 0",
               "##$NOISF2= 0","##$NSP= 0","##$NTH_PI= 0","##$NZP= 0",
               paste("##$OFFSET=",OFFSET),
               "##$PC= 1","##$PHC0= 0","##$PHC1= 0","##$PH_mod= 0",
               "##$PKNL= yes","##$PPARMOD= 0","##$PPDIAG= 0","##$PPIPTYP= 0","##$PPMPNUM= 100","##$PPRESOL= 1","##$PSCAL= 4","##$PSIGN= 0",
               "##$PYNMP= <proc.py>","##$RDF1= 0","##$RDF2= 0","##$RDFWHM= 0","##$RDINT= 0","##$RDPOS= 0","##$REVERSE= no",
               paste("##$SF=",SF),
               paste("##$SI=",SI),
               "##$SIGF1= 13.04768","##$SIGF2= -3.617876","##$SINO= 400","##$SIOLD= 1024","##$SPECTYP= <>","##$SREF_mod= 1",
               "##$SREGLST= <1H.H2O+D2O>","##$SSB= 0","##$STSI= 16384","##$STSR= 0",
               paste("##$SW_p=",SW_p),
               "##$SYMM= 3","##$S_DEV= 0","##$TDeff= 12288",
               "##$TDoff= 0","##$TI= <1	10680>","##$TILT= yes","##$TM1= 0","##$TM2= 0","##$TOPLEV= 0","##$USERP1= <user>","##$USERP2= <user>",
               "##$USERP3= <user>","##$USERP4= <user>","##$USERP5= <user>","##$WDW= 3","##$XDIM= 0","##$YMAX_p= 128479952","##$YMIN_p= -34434","##END="
  ), fileConn)
  close(fileConn)
  
}





write_used_from = function(address){
  txt=paste(address,"/used_from",sep="")
  fileConn<-file(txt)
  writeLines(c("##TITLE= Parameter file, TopSpin 4.1.3","##JCAMPDX= 5.0","##DATATYPE= Parameter Values",
               "##NPOINTS= 1	$$ modification sequence number","##ORIGIN= Bruker BioSpin GmbH","##OWNER= nmruser",
               "##$CUREXP= <row>","##$DATPATH= <Y:/AMIFLORENCE/nmr>","##$DETAILS= <Qui ci mettiamo un bel titolo>",
               "##$EXPNO= 2","##$NAME= <AMI-0001-3>","##$PROCNO= 1","##$PROCNO2= 1","##$PROCNO3= 256","##END="
  ), fileConn)
  close(fileConn)
  
}


alanine_allignment =
  function(spectra,range=c(1.45,1.55),MHz=500){
    
    
    u=spectra$Spectra[[1]]
    x=u[[2]]
    y=u[[1]]
    sel=x>range[1] & x<range[2]
    x1=x[sel]
    y1=y[sel]
    lx1=length(x1)
    z1=predict(loess(y1~x1,span = 0.05),x1)
    pp=pick.peaks(z1,50)
    
    M=as.matrix(dist(x1[pp]))
    
    W=which(M>(0.0115*600/MHz) & M<(0.01315*600/MHz),arr.ind=T)
    
    
    ss=rowMeans(matrix(x1[pp[W]],ncol=2))
    W=W[ss>(1.475) & ss<(1.52),]  
    
    
    pp=unique(pp[W])
    
    
    ma=matrix(nrow=length(pp),ncol=2)
    ma[,1]=pp-40
    ma[,2]=pp+40
    ma[ma>lx1]=lx1
    ma[ma<1]=1
    
    del=unique(as.numeric(apply(ma,1,function(x) x[1]:x[2])))
    
    x1_del=x1[-del]
    y1_del=y1[-del]
    z2=predict(loess(y1_del~x1_del,span = 0.5),x1)
    
    pp=pp[order(y1[pp]-z2[pp],decreasing=T)[1:2]]
    
    ma=matrix(nrow=length(pp),ncol=2)
    ma[,1]=pp-40
    ma[,2]=pp+40
    ma[ma>lx1]=lx1
    ma[ma<1]=1
    
    del=unique(as.numeric(apply(ma,1,function(x) x[1]:x[2])))
    
    x1_del=x1[-del]
    y1_del=y1[-del]
    z2=predict(loess(y1_del~x1_del,span = 0.5),x1)
    
    peaks_ppm=x1[pp]
    peaks_height=y1[pp]
    
    return(list(baseline=list(x=x1,y=z2),ppm=mean(x1[pp]),peaks=list(ppm=peaks_ppm,height=peaks_height)))
    
    
  }





glucose_allignment =
  function(spectra,range=c(5.20,5.40),MHz=500){
    
    x1=NA
    z2=NA
    peaks_ppm=NA
    peaks_height=NA
    u=spectra$Spectra[[1]]
    x=u[[2]]
    y=u[[1]]
    sel=x>range[1] & x<range[2]
    x1=x[sel]
    y1=y[sel]
    lx1=length(x1)
    z1=predict(loess(y1~x1,span = 0.02),x1)
    pp=pick.peaks(z1,50)
    
    
    
    if(length(pp)>=2){
      
      M=as.matrix(dist(x1[pp]))
      W=which(M>(0.005*600/MHz) & M<(0.007*600/MHz),arr.ind=T)
      
      
      ss=rowMeans(matrix(x1[pp[W]],ncol=2))
      W=W[ss>(5.22) & ss<(5.26),]  
      
      
      
      
      pp=unique(pp[W])
      if(length(pp)>=2){
        
        
        
        ma=matrix(nrow=length(pp),ncol=2)
        ma[,1]=pp-40
        ma[,2]=pp+40
        ma[ma>lx1]=lx1
        ma[ma<1]=1
        
        del=unique(as.numeric(apply(ma,1,function(x) x[1]:x[2])))
        
        x1_del=x1[-del]
        y1_del=y1[-del]
        z2=predict(loess(y1_del~x1_del,span = 0.5),x1)
        
        pp=pp[order(y1[pp]-z2[pp],decreasing=T)[1:2]]
        if(length(pp)>=2){
          ma=matrix(nrow=length(pp),ncol=2)
          ma[,1]=pp-40
          ma[,2]=pp+40
          ma[ma>lx1]=lx1
          ma[ma<1]=1
          
          del=unique(as.numeric(apply(ma,1,function(x) x[1]:x[2])))
          
          x1_del=x1[-del]
          y1_del=y1[-del]
          z2=predict(loess(y1_del~x1_del,span = 0.5),x1)
          
          
          
          peaks_ppm=x1[pp]
          peaks_height=y1[pp]
          
          #    plot(x1,y1,type="l")
          #    points(x1,z1,type="l",col=2)
          #    points(x1[pp],y1[pp],cex=3)
          #    points(x1,z2,type="l",col=3)
          
        }
      }
    }
    
    return(list(baseline=list(x=x1,y=z2),ppm=mean(x1[pp]),peaks=list(ppm=peaks_ppm,height=peaks_height)))
    
    
  }











TSP_ref =
  function(x,y){
    
    
    fr1 <- function(z) {   
      fwhmopt <- z[1]
      etaopt <- z[2]
      a=y1-TSP(x1,height,ppm,fwhmopt,etaopt)
      sum(a*a)
    }
    fr2 <- function(z) {   
      heightopt <- z[1]
      ppmopt <- z[2]
      a=y1-TSP(x1,heightopt,ppmopt,fwhm,eta)
      sum(a*a)
    }
    
    range=c(-0.1,0.1)
    sel=x>range[1] & x<range[2]
    x1=x[sel]
    y1=y[sel]
    pp=which.max(y1)
    
    ppm=x1[pp]
    height=y1[pp]
    
    range=c(-0.03,0.03)+ppm
    sel=x>range[1] & x<range[2]
    x1=x[sel]
    y1=y[sel]
    pp=which.max(y1)
    
    ppm=x1[pp]
    height=y1[pp]
    
    
    
    
    
    o=optim(c(0.002,1), fr1,method = "L-BFGS-B",
            lower = c(0.005,0.5), 
            upper = c(0.03,1.5))    
    
    
    fwhm=abs(o$par[1])
    eta=o$par[2]
    #  o=optim(c(height,ppm), fr2)
    
    o=optim(c(height,ppm), fr2,method = "L-BFGS-B",
            lower = c(height,ppm)-c(height*0.005,0.0001), 
            upper = c(height,ppm)+c(height*0.05,0.0001))
    
    
    
    height=abs(o$par[1])
    ppm=o$par[2]
    
    area=sum(TSP(x1,height,ppm,fwhm,eta))
    RMSD=sqrt(o$value/length(x1))
    return(list(ppm=ppm,height=height,fwhm=fwhm,eta=eta,area=area,RMSD=RMSD))
    
  }


tsp_allignment =
function(spectra,range=c(-0.2,0.2)){
  TSP=NULL
  range=sort(range)
  for(k in 1:spectra$n){    
    u=spectra$Spectra[[k]]
    x=u[[2]]
    y=u[[1]]
    sel=x>range[1] & x<range[2]
    x1=x[sel]
    y1=y[sel]
    TSP[k]=x1[which.max(y1)]
  } 
  TSP
}
TSP = function(x1,height,ppm,fwhm,eta) {
  fitting=voight(height,x1,ppm,fwhm,eta)
  fitting=fitting+voight(height/35,x1,ppm-0.005497886,fwhm,eta)
  fitting=fitting+voight(height/35,x1,ppm+0.005497886,fwhm,eta)
  fitting
}

 allignment_spectra = function(Source,diff){
  
  Experiments <- Source
  rbnmr <- function (i)
  { 
    setwd (Experiments[i])
    Procs <- readLines("procs", n = -1)
    pma=pmatch("##$OFFSET=",Procs)
    OFFSET <- as.numeric(gsub(".*[=]", "",Procs[pma]))
    Procs[pma]=paste("##$OFFSET=",round(OFFSET-diff[i],digits=6))
    writeLines(Procs,con="procs")
    i=i+1
    
  }
  Spec <- lapply (1:length(Experiments),rbnmr)
  
}
pick.peaks = function (x, span) 
{
  span.width <- span * 2 + 1
  loc.max <- span.width + 1 - apply(embed(x, span.width), 1, 
                                    which.max)
  loc.max[loc.max == 1 | loc.max == span.width] <- NA
  pks <- loc.max + 0:(length(loc.max) - 1)
  unique(pks[!is.na(pks)])
}
voight = function(A, x, x0, fwhm, eta){
  G= exp( -(x-x0)^2 * 4 * log(2)  / (fwhm^2))
  L= (0.5 * fwhm)^2 / ((x-x0)^2 + (0.5 * fwhm)^2)
  A*((1-eta) * G + eta * L)
}


