library("bio3d")
dcdfile <- "1074trj-cas100-500ns.dcd"
pdbfile <- "499CA.pdb"
dcd <- read.dcd(dcdfile)
pdb <- read.pdb(pdbfile)
ca.inds <- atom.select(pdb, elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
#md_msf and md_cof
rf <- rmsf(xyz[,ca.inds$xyz])
msf<-rf^2
msf_n<-msf*2000/2001
cof<-matrix(0,909,909)
for (i in 1:303){
  for (j in 1:303){
    cof[i,j]=rmsf_n[i]*rmsf_n[j]*cij[i,j]
  }
}
cij<-dccm(xyz[,ca.inds$xyz])