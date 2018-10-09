load("allpineGATK.Rdata")

library(related)

bi2di<-function(x) {
SNPvect<-array(NA, length(x)*2)
L<-length(SNPvect)
SNPvect[{1:L} %% 2==1]<-ifelse(is.na(x), 0, ifelse(x==0, 1, ifelse(x==1, 1, 2)))
SNPvect[{1:L} %% 2==0]<-ifelse(is.na(x), 0, ifelse(x==0, 1, ifelse(x==1, 2, 2)))
return(SNPvect)
}

allpine <- allpine[rowSums(!is.na(allpine))>200,]

# index for the different seedlings and parents
Par <-grep("^T|^BD|-AC|-Z|-X|-Y|^Root", rownames(allpine))
Ref <- grep("^T|^B|SAL|REI|PAS|-21-|-60-|^5......$|^5....$", rownames(allpine))
Klocke <- which(grepl("-03-|-07-|-22-|^Root", rownames(allpine)) & !grepl("^LJPy", rownames(allpine)))
Vast <- grep("-14-|-AC|-Z|-X|-Y", rownames(allpine))
Lill <- which(grepl("^1", rownames(allpine)) & nchar(rownames(allpine))<8)



off<-c(Klocke, Vast, Lill)
gendat<-t(apply(allpine, 1, bi2di))
gendata<-data.frame(rownames(gendat), gendat)
gendata[,1]<-as.character(gendata[,1])

L_parent<-split(sample(Par), ceiling(seq_along(Par)/24))

Freq<-readgenotypedata(gendata[Ref,])$freqs

dat <- readgenotypedata(gendata)


library(related)
# ritland 1996 best option for GBS according to Attard et al. 2017 Mol Ecol Res.
# method=2 gives 95% CI

for(i in 1:length(off)){
        id<-off[i]
                for(j in 1:length(L_parent)){
                        dat$gdata<-gendata[c(id,L_parent[[j]]),]
                        dat$ninds<-nrow(dat$gdata)
                        output<-coancestry(dat$gdata, error.rates =0.05, allele.freqs = Freq, ritland=2)[c("relatedness", "relatedness.ci95", "inbreeding", "inbreeding.ci95")]
                        nam<-gendata[off[i],1]
                        id_off<-which(output$relatedness$ind1.id==nam | output$relatedness$ind2.id==nam)
                        output$relatedness<-output$relatedness[id_off,c(2,3,9,10)]
                        output$relatedness.ci95<-output$relatedness.ci95[id_off,c(2,3,13,14)]
                        if(j>1) {
                                output$inbreeding<-output$inbreeding[-1,]
                                output$inbreeding.ci95<-output$inbreeding.ci95[-1,]
                        }
                        if(j==1)        Tut<-output else{
                                Tut <- mapply(rbind, Tut, output)}
                }
                Ord<-order(Tut$relatedness$ritland, decreasing=TRUE)
                Tut2<-list()
                Tut2$relatedness<-Tut$relatedness[Ord,]
                Tut2$relatedness.ci95<-Tut$relatedness.ci95[Ord,]
                Tut2$inbreeding<-Tut$inbreeding[c(1,Ord+1),]
                Tut2$inbreeding.ci95<-Tut$inbreeding.ci95[c(1,Ord+1),]
                Tut2$inbreeding.ci95$ind.id<-Tut2$inbreeding$ind.id
                idx<-which(Tut2$relatedness$ritland>0)
                Tut3<-append(lapply(Tut2[1:2], function(x) x[idx,]),
                        lapply(Tut2[3:4], function(x) x[c(1, idx+1),]))
                if(i==1)        Out<-Tut3 else{
                        Out <- mapply(rbind, Out, Tut3)}
        if(i %% 20 == 0) {
                print(i); Sys.time()
                save.image(paste("Offtemp_", i, ".Rdata", sep=""))
        }
}

save.image("Related_output.Rdata")
