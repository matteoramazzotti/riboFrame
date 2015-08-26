ranks<-c("domain","phylum","order","class","family","genus")
for (i in 1:length(ranks)) {
file<-paste("xxxxx/pyro_v1v3_count.",ranks[i],".cnt",sep="")
if (file.exists(file)==T) {
pyro1<-read.table(file,sep="\t",header=T)
}
file<-paste("xxxxx/pyro_v3v5_count.",ranks[i],".cnt",sep="")
if (file.exists(file)==T) {
pyro2<-read.table(file,sep="\t",header=T)
}
file<-paste("xxxxx/ilmn_v1v3_count.",ranks[i],".cnt",sep="")
if (file.exists(file)==T) {
ilmn1<-read.table(file,sep="\t",header=T)
}
file<-paste("xxxxx/ilmn_v3v5_count.",ranks[i],".cnt",sep="")
if (file.exists(file)==T) {
ilmn2<-read.table(file,sep="\t",header=T)
}
file<-paste("xxxxx/ilmn_full_count.",ranks[i],".cnt",sep="")
if (file.exists(file)==T) {
ilmn3<-read.table(file,sep="\t",header=T)
}
file<-paste("xxxxx/ilmn_all_count.",ranks[i],".cnt",sep="")
if (file.exists(file)==T) {
ilmn4<-read.table(file,sep="\t",header=T)
}

if(exists("pyro1") && exists("ilmn1")) {
list<-sort(unique(union(pyro1$Name,ilmn1$Name)))
x<-ifelse(match(list,pyro1$Name),pyro1$Perc[match(list,pyro1$Name)],0)
y<-ifelse(match(list,ilmn1$Name),ilmn1$Perc[match(list,ilmn1$Name)],0)
x<-ifelse(is.na(x),0,x)
y<-ifelse(is.na(y),0,y)
cat(ranks[i],"pyro_v1v3 vs.","ilmn_v1v3 R: ",cor(x,y),"D:",sqrt(sum((x-y)^2)),"slope:",coef(lm(x~y))[2],"\n")
print(data.frame(list,x,y))
}

if(exists("pyro1") && exists("ilmn2")) {
list<-sort(unique(union(pyro1$Name,ilmn2$Name)))
x<-ifelse(match(list,pyro1$Name),pyro1$Perc[match(list,pyro1$Name)],0)
y<-ifelse(match(list,ilmn2$Name),ilmn2$Perc[match(list,ilmn2$Name)],0)
x<-ifelse(is.na(x),0,x)
y<-ifelse(is.na(y),0,y)
cat(ranks[i],"pyro_v1v3 vs.","ilmnv_3v5 R: ",cor(x,y),"D:",sqrt(sum((x-y)^2)),"slope:",coef(lm(x~y))[2],"\n")
print(data.frame(list,x,y))
}
if(exists("pyro1") && exists("ilmn3")) {
list<-sort(unique(union(pyro1$Name,ilmn3$Name)))
x<-ifelse(match(list,pyro1$Name),pyro1$Perc[match(list,pyro1$Name)],0)
y<-ifelse(match(list,ilmn3$Name),ilmn3$Perc[match(list,ilmn3$Name)],0)
x<-ifelse(is.na(x),0,x)
y<-ifelse(is.na(y),0,y)
cat(ranks[i],"pyro_v1v3 vs.","ilmn_full R: ",cor(x,y),"D:",sqrt(sum((x-y)^2)),"slope:",coef(lm(x~y))[2],"\n")
print(data.frame(list,x,y))
}
if(exists("pyro1") && exists("ilmn4")) {
list<-sort(unique(union(pyro1$Name,ilmn4$Name)))
x<-ifelse(match(list,pyro1$Name),pyro1$Perc[match(list,pyro1$Name)],0)
y<-ifelse(match(list,ilmn4$Name),ilmn4$Perc[match(list,ilmn4$Name)],0)
x<-ifelse(is.na(x),0,x)
y<-ifelse(is.na(y),0,y)
cat(ranks[i],"pyro_v1v3 vs.","ilmn_all R: ",cor(x,y),"D:",sqrt(sum((x-y)^2)),"slope:",coef(lm(x~y))[2],"\n")
print(data.frame(list,x,y))
}

if(exists("pyro2") && exists("ilmn1")) {
list<-sort(unique(union(pyro2$Name,ilmn1$Name)))
x<-ifelse(match(list,pyro2$Name),pyro2$Perc[match(list,pyro2$Name)],0)
y<-ifelse(match(list,ilmn1$Name),ilmn1$Perc[match(list,ilmn1$Name)],0)
x<-ifelse(is.na(x),0,x)
y<-ifelse(is.na(y),0,y)
cat(ranks[i],"pyro_v3v5 vs.","ilmn_v1v3 R: ",cor(x,y),"D:",sqrt(sum((x-y)^2)),"slope:",coef(lm(x~y))[2],"\n")
print(data.frame(list,x,y))
}
if(exists("pyro2") && exists("ilmn2")) {
list<-sort(unique(union(pyro2$Name,ilmn2$Name)))
x<-ifelse(match(list,pyro2$Name),pyro2$Perc[match(list,pyro2$Name)],0)
y<-ifelse(match(list,ilmn2$Name),ilmn2$Perc[match(list,ilmn2$Name)],0)
x<-ifelse(is.na(x),0,x)
y<-ifelse(is.na(y),0,y)
cat(ranks[i],"pyro_v3v5 vs.","ilmn_v3v5 R: ",cor(x,y),"D:",sqrt(sum((x-y)^2)),"slope:",coef(lm(x~y))[2],"\n")
print(data.frame(list,x,y))
}
if(exists("pyro2") && exists("ilmn3")) {
list<-sort(unique(union(pyro2$Name,ilmn3$Name)))
x<-ifelse(match(list,pyro2$Name),pyro2$Perc[match(list,pyro2$Name)],0)
y<-ifelse(match(list,ilmn3$Name),ilmn3$Perc[match(list,ilmn3$Name)],0)
x<-ifelse(is.na(x),0,x)
y<-ifelse(is.na(y),0,y)
cat(ranks[i],"pyro_v3v5 vs.","ilmn_full R: ",cor(x,y),"D:",sqrt(sum((x-y)^2)),"slope:",coef(lm(x~y))[2],"\n")
print(data.frame(list,x,y))
}
if(exists("pyro2") && exists("ilmn4")) {
list<-sort(unique(union(pyro2$Name,ilmn4$Name)))
x<-ifelse(match(list,pyro2$Name),pyro2$Perc[match(list,pyro2$Name)],0)
y<-ifelse(match(list,ilmn4$Name),ilmn4$Perc[match(list,ilmn4$Name)],0)
x<-ifelse(is.na(x),0,x)
y<-ifelse(is.na(y),0,y)
cat(ranks[i],"pyrov_3v5 vs.","ilmn_all R: ",cor(x,y),"D:",sqrt(sum((x-y)^2)),"slope:",coef(lm(x~y))[2],"\n")
print(data.frame(list,x,y))
}
}
