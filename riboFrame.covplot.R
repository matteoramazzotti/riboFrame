options <- commandArgs(trailingOnly = TRUE)
cat(options[1], ", ",options[2], ", ",options[3], ", ",options[4],"\n")
varst<-c(69,137,433,576,822,986,1117,1243,1435)
varen<-c(99,242,497,682,879,1043,1173,1294,1465)
labels=c("V1","V2","V3","V4","V5","V6","V7","V8","V9")
pdf(file=paste(options[1],".pdf", sep=""), width=10, height=5)
max<-0
for (i in 1:length(options)) {
	m<-max(read.table(file=options[i],sep="\t", header=T)$cnt)
	if (m > max) {
		max<-m
	}
}
for (i in 1:length(options)) {
	coverage<-read.table(file=options[i],sep="\t", header=T)
	if (i == 1) {
		plot(coverage$pos,coverage$cnt,type="l",xlab="Position", ylab="Coverage", xaxt="n", ylim=c(0,max))
	} 
	else {
		lines(coverage$pos, coverage$cnt,lty=i)
	}
}
abline(v=c(varen,varst),lty=5)
axis(3,at=(varen-varst)/2+varst,labels=labels)
axis(1,at=seq(1,1550,by=100)-1, labels=seq(1,1550,by=100)-1)
rect(varst,min(coverage$cnt),varen,max(coverage$cnt), col=rgb(1,0,0,0.3), border=NA)
dev.off()
