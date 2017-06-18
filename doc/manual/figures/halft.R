pdf("halft.pdf", width=6, height=4)
opar <- par(lwd=1.5, mar=c(3,2,1,1)+0.1)
x <- seq(from=0, to=5.15, length=1001)
y1 <- dt(x, df=1)*2
y2 <- dt(x, df=2)*2
plot(x,y2,type="l",ylim=c(0,0.8), yaxs="i", yaxt="n", xaxs="i", xaxt="n",
     xlab="sigma",ylab="density")
axis(side=1, at=0:5, label=expression("0","S", 2 %*% S, 3 %*% S, 4 %*% S,
                                      5 %*% S))
ynorm <- dnorm(x)*2
#lines(x,y2,lty=2)
lines(x,ynorm,lty=2)

legend("topright", lty=1:2, legend=expression(t[2], t[infinity]))
par(opar)
dev.off()
