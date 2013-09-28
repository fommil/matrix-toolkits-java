par(cex = 1.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5, family="Palatino")

setPdfOut <- function(filename){
	if (interactive()) {
		quartz()
	} else {
		pdf(filename, width=11, height=8.5)
	}
}

logAxis = function(type, lims) {
	x1 <- floor(log10(lims))
	pow <- seq(x1[1], x1[2]+1)
	ticksat <- as.vector(sapply(pow, function(p) (1:10)*10^p))
	axis(type, 10^pow)
	axis(type, ticksat, labels=NA, tcl=-0.25, lwd=0, lwd.ticks=1)
}


data = read.csv("LinkedSparseMatrixPerf.csv",
				col.names=c("n", "m",
					 		"denseMem", "denseInit", "denseMult",
							 "sparseMem", "sparseInit", "sparseMult"))
data$denseInit = data$denseInit / 1000000000
data$sparseInit = data$sparseInit / 1000000000
data$denseMult = data$denseMult / 1000000000
data$sparseMult = data$sparseMult / 1000000000

ms = sort(unique(data$m))
ns = sort(unique(data$n))

xlim = c(min(ns), max(ns))
ylim = c(min(data$denseMult), max(data$denseMult))
setPdfOut("mult.pdf")
plot(c(), xlim=xlim, xlab="n", ylab="Time (seconds)", log="y", yaxt="n", ylim=ylim, type="n")
leg = matrix(nrow = 2 * length(ms), ncol = 3)
for (j in 1:length(ms)) {
	m = ms[j]
	print(m)
	set = data[data$m == m,]	
	avgs = matrix(nrow = length(ns), ncol = 3)
	for (i in 1:length(ns)) {
		avgs[i,1] = ns[i]
		avgs[i,2] = mean(set[set$n == ns[i],]$denseMult)
		avgs[i,3] = mean(set[set$n == ns[i],]$sparseMult)
	}
	col = heat.colors(length(ms))[j]
	lines(avgs[,1], avgs[,2], lty = 1, col = col)
	leg[2 * j - 1, 1] = paste("dense", m)
	leg[2 * j - 1, 2] = col
	leg[2 * j - 1, 3] = 1
	lines(avgs[,1], avgs[,3], lty = 2, col = col)
	leg[2 * j, 1] = paste("sparse", m)
	leg[2 * j, 2] = col
	leg[2 * j, 3] = 2
}
logAxis(2, ylim)

#legend("bottomright", legend=leg[,1], lty=leg[,3], col=leg[,2], bty="n")

ylim = c(min(data$denseInit), max(data$sparseInit))
setPdfOut("init.pdf")
plot(c(), xlim=xlim, xlab="n", ylab="Time (seconds)", log="y", yaxt="n", ylim=ylim, type="n")
leg = matrix(nrow = 2 * length(ms), ncol = 3)
for (j in 1:length(ms)) {
	m = ms[j]
	print(m)
	set = data[data$m == m,]	
	avgs = matrix(nrow = length(ns), ncol = 3)
	for (i in 1:length(ns)) {
		avgs[i,1] = ns[i]
		avgs[i,2] = mean(set[set$n == ns[i],]$denseInit)
		avgs[i,3] = mean(set[set$n == ns[i],]$sparseInit)
	}
	col = heat.colors(length(ms))[j]
	lines(avgs[,1], avgs[,2], lty = 1, col = col)
	leg[2 * j - 1, 1] = paste("dense", m)
	leg[2 * j - 1, 2] = col
	leg[2 * j - 1, 3] = 1
	lines(avgs[,1], avgs[,3], lty = 2, col = col)
	leg[2 * j, 1] = paste("sparse", m)
	leg[2 * j, 2] = col
	leg[2 * j, 3] = 2
}
logAxis(2, ylim)

ylim = c(min(data$sparseMem), max(data$denseMem))
setPdfOut("mem.pdf")
plot(c(), xlim=xlim, xlab="n", ylab="Memory (bytes)", log="y", yaxt="n", ylim=ylim, type="n")
leg = matrix(nrow = 2 * length(ms), ncol = 3)
for (j in 1:length(ms)) {
	m = ms[j]
	print(m)
	set = data[data$m == m,]	
	avgs = matrix(nrow = length(ns), ncol = 3)
	for (i in 1:length(ns)) {
		avgs[i,1] = ns[i]
		avgs[i,2] = mean(set[set$n == ns[i],]$denseMem)
		avgs[i,3] = mean(set[set$n == ns[i],]$sparseMem)
	}
	col = heat.colors(length(ms))[j]
	lines(avgs[,1], avgs[,2], lty = 1, col = col)
	leg[2 * j - 1, 1] = paste("dense", m)
	leg[2 * j - 1, 2] = col
	leg[2 * j - 1, 3] = 1
	lines(avgs[,1], avgs[,3], lty = 2, col = col)
	leg[2 * j, 1] = paste("sparse", m)
	leg[2 * j, 2] = col
	leg[2 * j, 3] = 2
}
logAxis(2, ylim)
