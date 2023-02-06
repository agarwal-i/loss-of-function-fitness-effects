### input files from https://github.com/zfuller5280/MutationSelection

library(data.table)

# par
par=fread("par_posteriors.csv", header=T)
par=par[,lapply(.SD, as.numeric)]
par=log10(par)

# X chromosome
post_xs=fread("s_posteriors.xchr.csv", header=T)
post_xh=fread("h_posteriors.xchr.csv", header=T)
names(post_xs)==names(post_xh)
post_xhs = post_xh*post_xs
post_x_sexavg = log10((0.5*post_xhs)+(0.5*post_xs))
post_xhs = log10(post_xhs)

# autosomes
post=fread("updated_posteriors.12_16.csv.gz", header=T)
post=post[,lapply(.SD, as.numeric)]
post=log10(post)
#post=post[ , which(sapply(post, function(x) all(is.na(x)))) := NULL]
post=post[complete.cases(post),]

###################################################################################
# calculate densities
###################################################################################

minl=log10(0.1*1e-6)
maxl=log10(1)

xax=vector('list',1)
xax[[1]]=density(as.numeric(unlist(post[,1])), kernel = c("gaussian"), n = 1024, from = minl, to = maxl, na.rm=T)$x
xax <- setNames(xax,"x_axis")
saveRDS(xax,"xax.Rda")

yax=vector('list',ncol(post))
for (i in 1:ncol(post)){
yax[[i]]=density(as.numeric(unlist(post[, ..i])), kernel = c("gaussian"), n = 1024, from = minl, to = maxl, na.rm=T)$y
}
yax <- setNames(yax,names(post))
saveRDS(yax,"yax_aut.Rda")

yax_xhs=vector('list',ncol(post_xhs))
for (i in 1:ncol(post_xhs)){
yax_xhs[[i]]=density(as.numeric(unlist(post_xhs[, ..i])), kernel = c("gaussian"), n = 1024, from = minl, to = maxl, na.rm=T)$y
}
yax_xhs <- setNames(yax_xhs,names(post_xhs))
saveRDS(yax_xhs,"yax_xhs.Rda")

yax_xsa=vector('list',ncol(post_x_sexavg))
for (i in 1:ncol(post_x_sexavg)){
yax_xsa[[i]]=density(as.numeric(unlist(post_x_sexavg[, ..i])), kernel = c("gaussian"), n = 1024, from = minl, to = maxl, na.rm=T)$y
}
yax_xsa <- setNames(yax_xsa,names(post_x_sexavg))
saveRDS(yax_xsa,"yax_xsa.Rda")

yax_par=vector('list',ncol(par))
for (i in 1:ncol(par)){
yax_par[[i]]=density(as.numeric(unlist(par[, ..i])), kernel = c("gaussian"), n = 1024, from = minl, to = maxl, na.rm=T)$y
}
yax_par <- setNames(yax_par,names(par))
saveRDS(yax_par,"yax_par.Rda")

