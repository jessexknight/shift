library('ggplot2')
library('reshape2')

plot.save = function(...,w=6,h=4,ext='.pdf'){
  fname = paste0(root.path('out','fig',...,create=TRUE),ext)
  cat(paste('saving:',fname,'\n'))
  ggsave(fname,w=w,h=h)
}

plot.clean = function(g,...){
  g = g + theme_light() + theme(...,
    strip.background=element_rect(fill='gray85'),
    strip.text.x=element_text(color='black'),
    strip.text.y=element_text(color='black'))
}

plot.rr = function(rr.,rr){
  melt.rr = function(rr){
    melt(rbind.lapply(names(rr.),function(e){
      cbind(e=e,as.data.frame(rr[[e]])) }),id=c('e','t'),measure='rr') }
  X. = melt.rr(rr.)
  X  = melt.rr(rr)
  g = ggplot(X,aes(x=t,y=value,color=e)) +
    geom_hline(yintercept=1,color='gray',linetype='33') +
    geom_line(alpha=.5) +
    geom_point(shape=16,size=1) +
    geom_point(aes(shape=e),size=2,data=X.) +
    scale_shape_manual(values=c(0,1,2,5,6)) +
    labs(x='Days',y='Relative Rate',color='Event',shape='Event')
  g = plot.clean(g)
}
