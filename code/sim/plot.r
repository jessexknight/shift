library('ggplot2')
library('reshape2')

plot.save = function(...,w=6,h=4,ext='.pdf'){
  fname = str(root.path('out','fig',...,create=TRUE),ext)
  status(3,'saving: ',fname)
  ggsave(fname,w=w,h=h)
}

plot.clean = function(g,...){
  g = g + theme_light() + theme(...,
    strip.background=element_rect(fill='gray85'),
    strip.text.x=element_text(color='black'),
    strip.text.y=element_text(color='black'))
}
