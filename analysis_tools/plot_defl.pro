
fname = 'defl_grid'
rows = file_lines(fname+'.dat')
in   = fltarr(2,rows)
openr,1,fname+'.dat'
readf,1,in
close,1







y   = 0.89
top = 0.01
bot = 0.1
;h = 5.6202
h   = 14.5502
w   = 14.5502
;w = 10.1165






!p.font = 1
set_plot,'PS'
device,filename=fname+'.ps',xsize=w,ysize=h,font_size=14,/tt_font,set_font='Times',/color,decomposed=1

range = 0.35

;========= Plot nodata ========
plot,[-range,range],[-range,range],/isotropic,/nodata,background=255,$
  xstyle=1,xtickinterval=0.2,xminor=5,xticklen=0.02,$
  ystyle=1,ytickinterval=0.2,yminor=5,yticklen=0.02


;========= Plot data ========
oplot,in[0,*],in[1,*],psym=3,symsize=0.05,color='3000cc'x

plots,[0,0],[-range,range]
plots,[-range,range],[0,0]

device,/close
set_plot,'X'


spawn,'convert -density 250x250 '+fname+'.ps '+fname+'.png'
spawn,'rm '+fname+'.ps'


end
