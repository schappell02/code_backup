;-

PRO ppf_star, star_fits, young_fits, old_fits, dir, guess

s_flux = readfits(dir+star_fits)
tmp_yflux = readfits(dir+young_fits)
tmp_oflux = readfits(dir+old_fits)

shead = headfits(dir+star_fits)
yhead = headfits(dir+young_fits)
ohead = headfits(dir+old_fits)

num = sxpar(shead,'NAXIS1')
wmin = sxpar(shead,'CRVAL1')
wsize = sxpar(shead,'CDELT1')
wave = findgen(num)*wsize + wmin

onum = sxpar(ohead,'NAXIS1')
owmin = sxpar(ohead,'CRVAL1')
owsize = sxpar(ohead,'CDELT1')
owave = findgen(onum)*owsize + owmin

ynum = sxpar(yhead,'NAXIS1')
ywmin = sxpar(yhead,'CRVAL1')
ywsize = sxpar(yhead,'CDELT1')
ywave = findgen(ynum)*ywsize + ywmin

y_flux = interpol(tmp_yflux,ywave,wave,/quadratic)
o_flux = interpol(tmp_oflux,owave,wave,/quadratic)

template = dblarr(num,2)
template[*,0]= y_flux/median(y_flux)
template[*,1]= o_flux/median(o_flux)

ppxf,template,s_flux,s_flux*0+1.0,34.6,guess,moments=2

save, filename='check.sav',/all

end
