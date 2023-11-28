pro st2,h,ss,k,negative=negative
;+
;PURPOSE:
;   S-transform ver2
;
;INPUT:
;   h:		time series, complex, n points
;
;OUTPUT:
;   ss:		S-transform spectrum, complex, n x n points
;
;OPTION:
;   k:		width of Gaussian window (default k=1)
;
;KEYWORD
;   negative:		calculate negative frequency
;
;ryuho@nict.go.jp, 2004/11/03
;ryuho@nict.go.jp, 2005/05/29 generalized k
;ryuho@nict.go.jp, 2006/02/03 ver2, negative frequency
;-

if n_elements(k) eq 0 then k=1.0

n=n_elements(h)
if n mod 2 ne 0 then h=h(0:n-2)
n=n_elements(h)

gg=fltarr(n)			;voice Gaussian
ss=complexarr(n,n/2)	;S-transform spectrum
hh=fft(h,-1)			;forward Fourier transform

IF KEYWORD_SET(negative) THEN BEGIN

  hh=reverse(hh)  ;negative frequency are first
  hh=shift(hh,1)  ;shift DC to index 0

ENDIF

for f=1.,n/2. do begin

  gg=exp(-2.*(!dpi*(findgen(n)-(n/2))*k/f)^2.)
  gg=shift(gg,-n/2)
  hh=shift(hh,-1)
  ss(*,f-1)=fft(hh*gg,1)	;inverse Fourier transform

endfor

IF KEYWORD_SET(negative) THEN ss=reverse(ss,1)

end
