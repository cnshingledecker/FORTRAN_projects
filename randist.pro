PRO randist


vals = READ_CSV('granvals.csv')
;PRINT, vals.FIELD1

pdf = HISTOGRAM(vals.FIELD1,LOCATIONS=xbin)
N = N_ELEMENTS(vals.FIELD1)
PRINT, N
pdf = FLOAT(pdf)
N = FLOAT(N)
pdf = pdf/N

;FOR i=0,N-1 DO BEGIN
;  pdf[i] = pdf[i]/N
;  PRINT, i
;ENDFOR

PLOT, xbin, pdf, TITLE='Histogram', XTITLE='PDF', YTITLE='Probability',PSYM='2wc kkvkkkkvikkkk', LINESTYLE=' '


END