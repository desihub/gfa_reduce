pro dl_compare_daily_summary_files_relative_transparency

  nersc = mrdfits('/Users/dan.li/ScienceData/desi/survey/GFA/offline_matched_coadd_ccds_SV3-thru_20240114_nersc.fits', 2)
  kpno = mrdfits('/Users/dan.li/ScienceData/desi/survey/GFA/offline_matched_coadd_ccds_SV3-thru_20240114_kpno.fits', 2)

  nights = [nersc.night, kpno.night]
  nights = nights[uniq(nights, sort(nights))]
  
  nightMin = 20230830l
  nightMax = 99999999l

  _index = where ( nersc.night ge nightMin and nersc.night le nightMax )
  nersc = nersc[_index]
  kpno = kpno[_index]

  nightMax = max(nersc.night)

  recordNumber = indgen(n_elements(nersc))

  nerscData = nersc.transparency
  kpnoData = kpno.transparency

  dataToPlot = (kpnoData-nerscData)/nerscData

  dataMax = max([dataToPlot, dataToPlot], min=dataMin)
  ;yrange = [dataMin,dataMax]
  yrange = [-0.04, 0.04]
  
  position = [0.16, 0.15, 0.95, 0.95] 
  charSize = 0.75
  lineThick = 2.0

  histYRange = [0, 1000]

   ; create EPS  
  !p.font = 0
  set_plot, 'ps'
  dl_loadct

  epsFileName = '/Users/dan.li/ScienceData/desi/survey/GFA/relative_transimission.eps'
  
  device, /enca, filename=epsFileName, decomposed=0, $
          xsize=4, ysize=2.5, /inches, /color, Bits_per_Pixel=8

  plot, [0], [0], /nodata, $
        xstyle=1, xrange=[0,max(recordNumber)], xtitle='Record #', xthick=lineThick, $
        yrange=yrange, /ynozero, ytitle='Transmission (KPNO-NERSC)/NERSC', ythick=lineThick, $
        position=position, charsize=charSize

  oplot, recordNumber, dataToPlot

  device, /close_file
  spawn, 'epstopdf '+epsFileName
  spawn, 'rm '+epsFileName

  epsFileName = '/Users/dan.li/ScienceData/desi/survey/GFA/transmission_diff_histogram.eps'
  
  device, /enca, filename=epsFileName, decomposed=0, $
          xsize=4, ysize=2.5, /inches, /color, Bits_per_Pixel=8

  min = -0.02
  max = 0.02
  binSize = 0.001
  histY = histogram(dataToPlot, min=min, max=max, binSize=binSize)
  histX = findgen(n_elements(histY))*binSize+min+0.5*binSize

  _ind = where(histX gt -0.001 and histX lt 0.001)
  print, total(histY[_ind])/total(histY)
  _ind = where(histX gt -0.002 and histX lt 0.002)
  print, total(histY[_ind])/total(histY)
  _ind = where(histX gt -0.003 and histX lt 0.003)
  print, total(histY[_ind])/total(histY)
  _ind = where(histX gt -0.004 and histX lt 0.004)
  print, total(histY[_ind])/total(histY)
  _ind = where(histX gt -0.005 and histX lt 0.005)
  print, total(histY[_ind])/total(histY)

  plot, [0], [0], /nodata, $
        xstyle=1, xrange=[min,max], xtitle='Relative difference in transmission', xthick=lineThick, $
        yrange=histYRange, ytitle='N', ythick=lineThick, $
        position=position, charsize=charSize

  oplot, histX, histY, psym=10, thick=lineThick

  device, /close_file
  spawn, 'epstopdf '+epsFileName
  spawn, 'rm '+epsFileName

  nerscData = nersc.FWHM_ASEC
  kpnoData = kpno.FWHM_ASEC

  dataToPlot = (kpnoData-nerscData)/nerscData
  dataMax = max([dataToPlot, dataToPlot], min=dataMin)
  ;yrange = [dataMin,dataMax]
  yrange = [-0.005,0.005]
  
  epsFileName = '/Users/dan.li/ScienceData/desi/survey/GFA/relative_fwhm.eps'
  
  device, /enca, filename=epsFileName, decomposed=0, $
          xsize=4, ysize=2.5, /inches, /color, Bits_per_Pixel=8

  plot, [0], [0], /nodata, $
        xstyle=1, xrange=[0,max(recordNumber)], xtitle='Record #', xthick=lineThick, $
        yrange=yrange, /ynozero, ytitle='FWHM (KPNO-NERSC)/NERSC', ythick=lineThick, $
        position=position, charsize=charSize

  oplot, recordNumber, dataToPlot

  device, /close_file
  spawn, 'epstopdf '+epsFileName
  spawn, 'rm '+epsFileName


  epsFileName = '/Users/dan.li/ScienceData/desi/survey/GFA/fwhm_diff_histogram.eps'
  
  device, /enca, filename=epsFileName, decomposed=0, $
          xsize=4, ysize=2.5, /inches, /color, Bits_per_Pixel=8

  min = -0.002
  max = 0.002
  binSize = 0.0001
  histY = histogram(dataToPlot, min=min, max=max, binSize=binSize)
  histX = findgen(n_elements(histY))*binSize+min+0.5*binSize

  plot, [0], [0], /nodata, $
        xstyle=1, xrange=[min,max], xtitle='Relative difference in FWHM', xthick=lineThick, $
        yrange=histYRange, ytitle='N', ythick=lineThick, $
        position=position, charsize=charSize

  oplot, histX, histY, psym=10, thick=lineThick

  device, /close_file
  spawn, 'epstopdf '+epsFileName
  spawn, 'rm '+epsFileName
  
  close_eps:
  ; close EPS
  loadct, 0, /silent
  device, decomposed=1
  device, /close_file
  set_plot, 'X'
  !p.font = -1

end
