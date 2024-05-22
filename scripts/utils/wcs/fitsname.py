

fitsname="dirty_image_20230601T104250099_real.fits"

dtm_utc=fitsname[12:30]
year=int(dtm_utc[0:4])
month=int(dtm_utc[4:6])
day=int(dtm_utc[6:8])
hour=int(dtm_utc[9:11])
min=int(dtm_utc[11:13])
sec=int(dtm_utc[13:15])
msec=int(dtm_utc[15:18])

# 20230601T104250099 -> 2023-08-22T17:18:41.1
# fits[0].header['DATE-OBS'] = '2023-06-01T10:42:50.1'
dtm_out=('%04d-%02d-%02dT%02d:%02d:%02d.%03d' % (year,month,day,hour,min,sec,msec))
print("%s -> %s" % (dtm_utc,dtm_out))

