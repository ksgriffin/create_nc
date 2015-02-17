;load "/lb11/kgriffin/maps/lib/populate_gfs_level.ncl"
;load "/lb11/kgriffin/maps/lib/update_gfs_level.ncl"
;load "/lb11/kgriffin/maps/lib/update_gfs_level_fhr.ncl"
load "/home/ksgriffin2/nc_create/add_update_gfs_fhr.ncl"
;load "/home/ksgriffin2/nc_create/update_gfs_fhr.ncl"
load "/home/ksgriffin2/nc_create/dynamaps.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

begin

  loc = "/nialldata/ksgriffin2/realtime/data/new_gfs_nc/"
  logdir = "/nialldata/ksgriffin2/dynamaps/flogs/"
  timeUnits = "hours since 1800-01-01 00:00:00.0"

  system("rm -f /home/ksgriffin2/derived_log.file")

  levels = (/10,20,50,100,150,200,250,300,500,700,850,925,1000/)
  times = ispan(6,6,6)
;  times = (/0,6/)

  nlevels = toint(dimsizes(levels))
  ntimes = toint(dimsizes(times))

  rawdate = calc_date_gfscycle(timeUnits)
  strdate = cd_string(rawdate,"%y%N%D_%H")
  ;strdate = "15021312"
  ;system("rm -f "+logdir+"f???.txt")


  do i=0,ntimes-1
      print("time: "+times(i))
      logfilename = logdir+strdate+"_f"+sprinti("%0.3i",times(i))+".txt"
      update_gfs_fhr(loc,levels,True,times(i))
      system("echo NetCDF file produced. `date` > "+logfilename)
  end do

end
