load "/home/ksgriffin2/scripts/ncl/lib/load.ncl"

;**********
; Passed variables:
; filename --- the name of the netCDF file with the FHR of interest
; date ------- date of the model cycle
; cycle ------ hour of the model cycle
; fhr -------- forecast hour of the file (not everything is done at all hours)

begin

  strfhr = sprinti("%0.3i",fhr)
  strcycle = sprinti("%0.2i",cycle)
  outloc = "/nialldata/ksgriffin2/realtime/data/new_gfs_nc/"
  outfilename = outloc+"derivedgfs"+date+"_"+strcycle+"_F"+strfhr+".nc"
  system("rm -f "+outfilename)
  moist_levels = (/1000,925,850,700,500,250/)
  temp_levels = (/1000,925,850,700,500,300,250,200,150,100,50,20,10/)
  jet_levels = (/300,250,200/)
  vort_levels = (/925,850,500,250/)
  geo_levels = (/700,500,300/)
  anom_levels = (/1000,850,500,250/)

  infile = addfile(filename,"r")
  outfile = addfile(outfilename,"c")
  nlats = dimsizes(infile->lat)
  nlons = dimsizes(infile->lon)
  ntimes = dimsizes(infile->time)

  ;setfileoption(outfile,"DefineMode",True)
;	dimNames = (/"time", "lat", "lon", "moist_levels","temp_levels","jet_levels","vort_levels","geo_levels"/)  
;	dimSizes = (/ -1   ,  nlats,  nlons, dimsizes(moist_levels),dimsizes(temp_levels),dimsizes(jet_levels),dimsizes(vort_levels),dimsizes(geo_levels)/) 
;	dimUnlim = (/ True , False, False, False,False,False,False,False/)   
;	filedimdef(outfile,dimNames,dimSizes,dimUnlim)
;
 ; filevardef(outfile, "time" ,typeof(infile->time),getvardims(infile->time))
;  filevardef(outfile, "lat" ,typeof(infile->lat),getvardims(infile->lat))
;  filevardef(outfile, "lon" ,typeof(infile->lon),getvardims(infile->lon))
;  filevardef(outfile, "moist_levels" ,typeof(moist_levels),getvardims(moist_levels))
;  filevardef(outfile, "temp_levels" ,typeof(temp_levels),getvardims(temp_levels))
;  filevardef(outfile, "jet_levels" ,typeof(jet_levels),getvardims(jet_levels))
;  filevardef(outfile, "vort_levels" ,typeof(vort_levels),getvardims(vort_levels))
;  filevardef(outfile, "geo_levels" ,typeof(geo_levels),getvardims(geo_levels))
;
;  setfileoption(outfile,"DefineMode",False)


  ; For potential temperature
  theta = new((/ntimes,dimsizes(temp_levels),nlats,nlons/),float)
  theta!1 = "temp_levels"
  theta&temp_levels = temp_levels
  do i=0,dimsizes(temp_levels)-1
    thistemp := infile->t(:,{temp_levels(i)},:,:)
    theta(:,i,:,:) = pot_temp(thistemp,temp_levels(i)*100)
  end do  
  delete(thistemp)

  outfile->theta = theta

  ; For equivalent potential temperature
  thetae = new((/ntimes,dimsizes(moist_levels),nlats,nlons/),float)
  thetae!1 = "moist_levels"
  thetae&moist_levels = moist_levels
  do i=0,dimsizes(moist_levels)-1
    thistemp := infile->t(:,{moist_levels(i)},:,:)
    thisrh := infile->rh(:,{moist_levels(i)},:,:)
    thisq := mixhum_ptrh(conform(thistemp,moist_levels(i),-1),thistemp,thisrh,2)
    thetae(:,i,:,:) = equiv_pot_temp(thistemp,thisq,moist_levels(i))
  end do
  delete(thistemp)
  delete(thisrh)
  delete(thisq)

  outfile->thetae = thetae

  print("Now jet levels")
  ; For PV
  pv = new((/ntimes,dimsizes(jet_levels),nlats,nlons/),float)
  u_all = infile->u(:,{jet_levels},:,:)
  v_all = infile->v(:,{jet_levels},:,:)
  t_all = infile->t(:,{jet_levels},:,:)
  jet_levels@units = "Pa"
  print("pv now") 
  pv_all = PotVortIsobaric(jet_levels*100,u_all,v_all,t_all,t_all&lat,1,0)
  pv_all!1 = "jet_levels"
  jet_levels@units = "hPa"
  ; For isotachs
  print("isotachs now")
  mag_wind = (u_all^2+v_all^2)^0.5
  copy_VarCoords(u_all,mag_wind)
  mag_wind!1 = "jet_levels"
  print("divg/irrwind now")
  divg = uv2dvF_Wrap(u_all,v_all)
  divg!1 = "jet_levels"
  irr2 = dv2uvF_Wrap(divg)
  irr2!2 = "jet_levels"
  u_irr = irr2(0,:,:,:,:)
  v_irr = irr2(1,:,:,:,:)
  delete(u_all)
  delete(v_all)
  delete(t_all)
  delete(irr2)

  outfile->pv = pv_all
  outfile->mag_wind = mag_wind
  outfile->divg = divg
  outfile->u_irr = u_irr
  outfile->v_irr = v_irr

  print("Now rvort")
  ; For relative vorticity
  u_all = infile->u(:,{vort_levels},:,:)
  v_all = infile->v(:,{vort_levels},:,:)
  rvort = uv2vrF_Wrap(u_all,v_all)
  rvort!1 = "vort_levels"
  llvort = dim_avg_n_Wrap(rvort(:,{925:850},:,:),1)
  delete(u_all)
  delete(v_all)
  outfile->rvort = rvort
  outfile->llvort = llvort

  print("Now sfc dewp")
  ; For sfc dewp  
  tmpsfc = infile->t_sfc
  rhsfc = infile->rh_sfc
  td_sfc = dewtemp_trh(tmpsfc,rhsfc)
  copy_VarCoords(tmpsfc,td_sfc)
  delete(tmpsfc)
  delete(rhsfc)

  outfile->td_sfc = td_sfc

  print("Now geo levels")
  ; For geostrophic wind
  g_all = infile->g(:,{geo_levels},:,:)
  ;print(num(ismissing(g_all))+"")
  wind = z2geouv(g_all,g_all&lat,g_all&lon,1)
  u_geo = wind(0,:,:,:,:)
  v_geo = wind(1,:,:,:,:)
  ;print(num(ismissing(u_geo))+" "+num(ismissing(v_geo)))
  ;geo_vort = uv2vrF_Wrap(u_geo,v_geo)
  u_geo!1 = "geo_levels"
  u_geo!0 = "time"
  u_geo!2 = "lat"
  u_geo!3 = "lon"
  u_geo&time = g_all&time
  u_geo&lat = g_all&lat
  u_geo&lon = g_all&lon
  copy_VarCoords(u_geo,v_geo)
  ;geo_vort!1 = "geo_levels"
  delete(g_all)
  delete(wind)

  outfile->u_geo = u_geo
  outfile->v_geo = v_geo
  ;outfile->geo_vort = geo_vort


  ;For combining PV hemispheres onto one level
  print("Now combining PV hemispheres")
  dt_pres_comb = infile->dt_pres(:,{2},:,:)
  sh_pv_pres = infile->dt_pres(:,{-2},{-90:-0.5},:)
  dt_pres_comb(:,{-90:-0.5},:) = sh_pv_pres
  dt_temp_comb = infile->dt_temp(:,{2},:,:)
  sh_pv_temp = infile->dt_temp(:,{-2},{-90:-0.5},:)
  dt_temp_comb(:,{-90:-0.5},:) = sh_pv_temp
  dt_theta_comb = pot_temp(dt_temp_comb,dt_pres_comb)
  dt_pres_comb = dt_pres_comb/100.
  delete(sh_pv_pres)
  delete(sh_pv_temp)
  u_dt_comb = infile->u_dt(:,{2},:,:)
  sh_u_dt = infile->u_dt(:,{-2},{-90:-0.5},:)
  u_dt_comb(:,{-90:-0.5},:) = sh_u_dt
  v_dt_comb = infile->v_dt(:,{2},:,:)
  sh_v_dt = infile->v_dt(:,{-2},{-90:-0.5},:)
  v_dt_comb(:,{-90:-0.5},:) = sh_v_dt
  delete(sh_u_dt)
  delete(sh_v_dt)

  outfile->dt_pres_comb = dt_pres_comb
  outfile->dt_temp_comb = dt_temp_comb
  outfile->dt_theta_comb = dt_theta_comb
  outfile->u_dt_comb = u_dt_comb
  outfile->v_dt_comb = v_dt_comb

  ; For vertical wind shear
  print("Now vertical wind shear")
  shear_levels = (/850,700,500,250/) 
  u_winds = infile->u(:,{shear_levels},:,:)
  v_winds = infile->v(:,{shear_levels},:,:)
  u_deepshear = u_winds(:,{250},:,:)-u_winds(:,{850},:,:)
  v_deepshear = v_winds(:,{250},:,:)-v_winds(:,{850},:,:)
  mag_deepshear = (u_deepshear^2+v_deepshear^2)^0.5
  u_shortshear = u_winds(:,{500},:,:)-u_winds(:,{850},:,:)
  v_shortshear = v_winds(:,{500},:,:)-v_winds(:,{850},:,:)
  mag_shortshear = (u_shortshear^2+v_shortshear^2)^0.5
  u_deepshear!0 = "time"
  u_deepshear!1 = "lat"
  u_deepshear!2 = "lon"
  u_deepshear&time = u_winds&time
  u_deepshear&lat = u_winds&lat
  u_deepshear&lon = u_winds&lon
  copy_VarCoords(u_deepshear,v_deepshear)
  copy_VarCoords(u_deepshear,u_shortshear)
  copy_VarCoords(u_deepshear,v_shortshear)
  copy_VarCoords(u_deepshear,mag_deepshear)
  copy_VarCoords(u_deepshear,mag_shortshear)
  delete(u_winds)
  delete(v_winds)

  outfile->u_deepshear = u_deepshear
  outfile->v_deepshear = v_deepshear
  outfile->u_shortshear = u_shortshear
  outfile->v_shortshear = v_shortshear
  outfile->mag_deepshear = mag_deepshear
  outfile->mag_shortshear = mag_shortshear

  if(fhr.ne.0)
  print("Now p-type precip")
  ;print(getfilevarnames(infile))
  precip = infile->pcp_6h(:,:,:)
  ;;;;;not responsible to do with pcp_6h, except for GFS!
  rain_pcp6h = precip * (infile->cat_rain)
  frzrn_pcp6h = precip * (infile->cat_frzrn)
  ice_pcp6h = precip * (infile->cat_ice)
  snow_pcp6h = precip * (infile->cat_snow)
  copy_VarCoords(precip,rain_pcp6h)
  copy_VarCoords(precip,frzrn_pcp6h)
  copy_VarCoords(precip,ice_pcp6h)
  copy_VarCoords(precip,snow_pcp6h)
  outfile->rain_pcp6h = rain_pcp6h
  outfile->frzrn_pcp6h = frzrn_pcp6h
  outfile->ice_pcp6h = ice_pcp6h
  outfile->snow_pcp6h = snow_pcp6h
  end if

  ; For std anomalies
  print("Now std. anomalies")
  climo_loc = "/nialldata/ksgriffin2/cfsr_climo/"
  ;pwat = infile->pwat
  ;pwat_anom = pwat - pwat_climo
  ;pwat_std_anom = pwat_anom / pwat_stdev
  ;delete(pwat)
  ;outfile->pwat_std_anom = pwat_std_anom

  g_all = infile->g(:,{anom_levels},:,:)
  g_anoms = calc_anoms_rt(climo_loc,g_all,"g",anom_levels)
  g_std_anoms = calc_std_anoms_rt(climo_loc,g_anoms,"g",anom_levels)
  printVarSummary(g_std_anoms)
  print(min(g_std_anoms)+" "+max(g_std_anoms))
  outfile->g_anoms = g_anoms
  outfile->g_std_anoms = g_std_anoms
  delete(g_all)
  delete(g_anoms)
  delete(g_std_anoms)

  t_all = infile->t(:,{anom_levels},:,:)
  t_anoms = calc_anoms_rt(climo_loc,t_all,"t",anom_levels)
  t_std_anoms = calc_std_anoms_rt(climo_loc,t_anoms,"t",anom_levels)
  outfile->t_anoms = t_anoms
  outfile->t_std_anoms = t_std_anoms
  delete(t_all)
  delete(t_anoms)
  delete(t_std_anoms)

  logfilename = "/nialldata/ksgriffin2/dynamaps/flogs/"+date+"_"+strcycle+"_f"+strfhr+".txt"
  system("echo Derived GFS file produced. `date` >> "+logfilename)
  print("Done and success.")

end
