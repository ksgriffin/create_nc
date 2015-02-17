undef("add_cyclic")
function add_cyclic(data)
begin

  ntimes = dimsizes(data&time)
  nlats = dimsizes(data&lat)
  nlons = dimsizes(data&lon)
  lonmin = min(data&lon)
  lonmax = max(data&lon)
  dx = data&lon(1)-data&lon(0) 

  if(lonmin.eq.0.and.lonmax.eq.(360-dx)) then
    newdata = new((/ntimes,nlats,nlons+1/),typeof(data))
    newdata(:,:,0:nlons-1) = data
    newdata(:,:,nlons) = data(:,:,0)
  else if(lonmax.eq.360.and.lonmin.eq.(0+dx)) then
    newdata = new((/nlats,nlons+1/),typeof(data))
    newdata(:,:,1:nlons) = data
    newdata(:,:,0) = data(:,:,nlons)
  else if(lonmin.eq.-180.and.lonmax.eq.(180-dx))then
    newdata = new((/nlats,nlons+1/),typeof(data))
    newdata(:,:,0:nlons-1) = data
    newdata(:,:,nlons) = data(:,:,0)
  else if(lonmax.eq.180.and.lonmin.eq.(-180+dx))then
    newdata = new((/nlats,nlons+1/),typeof(data))
    newdata(:,:,1:nlons) = data
    newdata(:,:,0) = data(:,:,nlons)
  end if
  end if
  end if
  end if

  newdata&time = data&time
  newdata&lat = data&lat
  newdata&lon = fspan(lonmin,lonmax+dx,nlons+1)

  return(newdata)

end

undef("smooth_gempak")
function smooth_gempak(data)
begin

  smooth_data = smth9_Wrap(data, 0.50,  0.25, True)

  return(smooth_data)

end

undef("create_field")
function create_field(a:numeric)
begin

  data = create "data" scalarFieldClass defaultapp
        "sfDataArray" : a
        "sfMissingValueV" : a@_FillValue
        "sfXCStartV" : min(a&lon)
        "sfXCEndV": max(a&lon)
        "sfYCStartV" : min(a&lat)
        "sfYCEndV": max(a&lat)
  end create

  return(data)

end

undef("create_vector_field")
function create_vector_field(u:numeric,v:numeric)
begin

  data = create "data" vectorFieldClass defaultapp
    "vfUDataArray" : u
    "vfVDataArray" : v
    "vfMissingUValueV" : u@_FillValue
    "vfMissingVValueV" : v@_FillValue
    "vfXCStartV" : min(u&lon)
    "vfXCEndV" : max(u&lon)
    "vfYCStartV" : min(u&lat)
    "vfYCEndV" : max(u&lat)
  end create

  return(data)

end

undef("create_map")
function create_map(wks)
begin
  map = create "map" mapPlotClass wks
    "mpFillOn" : True
    "mpLandFillColor" : "gray80"
    "mpOceanFillColor" : "white"
  end create
  return(map)
end

undef("create_contour_plot")
function create_contour_plot(wks,data)
begin

  plot = create "plot" contourPlotClass wks
    "cnScalarFieldData" : data
  end create

  return(plot)

end

undef("create_vector_plot")
function create_vector_plot(wks,data)
begin

  plot = create "plot" vectorPlotClass wks
    "vcVectorFieldData" : data
  end create

  return(plot)

end


undef("setval_obj")
procedure setval_obj(data,res_string,res_value)
begin

  str_xcpt = (/"cnLevelSelectionMode","mpDataBaseVersion","cnLineLabelPlacementMode","mpGridAndLimbDrawOrder","mpOutlineDrawOrder","mpFillDrawOrder","mpProjection","mpGeophysicalLineColor","mpNationalLineColor","mpUSStateLineColor","mpGridLineColor","mpOutlineBoundarySets","mpFillColor","mpLimitMode","cnFillDrawOrder","cnFillPalette","cnLineLabelFontColor","cnLineColor","lbTitlePosition","lbTitleDirection","vcGlyphStyle","mpLandFillColor","mpOceanFillColor","mpInlandWaterFillColor","mpDefaultFillColor"/)
  strarr_xcpt = (/"mpOutlineSpecifiers"/)

  custom_res = (/"cnFillOpacities","cnScalarSmoothMe"/)

  if(any(res_string.eq.custom_res))then
    if(res_string.eq."cnFillOpacities")then
      getvalues data
        "cnFillPalette" : col_list
      end getvalues
      ;print(col_list(:,0)+"-"+col_list(:,1)+"-"+col_list(:,2))
      ;print(res_value+"")
      col_list_alpha = new((/dimsizes(col_list(:,0)),4/),float)
      col_list_alpha(:,0:2)=col_list
      num_res_values = str_fields_count(res_value,",")
      if(num_res_values.gt.1)then
        do jjj=0,num_res_values-1
          col_list_alpha(jjj,3) = tofloat(str_get_field(res_value,jjj+1,","))
        end do
      else
        col_list_alpha(:,3) = tofloat(res_value)
      end if 
      ;print(col_list_alpha(:,0)+"-"+col_list_alpha(:,1)+"-"+col_list_alpha(:,2)+"-"+col_list_alpha(:,3))
      setvalues data
        "cnFillPalette" : col_list_alpha
      end setvalues
    else if(res_string.eq."cnScalarSmoothMe")
      if(res_value.eq."True") then
        getvalues data
          "cnScalarFieldData" : data_plot
        end getvalues
        getvalues data_plot
          "sfDataArray" : data_to_smooth
          "sfMissingValueV" : aFillValue
        end getvalues
        data_to_smooth@_FillValue = aFillValue
        smooth_data = smooth_gempak(data_to_smooth)
        setvalues data_plot
          "sfDataArray" : smooth_data
        end setvalues
      end if
    else
      print("Unsupported custom resource: "+res_string)
    end if
    end if
  else ; for custom res if statement

  getvalues data
    res_string : def_val
  end getvalues

;  print(res_string+"")
  
  ; some resources convert string inputs to numbers and must be converted
  if(any(res_string.eq.str_xcpt)) then
    def_val := "default_string"
  else if(any(res_string.eq.strarr_xcpt)) then
    def_val := new(2,string)
  end if
  end if

  thistype = typeof(def_val)
  if(dimsizes(def_val).gt.1) then ;if (is_array)
    numvals = str_fields_count(res_value,",")
    thisvalue = new(numvals,thistype)
    if(thistype.eq."float") then
      do iii=1,numvals
        thisvalue(iii-1) = tofloat(str_get_field(res_value,iii,","))
      end do
    else if(thistype.eq."integer") then
      do iii=1,numvals
        thisvalue(iii-1) = toint(str_get_field(res_value,iii,","))
      end do
    else ;string and logical
      do iii=1,numvals
        thisvalue(iii-1) = str_get_field(res_value,iii,",")
      end do
    end if
    end if
  else ; *** if NOT an aray ***
     if(thistype.eq."float") then
       thisvalue = tofloat(res_value)
     else if(thistype.eq."integer") then
       thisvalue = toint(res_value)
     else ; string and logical
       thisvalue = res_value
     end if
     end if
  end if

  setvalues data
    res_string : thisvalue
  end setvalues

  end if ;for custom things

end

undef("calc_date_gfscycle")
function calc_date_gfscycle(timeUnits)
begin
  if(systemfunc("uname").eq."Darwin")then
    datecomm="gdate"
  else
    datecomm="date"
  end if

  year = stringtointeger( systemfunc( datecomm+" -u +'%Y' -d today " ) )
  month = stringtointeger( systemfunc( datecomm+" -u +'%m' -d today " ) )
  day = stringtointeger( systemfunc( datecomm+" -u +'%d' -d today " ) )
  inHour = stringtointeger( systemfunc( datecomm+" -u +'%H' -d today " ) )

  ;Assumes script is being run at least :30 minutes after start of cycle
  if( ( inHour.lt.3 ) ) then
    cycle = 18
  end if
  if( ( inHour.ge.3 ).and.( inHour.lt.9 ) ) then
    cycle = 00
  end if
  if( ( inHour.ge.9 ).and.( inHour.lt.15 ) ) then
    cycle = 06
  end if
  if( ( inHour.ge.15 ).and.( inHour.lt.21 ) ) then
    cycle = 12
  end if
  if( ( inHour.ge.21 ) ) then
    cycle = 18
  end if

  ;timeUnits = "days since 1800-01-01 00:00:00"
  currDate = cd_inv_calendar( year, month, day, cycle, 0, 0, timeUnits, 0 )
  if( inHour.lt.3 ) then
    currDate = currDate - 24
  end if

  currDate@units = timeUnits

  return(currDate)

end

undef("check_active_logs")
function check_active_logs(log_list)
begin

  active_inds = ind(.not.ismissing(log_list))
  if(.not.all(ismissing(log_list)))then
  do i=0,dimsizes(active_inds)-1
    ;this_log := asciiread(log_list(active_inds(i)),-1,"string")
    ;last_line = str_strip(this_log(dimsizes(this_log)-1))
    last_line = systemfunc("tail -1 "+log_list(active_inds(i)))
    if(.not.ismissing(last_line))then
    if(last_line.eq."PROD-DONE.")then
      system("rm -f "+log_list(active_inds(i)))
      log_list(active_inds(i)) = log_list@_FillValue
    end if
    end if
  end do
    return(dimsizes(ind(.not.ismissing(log_list))))
  else
    return(0)
  end if
end

undef("read_product_file")
function read_product_file(filename,return_var)
begin
  dynamaps_root = "/Users/kgriffin/dynamaps/"
  small_list = (/"usalaska","usmidwest","usnortheast","usnorthplains","ussoutheast","ussouthplains","uswest"/)
  small_list = dynamaps_root+small_list+".dyr"

  if(fileexists(filename))then
    prod = asciiread(filename,-1,"string")
    if(return_var.eq."num_overlays")then
      return(toint(prod(0)))
    end if
    if(return_var.eq."plot_name")then
      return(str_get_field(prod(1),1,","))
    end if
    if(return_var.eq."start_fhr")then
      return(toint(str_get_field(prod(1),2,",")))
    end if
    if(return_var.eq."end_fhr")then
      return(toint(str_get_field(prod(1),3,",")))
    end if
    if(return_var.eq."file_type")then
      return(str_get_field(prod(1),4,","))
    end if
    if(return_var.eq."regions")then
      if(prod(2).eq."all") then
        regions_list = systemfunc("ls "+dynamaps_root+"*.dyr -1")
      else if(prod(2).eq."all_us".or.prod(2).eq."us_all")then
        regions_list_all = systemfunc("ls "+dynamaps_root+"us*.dyr -1")
        regions_list = new(dimsizes(regions_list_all)+1,string)
        regions_list(0:dimsizes(regions_list_all)-1) = regions_list_all
        regions_list(dimsizes(regions_list_all)) = dynamaps_root+"conus.dyr"
        delete(regions_list_all)
      else if(prod(2).eq."all_no_small") then
        regions_list_all = systemfunc("ls "+dynamaps_root+"*.dyr -1")
        ind_list = new(dimsizes(regions_list_all),logical)
        do ii=0,dimsizes(regions_list_all)-1
          if(all(regions_list_all(ii).ne.small_list)) then
            ind_list(ii) = True
          end if
        end do
        inds = ind(ind_list)
        regions_list = regions_list_all(inds)
        print(regions_list)
      else
        nregs = str_fields_count(prod(2),",")
        regions_list = new(nregs,"string")
        do rrr=1,nregs
          regions_list(rrr-1) = dynamaps_root+str_get_field(prod(2),rrr,",")+".dyr"
        end do
      end if
      end if
      end if
      return(regions_list)
    end if
    if(return_var.eq."resources")then
      return(prod)
    end if
  else
   return("Missing file")
  end if
end

undef("get_prodlistfile")
function get_prodlistfile(fileroot)
begin
  return fileroot+"unified_products.txt"
end

undef("get_proclogfile")
function get_proclogfile(fileroot)
begin
  return fileroot+"proc_log.txt"
end

undef("init_products")
procedure init_products(proc_id,dynamaps_root)
begin
  if(floor(tofloat(proc_id)).eq.0)then
    listfile = get_prodlistfile(dynamaps_root) 
    logfile = get_proclogfile(dynamaps_root)
    system("rm -f "+listfile)
    system("rm -f "+logfile)
    prod_list=systemfunc("ls -1 "+dynamaps_root+"*.dyn")
    fhr_list = ispan(0,192,6)
    importance = 5 
    do m=0,dimsizes(fhr_list)-1
      do mm=0,dimsizes(prod_list)-1
        system("echo "+prod_list(mm)+","+fhr_list(m)+","+importance+" >> "+listfile) 
      end do
    end do
  end if
  ; procedure, no return.
end

undef("take_product")
function take_product(proc_id,dynamaps_root)
begin
  listfile = get_prodlistfile(dynamaps_root) 
  logfile = get_proclogfile(dynamaps_root)
  ; if not the first product of the run
  if(fileexists(logfile))then
    ; does this proc have something?
    this_proc = systemfunc("cat "+logfile+" | grep proc"+proc_id)
    ; only get new product if this proc_id doesn't have something else
    if(ismissing(this_proc))then
      ; read next (top) line of product file
      next_prod = systemfunc("sed -n 1p "+listfile)
      if(ismissing(next_prod))then
        next_prod = "0" ;return a 0 when the list is empty
      end if
    else
      ; set next_prod as missing bc this proc hasn't returned last prod
      next_prod = new(1,string)
    end if
  else
    ; if it is the first product, just pull the first line and delete
    next_prod = systemfunc("sed -n 1p "+listfile)
    system("sed -i 1d "+listfile)
  end if
  if(.not.ismissing(next_prod))then
    ;if there is a product to return, log it with the proc_id
    ;then delete it from the list of things still to do
    system("echo proc"+proc_id+": "+next_prod+" >> "+logfile)
    system("sed -i 1d "+listfile)
    print(next_prod+"")
  end if
  return(next_prod)
end

undef("return_product")
procedure return_product(proc_id,ret_status,dynamaps_root)
begin
  ; ret_status 0 = success
  ; ret_status {log file number} = failure of some sort...
  listfile = get_prodlistfile(dynamaps_root)
  logfile = get_proclogfile(dynamaps_root)
  this_proc = systemfunc("cat "+logfile+" | grep proc"+proc_id)
  ;proc_entry = systemfunc("sed -nr \'/proc|:/p\' "+logfile)
  ret_product = str_sub_str(this_proc,"proc"+proc_id+": ","")
  ;**remove proc entry from the log file, process is done.
  system("sed -i '/proc"+proc_id+"/d' "+logfile)
  if(ret_status.ne.0)then
    ; if not successful, add back to top of product list to run again.
    system("sed -i '1i"+ret_product+"' "+listfile)
  end if
end

undef("parse_product")
function parse_product(the_prod,want)
begin
  ; current return options are:
  ; field, name
  ; 0, name
  ; 1, fhr
  ; 2, imp
  splitprod = str_split_csv(the_prod,",",3)
  if(want.eq."fhr")then
    return(toint(splitprod(0,1))) 
  end if
  if(want.eq."name")then
    return(splitprod(0,0))
  end if
  if(want.eq."imp")then
    return(toint(splitprod(0,2)))
  end if
end
