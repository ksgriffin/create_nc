#!/bin/tcsh

set log_loc = "/nialldata/ksgriffin2/realtime/flogs/"
#set data_loc = "/nialldata/ksgriffin2/data/new_gfs_nc/"
set data_loc = "/weather/data/nc/"

set fhr_list = {'000','006','012','018','024','030','036','042','048','054','060','066','072','078','084','090','096','102','108','114','120','126','132','138','144','150','156','162','168','174','180','186','192','198','204','210','216','222','228','234','240'}

set hour = `date -u +%k`

if ($hour < 3) then
  set cycle = 18;
else if ($hour < 9) then
  set cycle = 00;
else if ($hour < 15) then
  set cycle = 06;
else if ($hour < 21) then
  set cycle = 12;
else
  set cycle = 18;
endif

set fdate = `date -u +%y%m%d`
set rmdate = `date -u --date "3 days ago" +%y%m%d`
set rmanldate = `date -u --date "15 days ago" +%y%m%d`
set rmstring = "rm -f "
set rmlogstring = "rm -f"

#set rmdate = "150126"
#set cycle = "06"

echo Clearing log files in ${log_loc}
rm -f ${log_loc}/f*.txt

foreach f (${fhr_list})
  if($f == 000)then
    echo Clearing analysis file from ${rmanldate}_${cycle}
    rm -f ${data_loc}gfs${rmanldate}_${cycle}_F000.nc
    rm -f ${data_loc}derivedgfs${rmanldate}_${cycle}_F000.nc
  else
    echo Removing FHR ${f} from ${rmdate}_${cycle}
    rm -f ${data_loc}gfs${rmdate}_${cycle}_F${f}.nc ${data_loc}derivedgfs${rmdate}_${cycle}_F${f}.nc
  endif
  rm -f ${log_loc}${rmdate}_${cycle}_f${f}.txt
end
