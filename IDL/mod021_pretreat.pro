PRO mod021_pretreat
  COMPILE_OPT idl2
  ENVI, /RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT, LOG_FILE = 'batch_log.txt'
  
  file = file_search('C:\modis\source\*.hdf') ; search all hdf files which need pretreatment
  print,((size(file))[3]-1)
  for i = 0, ((size(file))[3]-1) do begin     ; batch processing
    print,i
    b = strsplit(file[i],'.',/extract)        ; get the file name
    c = strsplit(b[0],'\',/extract)
    name = c[3]
    dict_tiff_name = 'C:\modis\tiff\'  ; define where to read or save file
    dict_rgb_name = 'C:\modis\rgb\'
    dict_temp_name = 'C:\modis\temp\'
    source_dict = 'C:\modis\source\'
    modisname = source_dict + name + '.hdf'
    outname_radiation = dict_temp_name + name + '_radiation'
    outname_angle = dict_temp_name + name + '_angle'
    outname_data = dict_temp_name + name + '_data'
    GCP1 = dict_temp_name + name + '_gcp1.dat'
    GCP2 = dict_temp_name + name + '_gcp2.dat'
    relon = dict_temp_name + name + '_relon.dat'
    relat = dict_temp_name + name + '_relat.dat'
    outname_stack = dict_temp_name + name + '_stack'
    outfile = dict_temp_name + name + '_final'
    tiff_name = dict_tiff_name + name + '.tif'
    out_name = dict_temp_name + name + '_rgb'
    rgb_name = dict_rgb_name + name + '.tif'
    
    RADIATION_CORRECTION,modisname,outname_radiation
    MODIS_REGISTER,modisname,outname_radiation,outname_angle,outname_data,GCP1,GCP2,relon,relat,outname_stack
    convert,outname_stack, outfile
    
    ; save the file in tiff format
    ENVI_OPEN_GDB,outfile,r_fid=cut_fid
    ENVI_FILE_QUERY, cut_fid,dims=dims, nb = nb
    pos = LINDGEN(nb)
    ENVI_OUTPUT_TO_EXTERNAL_FORMAT, FID=cut_fid,Dims=dims,POS=pos,Out_Name=tiff_name,/TIFF
    
    a = rgb_extract(outfile,out_name,rgb_name)
    print,'finish'+name
  endfor
  
  print,'finish all'
  ENVI_BATCH_EXIT
END

;;--------------------------------rgb----------------------------
FUNCTION rgb_extract,outfile,out_name,rgb_name

  ENVI_OPEN_GDB,outfile,r_fid=data_fid
  ENVI_FILE_QUERY, data_fid,dims=dims,nb=nb
  OUT_BNAME = ['R(band1) [250 Aggr]', 'G(band4) [500 Aggr]', 'B(band3) [500 Aggr]']
  ; define which bands to stack and their corresponding band names
  ENVI_FILE_QUERY, data_fid,dims=dims,nb=nb
  ENVI_DOIT, 'cf_doit',fid=[data_fid,data_fid,data_fid],$
    pos=[0,3,2], dims=dims, OUT_BNAME=OUT_BNAME,$
    remove=0,r_fid=r_fid,out_name=out_name
  print, 'finish extract'

  ; save the file in tiff format
  ENVI_FILE_QUERY, r_fid,dims=dims, nb = nb
  pos = LINDGEN(nb)
  ENVI_OUTPUT_TO_EXTERNAL_FORMAT, FID=r_fid,Dims=dims,POS=pos,Out_Name=rgb_name,/TIFF
  return,1
END


;; ------------------------------------cut-----------------------------------
PRO convert,outname_stack, outfile
  COMPILE_OPT idl2
  ENVI, /RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT, LOG_FILE = 'batch_log.txt'
  
  ENVI_OPEN_GDB,outname_stack,r_fid=data_fid
  ENVI_FILE_QUERY, data_fid,dims=dims,nb=nb
  pos = LINDGEN(nb)

  xmap = [103.6667,105.6667]  ; longitude
  ymap = [30.6667,29.1667]    ; latitude

  envi_convert_file_coordinates,data_fid,xf,yf,xmap,ymap

  dims = [-1,long(xf[0]),long(xf[1])-1,long(yf[0]),long(yf[1])-1]

  OUT_BNAME=['1KM Reflectance (band1) [250 Aggr]','1KM Reflectance (band2) [250 Aggr]',$
    '1KM Reflectance (band3) [500 Aggr]','1KM Reflectance (band4) [500 Aggr]',$
    '1KM Reflectance (band5) [500 Aggr]','1KM Reflectance (band6) [500 Aggr]',$
    '1KM Reflectance (band7) [500 Aggr]','1KM Reflectance (band17)',$
    '1KM Reflectance (band18)','1KM Reflectance (band19)',$
    '1KM Emissive (band24)','1KM Emissive (band25)',$
    'Warp(SensorZenith)','Warp(SensorAzimuth)',$
    'Warp(SolarZenith)','Warp(SolarAzimuth)']
  ; cut the data based on the longitude and latitude of the target area
  envi_doit,'RESIZE_DOIT',fid = data_fid,dims = dims,$
    out_name = outfile,pos = pos,out_bname = OUT_BNAME,$
    r_fid = cut_fid,rfact = [1.0,1.0]
  print,'finish cut'
  ENVI_BATCH_EXIT
end


;;;----------------------------------------------------------------------------------
;;=======================================================================
;;;  register
;;;======================================================================
PRO MODIS_REGISTER,modisname,outname_radiation,outname_angle,outname_data,GCP1,GCP2,relon,relat,outname_stack
  COMPILE_OPT idl2
  ENVI, /RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT, LOG_FILE = 'batch_log.txt'
  
  ;read longitude and latitude information
  Lon=READ_DATASET(modisname,'Longitude')
  Lat=READ_DATASET(modisname,'Latitude')
  ; resample the longitude and latitude information
  ENVI_WRITE_ENVI_FILE,Lon[*,*],r_fid=fid_Lon,out_name=relon
  FIDlon= MODIS_RESIZE(fid_Lon)
  ENVI_FILE_QUERY, FIDlon,dims=dimslon
  print,dimslon
  datalon = ENVI_GET_DATA(fid=FIDlon, dims=dimslon, pos=0)

  ENVI_WRITE_ENVI_FILE,Lat[*,*],r_fid=fid_Lat,out_name=relat
  FIDlat= MODIS_RESIZE(fid_Lat)
  ENVI_FILE_QUERY, FIDlat,dims=dimslat
  datalat = ENVI_GET_DATA(fid=FIDlat, dims=dimslat, pos=0)
  ;;;=======================================================================
  ;;;  set GCP
  ;;;======================================================================
  OPENW,lun,GCP1,/Get_lun
  PRINTF,lun,' gcp[0,l] ','gcp[1,l] ','gcp[2,l]  ','gcp[2,l]  ','n  ','m  '
  length=0L
  ; set GCP in the shape of 51 * 51
  gcp= DBLARR(4,2601)
  rstr = ["Input File :" , "Output File :"]
  ENVI_REPORT_INIT,rstr,title="output the GCP", base=base,/INTERRUPT
  FOR i=3.5d,2030,40.53d DO BEGIN
    FOR j=3.5d,1354,27.01d DO BEGIN
      gcp[0,length]=datalon[j,i]
      gcp[1,length]=datalat[j,i]
      gcp[2,length]=j
      gcp[3,length]=i
      PRINTF,lun,gcp[0,length],gcp[1,length],gcp[2,length],gcp[3,length]
      length=length+1L
    ENDFOR
    ENVI_REPORT_STAT, base, i*1354, 2601
  ENDFOR
  ENVI_REPORT_INIT, base=base, /finish
  FREE_LUN,lun
  ;;;=======================================================================
  ;;;  GCP projection transformation
  ;;;======================================================================
  OPENW,lun,GCP2,/Get_lun
  PRINTF,lun,'ogcp[0,l]','ogcp[1,l] ','ogcp[2,l]','ogcp[2,l]'
  iproj= ENVI_PROJ_CREATE(/geographic)
  units = envi_translate_projection_units('Degrees')
  datum = 'WGS-84'
  oproj= ENVI_PROJ_CREATE(/utm,zone=48,datum=datum,units=units)  ; zone=48
  ENVI_CONVERT_PROJECTION_COORDINATES, gcp[0,*], gcp[1,*],iproj,outx, outy, oproj
  ogcp= DBLARR(4,2601)
  length=0L
  FOR i=3.5d,2030,40.53d DO BEGIN
    FOR j=3.5d,1354,27.01d DO BEGIN
      ogcp[0,length]=outx[[(i-3.5)/40.53]*51+[(j-3.5)/27.01]]
      ogcp[1,length]=outy[[(i-3.5)/40.53]*51+[(j-3.5)/27.01]]
      ogcp[2,length]=gcp[2,length]
      ogcp[3,length]=gcp[3,length]
      PRINTF,lun,ogcp[0,length],ogcp[1,length],ogcp[2,length],ogcp[3,length]
      length=length+1L
    ENDFOR
  ENDFOR
  FREE_LUN,lun
  ;;;=======================================================================
  ;;; start registering
  ;;;======================================================================
  ; angle registering
  OUT_BNAME=['Warp(SensorZenith)','Warp(SensorAzimuth)',$
    'Warp(SolarZenith)','Warp(SolarAzimuth)']
  anglefid=ANGLE_DATA(modisname) ; angle data correction
  ENVI_FILE_QUERY,anglefid,nb=nb_angle,dims=dims_angle
  reg_id1 = REGISTER(anglefid,nb_angle,dims_angle,ogcp,oproj,outname_angle,OUT_BNAME)

  ; data registering
  OUT_BNAME=['1KM Reflectance (band1) [250 Aggr]','1KM Reflectance (band2) [250 Aggr]',$
    '1KM Reflectance (band3) [500 Aggr]','1KM Reflectance (band4) [500 Aggr]',$
    '1KM Reflectance (band5) [500 Aggr]','1KM Reflectance (band6) [500 Aggr]',$
    '1KM Reflectance (band7) [500 Aggr]','1KM Reflectance (band17)',$
    '1KM Reflectance (band18)','1KM Reflectance (band19)',$
    '1KM Emissive (band24)','1KM Emissive (band25)']

  ;read radiation data
  IF (ENVI_IS_GDB(outname_radiation)) THEN BEGIN
    ENVI_OPEN_GDB,outname_radiation,r_fid=data_fid
  ENDIF ELSE BEGIN
    ENVI_OPEN_FILE,outname_radiation,r_fid=data_fid
  ENDELSE
  ENVI_FILE_QUERY, data_fid, dims=dims_data, nb=nb_data
  reg_id2 = REGISTER(data_fid,nb_data,dims_data,ogcp,oproj,outname_data,OUT_BNAME)

;-------------------layer stacking---------------------
  OUT_BNAME=['1KM Reflectance (band1) [250 Aggr]','1KM Reflectance (band2) [250 Aggr]',$
    '1KM Reflectance (band3) [500 Aggr]','1KM Reflectance (band4) [500 Aggr]',$
    '1KM Reflectance (band5) [500 Aggr]','1KM Reflectance (band6) [500 Aggr]',$
    '1KM Reflectance (band7) [500 Aggr]','1KM Reflectance (band17)',$
    '1KM Reflectance (band18)','1KM Reflectance (band19)',$
    '1KM Emissive (band24)','1KM Emissive (band25)',$
    'Warp(SensorZenith)','Warp(SensorAzimuth)',$
    'Warp(SolarZenith)','Warp(SolarAzimuth)']
  ENVI_FILE_QUERY, reg_id2,dims=dims,nb=nb
  ; define which bands to stack and their corresponding band names
  ENVI_DOIT, 'cf_doit',fid=[reg_id2,reg_id2,reg_id2,reg_id2,reg_id2,reg_id2,reg_id2,$
    reg_id2,reg_id2,reg_id2, $
    reg_id2,reg_id2,reg_id1,reg_id1,reg_id1,reg_id1],$
    pos=[0,1,2,3,4,5,6,7,8,9,10,11,0,1,2,3], dims=dims, OUT_BNAME=OUT_BNAME,$
    remove=0,r_fid=stack_fid,out_name=outname_stack
    
  print,'finish register'  
  ENVI_BATCH_EXIT
END



FUNCTION  ANGLE_DATA,modisname
  COMPILE_OPT idl2
  ;read angle information
  sensorzenith=READ_DATASET(modisname,'SensorZenith')
  sensorazimuth=READ_DATASET(modisname,'SensorAzimuth')
  solarzenith=READ_DATASET(modisname,'SolarZenith')
  solarazimuth=READ_DATASET(modisname,'SolarAzimuth')
  ;;;=======================================================================
  ;;;  angle information integration
  ;;;======================================================================
  ENVI_WRITE_ENVI_FILE,sensorzenith[*,*],r_fid=fid_sez,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,sensorazimuth[*,*],r_fid=fid_sea,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,solarzenith[*,*],r_fid=fid_soz,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,solarazimuth[*,*],r_fid=fid_soa,/IN_MEMORY
  ENVI_FILE_QUERY,fid_sez,dims=sezdims
  fid_sez_c=ANGLE_CALCULATE_Z(fid_sez) ; transform the angle data into the range of (0,1)
  fid_sea_c=ANGLE_CALCULATE_A(fid_sea)
  fid_soz_c=ANGLE_CALCULATE_Z(fid_soz)
  fid_soa_c=ANGLE_CALCULATE_A(fid_soa)
  ENVI_DOIT, 'cf_doit',fid=[fid_sez_c,fid_sea_c,fid_soz_c,fid_soa_c], $
    pos=[0,0,0,0], dims=sezdims, $
    remove=0,r_fid=info_fid ,/in_memory
  ; release memory
  ENVI_FILE_MNG, id=fid_sez, /remove
  ENVI_FILE_MNG, id=fid_sea, /remove
  ENVI_FILE_MNG, id=fid_soz, /remove
  ENVI_FILE_MNG, id=fid_soa, /remove
  angle_fid= MODIS_RESIZE(info_fid)  ; resize angle information

  RETURN,angle_fid
END


; Band Math to calculate and transform the angle information
FUNCTION ANGLE_CALCULATE_A,fid
  COMPILE_OPT idl2

  ENVI_FILE_QUERY, fid, dims=dims
  t_fid = [fid]
  pos = [0]
  ; expression
  exp='((b1*0.0100)+180)/360'
  ENVI_DOIT, 'math_doit', $
    fid=t_fid, pos=pos, dims=dims,$
    exp=exp,r_fid=r_fid, /IN_MEMORY
  ENVI_FILE_MNG, id=fid, /remove
  RETURN,r_fid
END

FUNCTION ANGLE_CALCULATE_Z,fid
  COMPILE_OPT idl2

  ENVI_FILE_QUERY, fid, dims=dims
  t_fid = [fid]
  pos = [0]
  ; expression
  exp='b1*0.0100/90'
  ENVI_DOIT, 'math_doit', $
    fid=t_fid, pos=pos, dims=dims,$
    exp=exp,r_fid=r_fid, /IN_MEMORY
  ENVI_FILE_MNG, id=fid, /remove
  RETURN,r_fid
END


FUNCTION READ_DATASET,modisname,dataset
  SD_id = HDF_SD_START(modisname)
  index = HDF_SD_NAMETOINDEX(SD_id,dataset)
  sds_id= HDF_SD_SELECT(SD_id,index)
  HDF_SD_GETINFO,sds_id,DIMS=dim
  HDF_SD_GETDATA,sds_id,data
  HDF_SD_ENDACCESS,sds_id
  RETURN,data
END


FUNCTION MODIS_RESIZE,fid
  COMPILE_OPT IDL2
  ; Open the input file
  IF (fid EQ -1) THEN BEGIN
    RETURN,-1
  ENDIF
  ENVI_FILE_QUERY, fid, dims=dims, nb=nb
  pos = LINDGEN(nb)
  out_name = 'testimg'
  ;resample the data information to 1354*2030
  ENVI_DOIT, 'resize_doit', $
    fid=fid, pos=pos, dims=dims, $
    interp=1, rfact=[.20015,.2], $
    out_name=out_name, r_fid=r_fid
  RETURN,r_fid
END

FUNCTION REGISTER,fid,nb,dims,ogcp,oproj,out_name,OUT_BNAME
  COMPILE_OPT idl2
  pos = LINDGEN(nb)
  pixel_size = [0.008333,0.008333]
  ENVI_DOIT, 'envi_register_doit', w_fid=fid, w_pos=pos,w_dims=dims,method=7,$
    out_name=out_name,pts=ogcp, proj=oproj,r_fid=r_fid,pixel_size=pixel_size,OUT_BNAME=OUT_BNAME
  ;ENVI_FILE_MNG, id=r_fid3, /remove
  ENVI_FILE_MNG, id=fid, /remove
  RETURN,r_fid
END




;;;--------------------------------------------------------------------------------------
;;;=======================================================================
;;; Radition Correction of specified bands data
;   Band：  1-7    Aerosol
;          17-19  Water Vapor
;          24-25  Temperature
;;;======================================================================
PRO RADIATION_CORRECTION,modisname,outname_radiation
  COMPILE_OPT idl2
  ENVI, /RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT, LOG_FILE = 'batch_log.txt'

  ref250=READ_DATASET_1(modisname,'EV_250_Aggr1km_RefSB');read band1 and band2
  ref500=READ_DATASET_1(modisname,'EV_500_Aggr1km_RefSB');read band3-band7
  ref1000=READ_DATASET_1(modisname,'EV_1KM_RefSB')       ;read band8-band19 and band26
  Emi1000=READ_DATASET_1(modisname,'EV_1KM_Emissive')    ;read band20-band25 and band27-band36
  ;;;=======================================================================
  ;;; radiation information intergration 
  ;;;======================================================================
  ENVI_WRITE_ENVI_FILE,ref250[*,*,0],r_fid=fid_ref250_1,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,ref250[*,*,1],r_fid=fid_ref250_2,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,ref500[*,*,0],r_fid=fid_ref500_3,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,ref500[*,*,1],r_fid=fid_ref500_4,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,ref500[*,*,2],r_fid=fid_ref500_5,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,ref500[*,*,3],r_fid=fid_ref500_6,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,ref500[*,*,4],r_fid=fid_ref500_7,/IN_MEMORY

  ENVI_WRITE_ENVI_FILE,ref1000[*,*,11],r_fid=fid_ref1000_17,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,ref1000[*,*,12],r_fid=fid_ref1000_18,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,ref1000[*,*,13],r_fid=fid_ref1000_19,/IN_MEMORY

  ENVI_WRITE_ENVI_FILE,Emi1000[*,*,4],r_fid=fid_Emi1000_24,/IN_MEMORY
  ENVI_WRITE_ENVI_FILE,Emi1000[*,*,5],r_fid=fid_Emi1000_25,/IN_MEMORY

  ;EV_250_Aggr1km_RefSB
  r_fid_ref250_1=REFLECT_RADIANCE_CALCULATE(modisname,fid_ref250_1,6,8,9,0,0)
  r_fid_ref250_2=REFLECT_RADIANCE_CALCULATE(modisname,fid_ref250_2,6,8,9,1,1)
  ;EV_500_Aggr1km_RefSB
  r_fid_ref500_3=REFLECT_RADIANCE_CALCULATE(modisname,fid_ref500_3,9,8,9,0,0)
  r_fid_ref500_4=REFLECT_RADIANCE_CALCULATE(modisname,fid_ref500_4,9,8,9,1,1)
  r_fid_ref500_5=REFLECT_RADIANCE_CALCULATE(modisname,fid_ref500_5,9,8,9,2,2)
  r_fid_ref500_6=REFLECT_RADIANCE_CALCULATE(modisname,fid_ref500_6,9,8,9,3,3)
  r_fid_ref500_7=REFLECT_RADIANCE_CALCULATE(modisname,fid_ref500_7,9,8,9,4,4)

  ;EV_1KM_RefSB
  r_fid_ref1000_17=REFLECT_RADIANCE_CALCULATE(modisname,fid_ref1000_17,2,8,9,11,11)
  r_fid_ref1000_18=REFLECT_RADIANCE_CALCULATE(modisname,fid_ref1000_18,2,8,9,12,12)
  r_fid_ref1000_19=REFLECT_RADIANCE_CALCULATE(modisname,fid_ref1000_19,2,8,9,13,13)

  ;EV_1KM_Emissive
  r_fid_Emi1000_24=REFLECT_RADIANCE_CALCULATE(modisname,fid_Emi1000_24,4,5,6,4,4)
  r_fid_Emi1000_25=REFLECT_RADIANCE_CALCULATE(modisname,fid_Emi1000_25,4,5,6,5,5)

  ENVI_FILE_QUERY, r_fid_ref250_1,dims=dims
  ; layer stack
  OUT_BNAME=['1KM Reflectance (band1) [250 Aggr]','1KM Reflectance (band2) [250 Aggr]',$
    '1KM Reflectance (band3) [500 Aggr]','1KM Reflectance (band4) [500 Aggr]',$
    '1KM Reflectance (band5) [500 Aggr]','1KM Reflectance (band6) [500 Aggr]',$
    '1KM Reflectance (band7) [500 Aggr]','1KM Reflectance (band17)',$
    '1KM Reflectance (band18)','1KM Reflectance (band19)',$
    '1KM Emissive (band24)','1KM Emissive (band25)']
  ENVI_DOIT, 'cf_doit',fid=[r_fid_ref250_1,r_fid_ref250_2,r_fid_ref500_3,r_fid_ref500_4,r_fid_ref500_5,r_fid_ref500_6,r_fid_ref500_7,$
    r_fid_ref1000_17,r_fid_ref1000_18,r_fid_ref1000_19,r_fid_Emi1000_24,r_fid_Emi1000_25], $
    pos=[0,0,0,0,0,0,0,0,0,0,0,0], dims=dims, OUT_BNAME=OUT_BNAME,$
    remove=0,r_fid=data_fid,out_name=outname_radiation
  print,'finish calculate'  
  ENVI_BATCH_EXIT
END

FUNCTION READ_DATASET_1,modisname,dataset
  SD_id = HDF_SD_START(modisname)
  index = HDF_SD_NAMETOINDEX(SD_id,dataset)
  sds_id= HDF_SD_SELECT(SD_id,index)
  HDF_SD_GETDATA,sds_id,data
  HDF_SD_ENDACCESS,sds_id
  RETURN,data
END

FUNCTION REFLECT_RADIANCE_CALCULATE,modisname,fid,dataset_number,scalesattr_index,offsetsattr_index,scales_index,offsets_index
  COMPILE_OPT idl2
  SD_id = HDF_SD_START(modisname,/Read)
  SdsID=HDF_SD_SELECT(SD_id,dataset_number)
  HDF_SD_ATTRINFO,SdsID,scalesattr_index,NAME=re_scales,DATA=scalesData
  HDF_SD_ATTRINFO,SdsID,offsetsattr_index,NAME=re_offsets,DATA=offsetsData

  ;Reflectance: R=reflectance_scale*(DN-reflectance _offset)
  ;Radiance:    R=radiance_scale*(DN-radiance _offset)
  ENVI_FILE_QUERY, fid, dims=dims
  t_fid = [fid]
  pos = [0]
  ;expression
  exp=STRTRIM(scalesData[scales_index],2)+'*[b1-('+STRTRIM(offsetsData[offsets_index],2)+')]'
  ENVI_DOIT, 'math_doit', $
    fid=t_fid, pos=pos, dims=dims,$
    exp=exp,r_fid=r_fid, /IN_MEMORY
  ENVI_FILE_MNG, id=fid, /remove
  RETURN,r_fid
END


