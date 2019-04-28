PRO mod04_pretreat
  COMPILE_OPT idl2
  ENVI, /RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT, LOG_FILE = 'batch_log.txt'
  
  file = file_search('C:\test\source\*.hdf')  ; search all hdf files which need pretreatment
  print,((size(file))[3]-1)
  for i = 0, ((size(file))[3]-1) do begin     ; batch processing
    print,i
    b = strsplit(file[i],'.',/extract) ; get the name of the file
    c = strsplit(b[0],'\',/extract)
    name = c[3]
    dict_tiff_name = 'C:\test\near\'   ; the dict to store tiff files
    dict_temp_name = 'C:\test\temp\'   ; the dict to store mid files
    source_dict = 'C:\test\source\'    ; the dict to read source file
    ; which file to save and where to save
    modisname = source_dict + name + '.hdf'
    outname_aod = dict_temp_name + name + '_aod'
    outname_r_aod = dict_temp_name + name + '_r_aod'
    GCP1 = dict_temp_name + name + '_gcp1.dat'
    GCP2 = dict_temp_name + name + '_gcp2.dat'
    relon = dict_temp_name + name + '_relon.dat'
    relat = dict_temp_name + name + '_relat.dat'
    outfile = dict_temp_name + name + '_final'
    tiff_name = dict_tiff_name + name + '.tif'

    READ_MOD04,modisname,outname_aod   ; read file
    ; register
    MODIS_REGISTER,modisname,outname_aod,outname_r_aod,GCP1,GCP2,relon,relat
    convert,outname_r_aod, outfile     ; cut the file

    ; save the file as tiff format
    ENVI_OPEN_GDB,outfile,r_fid=cut_fid
    ENVI_FILE_QUERY, cut_fid,dims=dims, nb = nb
    pos = LINDGEN(nb)
    ENVI_OUTPUT_TO_EXTERNAL_FORMAT, FID=cut_fid,Dims=dims,POS=pos,Out_Name=tiff_name,/TIFF

    print,'finish'+name
  endfor

  print,'finish all'
  ENVI_BATCH_EXIT
END

;; ------------------------------------cut-----------------------------------
PRO convert,outname_r_aod, outfile
  COMPILE_OPT idl2
  ENVI, /RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT, LOG_FILE = 'batch_log.txt'

  ENVI_OPEN_GDB,outname_r_aod,r_fid=data_fid
  ENVI_FILE_QUERY, data_fid,dims=dims,nb=nb
  pos = LINDGEN(nb)

  xmap = [103.6667,105.6667]   ; the area longitude
  ymap = [30.6667,29.1667]     ; the area latitude

  envi_convert_file_coordinates,data_fid,xf,yf,xmap,ymap

  dims = [-1,long(xf[0]),long(xf[1])-1,long(yf[0]),long(yf[1])-1]

  OUT_BNAME=['AOD']
  ; cut the data according to the longitude and latitude of the area
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
PRO MODIS_REGISTER,modisname,outname_aod,outname_r_aod,GCP1,GCP2,relon,relat
  COMPILE_OPT idl2
  ENVI, /RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT, LOG_FILE = 'batch_log.txt'

  ; read the longitude and latitude data
  Lon=READ_DATASET(modisname,'Longitude')
  Lat=READ_DATASET(modisname,'Latitude')
  ENVI_WRITE_ENVI_FILE,Lon[*,*],r_fid=fid_Lon,out_name=relon
  ENVI_FILE_QUERY, fid_Lon,dims=dimslon
  datalon = ENVI_GET_DATA(fid=fid_Lon, dims=dimslon, pos=0)
  ENVI_WRITE_ENVI_FILE,Lat[*,*],r_fid=fid_Lat,out_name=relat
  ENVI_FILE_QUERY, fid_Lat,dims=dimslat
  datalat = ENVI_GET_DATA(fid=fid_Lat, dims=dimslat, pos=0)
  ;;;=======================================================================
  ;;;  set GCP
  ;;;======================================================================
  OPENW,lun,GCP1,/Get_lun
  PRINTF,lun,' gcp[0,l] ','gcp[1,l] ','gcp[2,l]  ','gcp[2,l]  ','n  ','m  '
  length=0L
  ; set GCP in the shape of  51 * 51
  gcp= DBLARR(4,2601)
  rstr = ["Input File :" , "Output File :"]
  ENVI_REPORT_INIT,rstr,title="output the GCP", base=base,/INTERRUPT
  FOR i=0.0d,675,13.5d DO BEGIN
    FOR j=0.0d,450,9.0d DO BEGIN
      gcp[0,length]=datalon[j,i]
      gcp[1,length]=datalat[j,i]
      gcp[2,length]=j
      gcp[3,length]=i
      PRINTF,lun,gcp[0,length],gcp[1,length],gcp[2,length],gcp[3,length]
      length=length+1L
    ENDFOR
    ENVI_REPORT_STAT, base, i*451, 2601
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
  FOR i=1.0d,676,13.5d DO BEGIN
    FOR j=1.0d,451,9.0d DO BEGIN
      ogcp[0,length]=outx[[(i-1.0)/13.5]*51+[(j-1.0)/9.0]]
      ogcp[1,length]=outy[[(i-1.0)/13.5]*51+[(j-1.0)/9.0]]
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
  OUT_BNAME=['AOD']

  ; read data
  IF (ENVI_IS_GDB(outname_aod)) THEN BEGIN
    ENVI_OPEN_GDB,outname_aod,r_fid=data_fid
  ENDIF ELSE BEGIN
    ENVI_OPEN_FILE,outname_aod,r_fid=data_fid
  ENDELSE

  ENVI_FILE_QUERY, data_fid, dims=dims_data, nb=nb_data
  reg_id2 = REGISTER(data_fid,nb_data,dims_data,ogcp,oproj,outname_r_aod,OUT_BNAME)

  print,'finish register'
  ENVI_BATCH_EXIT
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

FUNCTION REGISTER,fid,nb,dims,ogcp,oproj,out_name,OUT_BNAME
  COMPILE_OPT idl2
  pos = LINDGEN(nb)
  pixel_size = [0.008333,0.008333] ; define the pixel size( here, the unit is degree)
  ENVI_DOIT, 'envi_register_doit', w_fid=fid, w_pos=pos,w_dims=dims,method=6,$
    out_name=out_name,pts=ogcp, proj=oproj,r_fid=r_fid,pixel_size=pixel_size,OUT_BNAME=OUT_BNAME
  ENVI_FILE_MNG, id=fid, /remove
  RETURN,r_fid
END




;;;--------------------------------------------------------------------------------------
;;;=======================================================================
;;; read mod04_1K file
;;;======================================================================
PRO READ_MOD04,modisname,outname_aod
  COMPILE_OPT idl2
  ENVI, /RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT, LOG_FILE = 'batch_log.txt'

  ; read data
  aod=READ_DATASET_1(modisname,'Optical_Depth_Land_And_Ocean')
  ENVI_WRITE_ENVI_FILE,aod[*,*],r_fid=fid_aod,/IN_MEMORY

  ENVI_FILE_QUERY, fid_aod,dims=dims
  ; stack and save data
  OUT_BNAME=['AOD']
  ENVI_DOIT, 'cf_doit',fid=fid_aod, $
    pos=[0], dims=dims, OUT_BNAME=OUT_BNAME,$
    remove=0,r_fid=data_fid,out_name=outname_aod
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


