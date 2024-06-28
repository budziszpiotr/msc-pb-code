echo off
set CORE=C:\STP_1.6\src\stp-core
set THREADS=6

echo Start time: %time%

REM ####### SET RELEVANT PARAMETERS FROM HERE #######

REM ####### Copy INPUT from E:/TDF_preview/TDF_preview/20215802_Dahl  #######

set INPUT=E:\STP\20215802_Dahl\TDF
set TMP=E:\STP\20215802_Dahl\tmp
set LOG=E:\STP\20215802_Dahl\LOG
set RECON=E:\STP\20215802_Dahl\RECONS

set FILENAMES[2]=rat_77_full_Z-3.30_X-3.00_Y1.83
set OFFSETS[2]=-1
set FILENAMES[3]=rat_77_full_Z-3.30_X-3.00_Y4.00
set OFFSETS[3]=0
set FILENAMES[4]=rat_77_full_Z-3.30_X-1.17_Y-2.50
set OFFSETS[4]=-2
set FILENAMES[5]=rat_77_full_Z-3.30_X-1.17_Y-0.33
set OFFSETS[5]=-1
set FILENAMES[6]=rat_77_full_Z-3.30_X-1.17_Y1.83
set OFFSETS[6]=-1
set FILENAMES[7]=rat_77_full_Z-3.30_X-1.17_Y4.00
set OFFSETS[7]=0
set FILENAMES[8]=rat_77_full_Z-3.30_X0.67_Y-2.50
set OFFSETS[8]=-2
set FILENAMES[9]=rat_77_full_Z-3.30_X0.67_Y-0.33
set OFFSETS[9]=-1
set FILENAMES[10]=rat_77_full_Z-3.30_X0.67_Y1.83
set OFFSETS[10]=0
set FILENAMES[11]=rat_77_full_Z-3.30_X0.67_Y4.00
set OFFSETS[11]=1
set FILENAMES[12]=rat_77_full_Z-3.30_X2.50_Y-2.50
set OFFSETS[12]=-2
set FILENAMES[13]=rat_77_full_Z-3.30_X2.50_Y-0.33
set OFFSETS[13]=-1
set FILENAMES[0]=rat_77_full_Z-3.30_X2.50_Y1.83
set OFFSETS[0]=0
set FILENAMES[1]=rat_77_full_Z-3.30_X2.50_Y4.00
set OFFSETS[1]=1

REM #################### TO HERE ####################
set START_PROJ=0
set END_PROJ=1799

REM PREPROCESSING
set AIR_LEFT=0
set AIR_RIGHT=0
REM TOMO360 
set EXT_FOV=False
set EXT_FOV_RIGHT=False
set EXT_FOV_OVERLAP=0
set LINE_BY_LINE=False
set AVG_OVERLAP=False
set EXT_FOV_PROJ=3600
REM RING REMOVAL 
set RR_METHOD="none:0;0"
REM "none", "rivers:11;0", "muench:3;1.8", ...
REM NORMALIZATION
set DYN_FF=False

REM PHASE RETRIEVAL
set PR_METHOD=0
REM 0: TIE-Paganin 2002, 1: Gen. TIE-Paganin 2020, 2: CTF-Moosmann 2011
set BETA=1.0E-09
set DELTA=0.8E-07
set ENERGY=20.7
set DISTANCE=180
set PIXEL=1.6
set OVERPAD=True

REM RECONSTRUCTION
set ANG_RANGE=3.14159265358979
set FILTER=ram-lak
set SCALE=1
set OVERPADDING=True
set LOGTRSF=False
set CIRCLE=True
set RECON_METHOD=FBP_CUDA
set DECIM_FACTOR=1
set DOWNSC_FACTOR=1
set ROLLING=True
set ROLL_SHIFT=0

REM POSTPROCESSING
echo Start processing: %time%

set POSTPROC=False
set POLAR_FILT="homomorphic:0.8;0.2"
set CONVERT_OPT="linear8:-0.01;0.01"
set CROP_OPT="0:0:0:0"


set x=0
:MainLoop
if not defined FILENAMES[%x%] GOTO :EndLoop
call set FILENAME=%%FILENAMES[%x%]%%
echo Starting file: %FILENAME%
if not defined OFFSETS[%x%] GOTO :FailedLoop
call set OFFSET=%%OFFSETS[%x%]%%
SET PADDED=0%X%
SET ID=%PADDED:~-2%
echo Preprocessing...start_time:%time%
C:\STP_1.6\python\python-2.7.10.amd64\python.exe %CORE%\exec_preprocessing.py %START_PROJ% %END_PROJ% %INPUT%\%FILENAME%.tdf %TMP%\%FILENAME%_corr.tdf %AIR_LEFT% %AIR_RIGHT% False False 0 %EXT_FOV% %EXT_FOV_RIGHT% %EXT_FOV_OVERLAP% %LINE_BY_LINE% %AVG_OVERLAP% %EXT_FOV_PROJ% %RR_METHOD% %DYN_FF% %THREADS% %LOG%\%FILENAME%_corr_log.txt
echo Phase retrieval...start_time:%time%
C:\STP_1.6\python\python-2.7.10.amd64\python.exe %CORE%\exec_phaseretrieval.py %START_PROJ% %END_PROJ% %TMP%\%FILENAME%_corr.tdf %TMP%\%FILENAME%_corr_phrt.tdf %PR_METHOD% %BETA% %DELTA% %ENERGY% %DISTANCE% %PIXEL% %OVERPAD% %THREADS% %LOG%\%FILENAME%_phrt_log.txt
del %TMP%\%FILENAME%_corr.tdf
echo Reconstruction...start_time:%time%
C:\STP_1.6\python\python-2.7.10.amd64\python.exe %CORE%\exec_reconstruct.py %START_PROJ% 2048 %TMP%\%FILENAME%_corr_phrt.tdf %RECON%\%FILENAME%\slices %ANG_RANGE% %OFFSET% %FILTER% %SCALE% %OVERPADDING% %LOGTRSF% %CIRCLE% slice False False False 0 False False 0 False False 0 0 0 False "rivers:3;0" False 0 %RECON_METHOD% %DECIM_FACTOR% %DOWNSC_FACTOR% %POSTPROC% %POLAR_FILT% %CONVERT_OPT% %CROP_OPT% %START_PROJ% %END_PROJ% %ROLLING% %ROLL_SHIFT% False %THREADS% %LOG%\%FILENAME%_recon_log_%ID%.txt
echo Reconstruction...finish_time:%time%
del %TMP%\%FILENAME%_corr_phrt.tdf
set /a x+=1
GOTO :MainLoop
:FailedLoop
echo "*** Offset not defined for file %FILENAME% ***"
GOTO :end

:endLoop
echo "Finished!"
echo Finished time: %time%
GOTO :end

:end
PAUSE