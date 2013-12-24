@ECHO OFF
SETLOCAL ENABLEDELAYEDEXPANSION
SET ROOT_PATH=build\vs2010\x64
SET BINS=32 64 96 128
SET ZS=0 1
if "%1" == "F" (
	SET ROOT_PATH=!ROOT_PATH!\SATFile
	SET GS=false
	SET MS=0
) ELSE (
	if "%1" == "S" (
		SET GS=true false
		SET MS=8192 4096 2048 0
		SET BINS=128
		SET ZS=1
	) ELSE (	
		SET GS=true false
		SET MS=8192 0
	)
)

SET ENCODER=!ROOT_PATH!\SimpleNDEncoder\Release\SimpleNDEncoder.exe
SET OCEAN_ENCODER=!ROOT_PATH!\Ocean\Release\Ocean.exe
SET MJO_ENCODER=!ROOT_PATH!\MJO\Release\MJO.exe
for %%g in (!GS!) do (
	for %%z in (!ZS!) do (
		for %%m in (!MS!) do (
			SET ARGLIST=--netcdf-deflate-level %%z --is-using-gpus %%g --size-of-full-arrays %%m
		
			IF NOT "%1" == "S" (
				ECHO ====================================================
				SET i=simulated
				SET j=hydrogenAtom
				ECHO !j!: !ARGLIST!
				for %%b in (!BINS!) do (
					SET TEST_ARGLIST=--n-bins %%b --vol-filepath D:\data\volumes\rawfiles\!i!\!j!.nhdr   --nc-filepath-prefix D:\projects\WaveletSAT\build\vs2010\x64\output\!j!.b_%%b.z_%%z !ARGLIST!
					ECHO !TEST_ARGLIST!
					%ENCODER% !TEST_ARGLIST!
				)
			)

			ECHO ====================================================
			SET i=medical
			SET j=foot
			ECHO !j!: !ARGLIST!
			for %%b in (!BINS!) do (
				SET TEST_ARGLIST=--n-bins %%b --vol-filepath D:\data\volumes\rawfiles\!i!\!j!.nhdr   --nc-filepath-prefix D:\projects\WaveletSAT\build\vs2010\x64\output\!j!.b_%%b.z_%%z !ARGLIST!
				ECHO !TEST_ARGLIST!
				%ENCODER% !TEST_ARGLIST!
			)

			ECHO ====================================================
			ECHO Ocean: !ARGLIST!
			for %%b in (!BINS!) do (
				SET TEST_ARGLIST=--n-bins %%b --nc-filepath-prefix D:\projects\WaveletSAT\build\vs2010\x64\output\Ocean.20010201.b_%%b.z_%%z --vec-filepaths D:\data\esg_ocean UVEL.t.20c.20010201.nc UVEL VVEL.t.20c.20010201.nc VVEL --depth 4 --timing-printing-level 1 !ARGLIST!
				ECHO !TEST_ARGLIST!
				%OCEAN_ENCODER% !TEST_ARGLIST!
			)

			IF NOT "%1" == "S" (
				ECHO ====================================================
				ECHO MJO: !ARGLIST!
				for %%b in (!BINS!) do (
					SET TEST_ARGLIST=--n-bins %%b --nc-filepath-prefix D:\projects\WaveletSAT\build\vs2010\x64\output\MJO.2008-01-01_00-00-00.QVAPOR.b_%%b.z_%%z  --data D:\data\flow\datasets\smhagos\1 wrfout_d01_2008-01-01_00-00-00 QVAPOR --timing-printing-level 1 !ARGLIST!
					ECHO !TEST_ARGLIST!
					%MJO_ENCODER% !TEST_ARGLIST!
				)
			)
		)
	)
)
ENDLOCAL

REM 
REM Copyright (c) 2013 Teng-Yok Lee
REM
REM See the file LICENSE.txt for copying permission.
REM
