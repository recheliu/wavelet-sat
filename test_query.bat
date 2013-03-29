@ECHO OFF
SETLOCAL ENABLEDELAYEDEXPANSION
SET ROOT_PATH=build\vs2010\x64
if "%1" == "F" (
	SET ROOT_PATH=!ROOT_PATH!\SATFile
	SET MS=8192 0
	SET ZS=1 0
	SET EXT=nc4
	SET WS=1
) ELSE (
	SET MS=8192
	SET ZS=1 0
	SET EXT=wnc4
	IF "%1" == "W" (
		SET WS=1 4 16
	) ELSE (
		SET WS=1
	)
)
SET BINS=32 64 96 128
SET DECODER=!ROOT_PATH!\SimpleNDQuery\Release\SimpleNDQuery.exe
for %%w in (!WS!) do (
	for %%z in (!ZS!) do (
		for %%m in (!MS!) do (
			SET ARGLIST=--size-of-full-arrays %%m --is-testing-query true --n-testing-values 8192 --query-win-length %%w
			for /F "eol=# tokens=1,2" %%i in (benchmark.txt) do (
				ECHO ====================================================
				ECHO %%j: !ARGLIST!
				for %%b in (!BINS!) do (
					SET TEST_ARGLIST=--nc-filepath D:\projects\WaveletSAT\build\vs2010\x64\output\%%j.b_%%b.z_%%z.!EXT! !ARGLIST!
					ECHO !TEST_ARGLIST!
					%DECODER% !TEST_ARGLIST!
				)
			)
			
			ECHO ====================================================
			ECHO Ocean: !ARGLIST!
			for %%b in (!BINS!) do (
				SET TEST_ARGLIST=--nc-filepath D:\projects\WaveletSAT\build\vs2010\x64\output\Ocean.20010201.b_%%b.z_%%z.!EXT! !ARGLIST!
				ECHO !TEST_ARGLIST!
				%DECODER% !TEST_ARGLIST!
			)

			ECHO ====================================================
			ECHO MJO: !ARGLIST!
			for %%b in (!BINS!) do (
				SET TEST_ARGLIST=--nc-filepath D:\projects\WaveletSAT\build\vs2010\x64\output\MJO.2008-01-01_00-00-00.QVAPOR.b_%%b.z_%%z.!EXT! !ARGLIST!
				ECHO !TEST_ARGLIST!
				%DECODER% !TEST_ARGLIST!
			)
		)
	)
)
ENDLOCAL