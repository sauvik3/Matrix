@echo off

set mytarget=%1
if /I "%1"=="" (
	set mytarget=clean
) else (
	if /I "%1"=="all" (
		set mytarget=allclean
	) else (
		echo I'm only a small script...
		echo.&pause&goto:eof
	)
)

echo Cleaning...

del /f /q .\obj\*.obj >nul 2>&1
del /f /q .\rc\*.res >nul 2>&1


if /I "%mytarget%"=="allclean" (
	del /f /q .\lib\matrix.lib >nul 2>&1
	del /f /q .\bin\matrix.lib >nul 2>&1
	del /f /q .\bin\matrix.exp >nul 2>&1
	del /f /q .\bin\matrix.dll >nul 2>&1
)
echo Complete...