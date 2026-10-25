@echo off
rem run_flu.bat - launcher for flu_pipeline.sh inside WSL
rem   run_flu.bat setup
rem   run_flu.bat C:\Chim_nanopore\barcode21 -e your@email.com
setlocal
if "%~1"=="" goto usage
if /i "%~1"=="setup" goto setup
set "IN=%~1"
set "IN=%IN:\=/%"
shift
wsl --cd "%~dp0." -e bash flu_pipeline.sh -i "%IN%" %1 %2 %3 %4 %5 %6 %7 %8 %9
exit /b %errorlevel%
:setup
wsl --cd "%~dp0." -e bash flu_pipeline.sh -s
exit /b %errorlevel%
:usage
echo Usage: run_flu.bat FOLDER [options]   or   run_flu.bat setup
exit /b 1
