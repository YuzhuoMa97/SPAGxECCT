@echo off
setlocal
cd /d "%~dp0..\.."
call examples\tutorials\01_spagxe.cmd
if errorlevel 1 exit /b 1
call examples\tutorials\02_spagxe_plus.cmd
if errorlevel 1 exit /b 1
call examples\tutorials\03_spagxemix.cmd
if errorlevel 1 exit /b 1
call examples\tutorials\04_spagxemix_local.cmd
if errorlevel 1 exit /b 1
call examples\tutorials\05_global_local_cct.cmd
if errorlevel 1 exit /b 1
call examples\tutorials\06_local_packing.cmd
if errorlevel 1 exit /b 1
exit /b 0
