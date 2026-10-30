@echo off
setlocal DisableDelayedExpansion
rem Works at the bundle root or from cpp\examples in a source checkout.
set "CPP_DIR=%~dp0cpp"
if exist "%CPP_DIR%\examples\tutorials\run_all.cmd" goto cpp_located
set "CPP_DIR=%~dp0.."
if exist "%CPP_DIR%\examples\tutorials\run_all.cmd" goto cpp_located
echo Cannot locate cpp\examples\tutorials\run_all.cmd. Extract the complete bundle. 1>&2
exit /b 1
:cpp_located
for %%I in ("%CPP_DIR%") do set "CPP_DIR=%%~fI"
if not exist "%CPP_DIR%\SPAGxE.exe" goto source_checkout
set "SPAGXE_EXE=%CPP_DIR%\SPAGxE.exe"
goto executable_selected
:source_checkout
if defined SPAGXE_EXE goto executable_selected
set "SPAGXE_EXE=%CPP_DIR%\build\bin\SPAGxE.exe"
:executable_selected
for %%I in ("%SPAGXE_EXE%") do set "SPAGXE_EXE=%%~fI"
if exist "%SPAGXE_EXE%" goto run_examples
echo SPAGxE executable is missing: "%SPAGXE_EXE%" 1>&2
exit /b 1
:run_examples
call "%CPP_DIR%\examples\tutorials\run_all.cmd"
exit /b %errorlevel%
