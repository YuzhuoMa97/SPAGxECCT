@echo off
setlocal DisableDelayedExpansion
cd /d "%~dp0..\.."
if not defined SPAGXE_EXE set "SPAGXE_EXE=build\bin\SPAGxE.exe"
if not exist results\tutorials mkdir results\tutorials
rem C++ count packing; refuses existing outputs. Run in a fresh results folder.
"%SPAGXE_EXE%" --pack-local --local-manifest examples/local/binary_demo/cohort.local.tsv --local-bim examples/local/binary_demo/local.localbim --local-fam examples/local/binary_demo/local.localfam --out results/tutorials/packed_demo
if errorlevel 1 exit /b 1
rem Same published local method, reading the newly packed files.
"%SPAGXE_EXE%" --method SPAGxEmix_local --local-prefix results/tutorials/packed_demo --pheno examples/local/binary_demo/samples.pheno --resid-name R --envir-name E --mac 0 --maf 0 --geno 1 --out results/tutorials/local_packed
exit /b %errorlevel%
