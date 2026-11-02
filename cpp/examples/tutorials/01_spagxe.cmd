@echo off
setlocal
cd /d "%~dp0..\.."
if not defined SPAGXE_EXE set "SPAGXE_EXE=build\bin\SPAGxE.exe"
if not exist results\tutorials mkdir results\tutorials
rem SPAGxE: unrelated subjects; paper weights with guarded partial-normal SPA.
rem Fit four traits from PLINK BED. Include E in the null-model covariates.
"%SPAGXE_EXE%" --method SPAGxE --bfile examples/synthetic/synthetic --pheno examples/synthetic/synthetic.pheno --pheno-name Quantitative,Binary,Time:Event,Ordinal --covar-name X1,E,PC1,PC2 --envir-name E --out results/tutorials/spagxe
if errorlevel 1 exit /b 1
exit /b 0
